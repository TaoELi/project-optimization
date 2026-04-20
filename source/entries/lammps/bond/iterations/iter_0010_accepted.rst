Iteration 0010 — cf4ba8e7204d (accepted)
========================================


GitHub commit: `cf4ba8e7204d <iter-0010-page-head-cf4ba8e7204d_>`_
Published branch: `fermilink-optimize/lammps-bond <https://github.com/skilled-scipkg/lammps/tree/fermilink-optimize%2Flammps-bond>`_

Change summary
--------------


Combine the correctness-safe `bond_class2` one-type two-way unrolled hot kernel with the correctness-safe `angle_harmonic` one-type dispatch and constant hoisting on top of incumbent `4a79c588ea08`.

Acceptance rationale
--------------------


Correctness passed and weighted_median_bond_seconds improved from 0.27301 to 0.26621 (+2.49% vs incumbent), clearing the 2% promotion threshold.

Guardrails & metrics
--------------------


+------------------+----------------------------------------------------------------------------------------------------------------------+
| field            | value                                                                                                                |
+==================+======================================================================================================================+
| decision         | ACCEPTED                                                                                                             |
+------------------+----------------------------------------------------------------------------------------------------------------------+
| correctness      | ok                                                                                                                   |
+------------------+----------------------------------------------------------------------------------------------------------------------+
| correctness mode | field_tolerances                                                                                                     |
+------------------+----------------------------------------------------------------------------------------------------------------------+
| hard reject      | no                                                                                                                   |
+------------------+----------------------------------------------------------------------------------------------------------------------+
| guardrail errors | 0                                                                                                                    |
+------------------+----------------------------------------------------------------------------------------------------------------------+
| incumbent commit | `4a79c588ea08 <iter-0010-guardrails-incumbent-4a79c588ea0869e5216d32c826251941d28332ac_>`_                           |
+------------------+----------------------------------------------------------------------------------------------------------------------+
| candidate commit | `cf4ba8e7204d <iter-0010-guardrails-candidate-cf4ba8e7204d7484a30eb354b07f05005552dd43_>`_                           |
+------------------+----------------------------------------------------------------------------------------------------------------------+
| incumbent metric | 0.27301                                                                                                              |
+------------------+----------------------------------------------------------------------------------------------------------------------+
| candidate metric | 0.26621                                                                                                              |
+------------------+----------------------------------------------------------------------------------------------------------------------+
| baseline metric  | 0.32627                                                                                                              |
+------------------+----------------------------------------------------------------------------------------------------------------------+
| Δ vs incumbent   | +2.491% (lower-is-better sign)                                                                                       |
+------------------+----------------------------------------------------------------------------------------------------------------------+
| changed files    | src/CLASS2/bond_class2.cpp, src/CLASS2/bond_class2.h, src/MOLECULE/angle_harmonic.cpp, src/MOLECULE/angle_harmonic.h |
+------------------+----------------------------------------------------------------------------------------------------------------------+


Diffstat
--------


.. code-block:: text

    src/CLASS2/bond_class2.cpp      | 102 +++++++++++++++++++++++++++++++
    src/CLASS2/bond_class2.h        |   1 +
    src/MOLECULE/angle_harmonic.cpp | 130 +++++++++++++++++++++++++++++++++++++---
    src/MOLECULE/angle_harmonic.h   |   1 +
    4 files changed, 226 insertions(+), 8 deletions(-)

Diff
----


:download:`download full diff <_diffs/iter_0010_cf4ba8e7204d.diff>`

.. code-block:: diff

   diff --git a/src/CLASS2/bond_class2.cpp b/src/CLASS2/bond_class2.cpp
   index 4eb7c0fb39..0984b0a7e1 100644
   --- a/src/CLASS2/bond_class2.cpp
   +++ b/src/CLASS2/bond_class2.cpp
   @@ -72,6 +72,11 @@ void BondClass2::compute(int eflag, int vflag)
    
      if (nbondlist <= 0) return;
    
   +  if (one_type && !evflag && newton_bond) {
   +    eval_one_type_hot(nbondlist);
   +    return;
   +  }
   +
      if (one_type) {
        if (evflag) {
          if (eflag) {
   @@ -101,6 +106,103 @@ void BondClass2::compute(int eflag, int vflag)
      }
    }
    
   +void BondClass2::eval_one_type_hot(int nbondlist)
   +{
   +  const auto *_noalias const x = (dbl3_t *) atom->x[0];
   +  auto *_noalias const f = (dbl3_t *) atom->f[0];    // NOLINT
   +  const int3_t *_noalias const bondlist = (int3_t *) neighbor->bondlist[0];
   +  const double coeff_r0 = r0[1];
   +  const double coeff_2k2 = 2.0 * k2[1];
   +  const double coeff_3k3 = 3.0 * k3[1];
   +  const double coeff_4k4 = 4.0 * k4[1];
   +
   +  int n = 0;
   +
   +  // Two-way unrolling gives the compiler more independent exact work around the sqrt/div pair.
   +  for (; n + 1 < nbondlist; n += 2) {
   +    const int i1_0 = bondlist[n].a;
   +    const int i2_0 = bondlist[n].b;
   +    const int i1_1 = bondlist[n + 1].a;
   +    const int i2_1 = bondlist[n + 1].b;
   +
   +    const double delx0 = x[i1_0].x - x[i2_0].x;
   +    const double dely0 = x[i1_0].y - x[i2_0].y;
   +    const double delz0 = x[i1_0].z - x[i2_0].z;
   +    const double delx1 = x[i1_1].x - x[i2_1].x;
   +    const double dely1 = x[i1_1].y - x[i2_1].y;
   +    const double delz1 = x[i1_1].z - x[i2_1].z;
   +
   +    const double rsq0 = delx0 * delx0 + dely0 * dely0 + delz0 * delz0;
   +    const double rsq1 = delx1 * delx1 + dely1 * dely1 + delz1 * delz1;
   +    const double r0 = sqrt(rsq0);
   +    const double r1 = sqrt(rsq1);
   +    const double dr0 = r0 - coeff_r0;
   +    const double dr1 = r1 - coeff_r0;
   +    const double dr20 = dr0 * dr0;
   +    const double dr21 = dr1 * dr1;
   +    const double dr30 = dr20 * dr0;
   +    const double dr31 = dr21 * dr1;
   +
   +    const double de_bond0 = coeff_2k2 * dr0 + coeff_3k3 * dr20 + coeff_4k4 * dr30;
   +    const double de_bond1 = coeff_2k2 * dr1 + coeff_3k3 * dr21 + coeff_4k4 * dr31;
   +    double fbond0, fbond1;
   +    if (r0 > 0.0) fbond0 = -de_bond0 / r0;
   +    else fbond0 = 0.0;
   +    if (r1 > 0.0) fbond1 = -de_bond1 / r1;
   +    else fbond1 = 0.0;
   +
   +    const double fx0 = delx0 * fbond0;
   +    const double fy0 = dely0 * fbond0;
   +    const double fz0 = delz0 * fbond0;
   +    f[i1_0].x += fx0;
   +    f[i1_0].y += fy0;
   +    f[i1_0].z += fz0;
   +    f[i2_0].x -= fx0;
   +    f[i2_0].y -= fy0;
   +    f[i2_0].z -= fz0;
   +
   +    const double fx1 = delx1 * fbond1;
   +    const double fy1 = dely1 * fbond1;
   +    const double fz1 = delz1 * fbond1;
   +    f[i1_1].x += fx1;
   +    f[i1_1].y += fy1;
   +    f[i1_1].z += fz1;
   +    f[i2_1].x -= fx1;
   +    f[i2_1].y -= fy1;
   +    f[i2_1].z -= fz1;
   +  }
   +
   +  for (; n < nbondlist; n++) {
   +    const int i1 = bondlist[n].a;
   +    const int i2 = bondlist[n].b;
   +
   +    const double delx = x[i1].x - x[i2].x;
   +    const double dely = x[i1].y - x[i2].y;
   +    const double delz = x[i1].z - x[i2].z;
   +
   +    const double rsq = delx * delx + dely * dely + delz * delz;
   +    const double r = sqrt(rsq);
   +    const double dr = r - coeff_r0;
   +    const double dr2 = dr * dr;
   +    const double dr3 = dr2 * dr;
   +
   +    const double de_bond = coeff_2k2 * dr + coeff_3k3 * dr2 + coeff_4k4 * dr3;
   +    double fbond;
   +    if (r > 0.0) fbond = -de_bond / r;
   +    else fbond = 0.0;
   +
   +    const double fx = delx * fbond;
   +    const double fy = dely * fbond;
   +    const double fz = delz * fbond;
   +    f[i1].x += fx;
   +    f[i1].y += fy;
   +    f[i1].z += fz;
   +    f[i2].x -= fx;
   +    f[i2].y -= fy;
   +    f[i2].z -= fz;
   +  }
   +}
   +
    template <int EVFLAG, int EFLAG, int NEWTON_BOND>
    void BondClass2::eval(int nbondlist)
    {
   diff --git a/src/CLASS2/bond_class2.h b/src/CLASS2/bond_class2.h
   index de72ebe03b..d8ae006de7 100644
   --- a/src/CLASS2/bond_class2.h
   +++ b/src/CLASS2/bond_class2.h
   @@ -42,6 +42,7 @@ class BondClass2 : public Bond {
      double *r0, *k2, *k3, *k4;
    
      virtual void allocate();
   +  void eval_one_type_hot(int);
      template <int EVFLAG, int EFLAG, int NEWTON_BOND> void eval(int);
      template <int EVFLAG, int EFLAG, int NEWTON_BOND> void eval_one_type(int);
    };
   diff --git a/src/MOLECULE/angle_harmonic.cpp b/src/MOLECULE/angle_harmonic.cpp
   index 49f2849985..1c2eec9828 100644
   --- a/src/MOLECULE/angle_harmonic.cpp
   +++ b/src/MOLECULE/angle_harmonic.cpp
   @@ -67,20 +67,36 @@ void AngleHarmonic::compute(int eflag, int vflag)
    
      int nanglelist = neighbor->nanglelist;
      int newton_bond = force->newton_bond;
   +  const int one_type = atom->nangletypes == 1;
    
      if (nanglelist <= 0) return;
    
   -  if (evflag) {
   -    if (eflag) {
   -      if (newton_bond) eval<1,1,1>(nanglelist);
   -      else eval<1,1,0>(nanglelist);
   +  if (one_type) {
   +    if (evflag) {
   +      if (eflag) {
   +        if (newton_bond) eval_one_type<1,1,1>(nanglelist);
   +        else eval_one_type<1,1,0>(nanglelist);
   +      } else {
   +        if (newton_bond) eval_one_type<1,0,1>(nanglelist);
   +        else eval_one_type<1,0,0>(nanglelist);
   +      }
        } else {
   -      if (newton_bond) eval<1,0,1>(nanglelist);
   -      else eval<1,0,0>(nanglelist);
   +      if (newton_bond) eval_one_type<0,0,1>(nanglelist);
   +      else eval_one_type<0,0,0>(nanglelist);
        }
      } else {
   -    if (newton_bond) eval<0,0,1>(nanglelist);
   -    else eval<0,0,0>(nanglelist);
   +    if (evflag) {
   +      if (eflag) {
   +        if (newton_bond) eval<1,1,1>(nanglelist);
   +        else eval<1,1,0>(nanglelist);
   +      } else {
   +        if (newton_bond) eval<1,0,1>(nanglelist);
   +        else eval<1,0,0>(nanglelist);
   +      }
   +    } else {
   +      if (newton_bond) eval<0,0,1>(nanglelist);
   +      else eval<0,0,0>(nanglelist);
   +    }
      }
    }
    
   @@ -183,6 +199,104 @@ void AngleHarmonic::eval(int nanglelist)
      }
    }
    
   +template <int EVFLAG, int EFLAG, int NEWTON_BOND>
   +void AngleHarmonic::eval_one_type(int nanglelist)
   +{
   +  int i1, i2, i3;
   +  double delx1, dely1, delz1, delx2, dely2, delz2;
   +  double eangle, f1[3], f3[3];
   +  double dtheta, tk;
   +  double rsq1, rsq2, r1, r2, c, s, a, a11, a12, a22;
   +
   +  const auto *_noalias const x = (dbl3_t *) atom->x[0];
   +  auto *_noalias const f = (dbl3_t *) atom->f[0];    // NOLINT
   +  const int4_t *_noalias const anglelist = (int4_t *) neighbor->anglelist[0];
   +  const int nlocal = atom->nlocal;
   +  const double coeff_k = k[1];
   +  const double coeff_theta0 = theta0[1];
   +
   +  eangle = 0.0;
   +
   +  for (int n = 0; n < nanglelist; n++) {
   +    i1 = anglelist[n].a;
   +    i2 = anglelist[n].b;
   +    i3 = anglelist[n].c;
   +
   +    // 1st bond
   +
   +    delx1 = x[i1].x - x[i2].x;
   +    dely1 = x[i1].y - x[i2].y;
   +    delz1 = x[i1].z - x[i2].z;
   +
   +    rsq1 = delx1 * delx1 + dely1 * dely1 + delz1 * delz1;
   +    r1 = sqrt(rsq1);
   +
   +    // 2nd bond
   +
   +    delx2 = x[i3].x - x[i2].x;
   +    dely2 = x[i3].y - x[i2].y;
   +    delz2 = x[i3].z - x[i2].z;
   +
   +    rsq2 = delx2 * delx2 + dely2 * dely2 + delz2 * delz2;
   +    r2 = sqrt(rsq2);
   +
   +    // angle (cos and sin)
   +
   +    c = delx1 * delx2 + dely1 * dely2 + delz1 * delz2;
   +    c /= r1 * r2;
   +
   +    if (c > 1.0) c = 1.0;
   +    if (c < -1.0) c = -1.0;
   +
   +    s = sqrt(1.0 - c * c);
   +    if (s < SMALL) s = SMALL;
   +    s = 1.0 / s;
   +
   +    // force & energy
   +
   +    dtheta = acos(c) - coeff_theta0;
   +    tk = coeff_k * dtheta;
   +
   +    if (EFLAG) eangle = tk * dtheta;
   +
   +    a = -2.0 * tk * s;
   +    a11 = a * c / rsq1;
   +    a12 = -a / (r1 * r2);
   +    a22 = a * c / rsq2;
   +
   +    f1[0] = a11 * delx1 + a12 * delx2;
   +    f1[1] = a11 * dely1 + a12 * dely2;
   +    f1[2] = a11 * delz1 + a12 * delz2;
   +    f3[0] = a22 * delx2 + a12 * delx1;
   +    f3[1] = a22 * dely2 + a12 * dely1;
   +    f3[2] = a22 * delz2 + a12 * delz1;
   +
   +    // apply force to each of 3 atoms
   +
   +    if (NEWTON_BOND || i1 < nlocal) {
   +      f[i1].x += f1[0];
   +      f[i1].y += f1[1];
   +      f[i1].z += f1[2];
   +    }
   +
   +    if (NEWTON_BOND || i2 < nlocal) {
   +      f[i2].x -= f1[0] + f3[0];
   +      f[i2].y -= f1[1] + f3[1];
   +      f[i2].z -= f1[2] + f3[2];
   +    }
   +
   +    if (NEWTON_BOND || i3 < nlocal) {
   +      f[i3].x += f3[0];
   +      f[i3].y += f3[1];
   +      f[i3].z += f3[2];
   +    }
   +
   +    if (EVFLAG)
   +      ev_tally(i1, i2, i3, nlocal, NEWTON_BOND, eangle, f1, f3, delx1, dely1, delz1, delx2, dely2,
   +               delz2);
   +  }
   +}
   +
    /* ---------------------------------------------------------------------- */
    
    void AngleHarmonic::allocate()
   diff --git a/src/MOLECULE/angle_harmonic.h b/src/MOLECULE/angle_harmonic.h
   index 0c58059539..104b1913ae 100644
   --- a/src/MOLECULE/angle_harmonic.h
   +++ b/src/MOLECULE/angle_harmonic.h
   @@ -43,6 +43,7 @@ class AngleHarmonic : public Angle {
    
      virtual void allocate();
      template <int EVFLAG, int EFLAG, int NEWTON_BOND> void eval(int);
   +  template <int EVFLAG, int EFLAG, int NEWTON_BOND> void eval_one_type(int);
    };
    
    }    // namespace LAMMPS_NS


.. _iter-0010-page-head-cf4ba8e7204d: https://github.com/skilled-scipkg/lammps/commit/cf4ba8e7204d
.. _iter-0010-guardrails-incumbent-4a79c588ea0869e5216d32c826251941d28332ac: https://github.com/skilled-scipkg/lammps/commit/4a79c588ea0869e5216d32c826251941d28332ac
.. _iter-0010-guardrails-candidate-cf4ba8e7204d7484a30eb354b07f05005552dd43: https://github.com/skilled-scipkg/lammps/commit/cf4ba8e7204d7484a30eb354b07f05005552dd43