Iteration 0011 — 672cbe901cae (accepted)
========================================


GitHub commit: `672cbe901cae <iter-0011-page-head-672cbe901cae_>`_
Published branch: `fermilink-optimize/lammps-bond <https://github.com/skilled-scipkg/lammps/tree/fermilink-optimize%2Flammps-bond>`_

Change summary
--------------


Add a dedicated `angle_harmonic` one-type no-`evflag`, Newton-on hot path with scalar force temporaries and two-way exact unrolling.

Acceptance rationale
--------------------


Correctness passed with no hard rejects, and the dedicated `angle_harmonic` hot path improved `weighted_median_bond_seconds` from 0.26621 to 0.24606 (+7.57% vs incumbent), clearing the 2% promotion threshold.

Guardrails & metrics
--------------------


+------------------+--------------------------------------------------------------------------------------------+
| field            | value                                                                                      |
+==================+============================================================================================+
| decision         | ACCEPTED                                                                                   |
+------------------+--------------------------------------------------------------------------------------------+
| correctness      | ok                                                                                         |
+------------------+--------------------------------------------------------------------------------------------+
| correctness mode | field_tolerances                                                                           |
+------------------+--------------------------------------------------------------------------------------------+
| hard reject      | no                                                                                         |
+------------------+--------------------------------------------------------------------------------------------+
| guardrail errors | 0                                                                                          |
+------------------+--------------------------------------------------------------------------------------------+
| incumbent commit | `cf4ba8e7204d <iter-0011-guardrails-incumbent-cf4ba8e7204d7484a30eb354b07f05005552dd43_>`_ |
+------------------+--------------------------------------------------------------------------------------------+
| candidate commit | `672cbe901cae <iter-0011-guardrails-candidate-672cbe901cae5b8e7041fb6f4c3ab883fd8e6576_>`_ |
+------------------+--------------------------------------------------------------------------------------------+
| incumbent metric | 0.26621                                                                                    |
+------------------+--------------------------------------------------------------------------------------------+
| candidate metric | 0.24606                                                                                    |
+------------------+--------------------------------------------------------------------------------------------+
| baseline metric  | 0.32627                                                                                    |
+------------------+--------------------------------------------------------------------------------------------+
| Δ vs incumbent   | +7.569% (lower-is-better sign)                                                             |
+------------------+--------------------------------------------------------------------------------------------+
| changed files    | src/MOLECULE/angle_harmonic.cpp, src/MOLECULE/angle_harmonic.h                             |
+------------------+--------------------------------------------------------------------------------------------+


Diffstat
--------


.. code-block:: text

    src/MOLECULE/angle_harmonic.cpp | 157 ++++++++++++++++++++++++++++++++++++++++
    src/MOLECULE/angle_harmonic.h   |   1 +
    2 files changed, 158 insertions(+)

Diff
----


:download:`download full diff <_diffs/iter_0011_672cbe901cae.diff>`

.. code-block:: diff

   diff --git a/src/MOLECULE/angle_harmonic.cpp b/src/MOLECULE/angle_harmonic.cpp
   index 1c2eec9828..44458e2813 100644
   --- a/src/MOLECULE/angle_harmonic.cpp
   +++ b/src/MOLECULE/angle_harmonic.cpp
   @@ -71,6 +71,11 @@ void AngleHarmonic::compute(int eflag, int vflag)
    
      if (nanglelist <= 0) return;
    
   +  if (one_type && !evflag && newton_bond) {
   +    eval_one_type_hot(nanglelist);
   +    return;
   +  }
   +
      if (one_type) {
        if (evflag) {
          if (eflag) {
   @@ -100,6 +105,158 @@ void AngleHarmonic::compute(int eflag, int vflag)
      }
    }
    
   +void AngleHarmonic::eval_one_type_hot(int nanglelist)
   +{
   +  const auto *_noalias const x = (dbl3_t *) atom->x[0];
   +  auto *_noalias const f = (dbl3_t *) atom->f[0];    // NOLINT
   +  const int4_t *_noalias const anglelist = (int4_t *) neighbor->anglelist[0];
   +  const double coeff_k = k[1];
   +  const double coeff_theta0 = theta0[1];
   +
   +  int n = 0;
   +
   +  // Two-way unrolling exposes more independent exact work around the sqrt/acos latency.
   +  for (; n + 1 < nanglelist; n += 2) {
   +    const int i1_0 = anglelist[n].a;
   +    const int i2_0 = anglelist[n].b;
   +    const int i3_0 = anglelist[n].c;
   +    const int i1_1 = anglelist[n + 1].a;
   +    const int i2_1 = anglelist[n + 1].b;
   +    const int i3_1 = anglelist[n + 1].c;
   +
   +    const double delx1_0 = x[i1_0].x - x[i2_0].x;
   +    const double dely1_0 = x[i1_0].y - x[i2_0].y;
   +    const double delz1_0 = x[i1_0].z - x[i2_0].z;
   +    const double delx2_0 = x[i3_0].x - x[i2_0].x;
   +    const double dely2_0 = x[i3_0].y - x[i2_0].y;
   +    const double delz2_0 = x[i3_0].z - x[i2_0].z;
   +    const double delx1_1 = x[i1_1].x - x[i2_1].x;
   +    const double dely1_1 = x[i1_1].y - x[i2_1].y;
   +    const double delz1_1 = x[i1_1].z - x[i2_1].z;
   +    const double delx2_1 = x[i3_1].x - x[i2_1].x;
   +    const double dely2_1 = x[i3_1].y - x[i2_1].y;
   +    const double delz2_1 = x[i3_1].z - x[i2_1].z;
   +
   +    const double rsq1_0 = delx1_0 * delx1_0 + dely1_0 * dely1_0 + delz1_0 * delz1_0;
   +    const double rsq2_0 = delx2_0 * delx2_0 + dely2_0 * dely2_0 + delz2_0 * delz2_0;
   +    const double rsq1_1 = delx1_1 * delx1_1 + dely1_1 * dely1_1 + delz1_1 * delz1_1;
   +    const double rsq2_1 = delx2_1 * delx2_1 + dely2_1 * dely2_1 + delz2_1 * delz2_1;
   +    const double r1_0 = sqrt(rsq1_0);
   +    const double r2_0 = sqrt(rsq2_0);
   +    const double r1_1 = sqrt(rsq1_1);
   +    const double r2_1 = sqrt(rsq2_1);
   +
   +    double c0 = delx1_0 * delx2_0 + dely1_0 * dely2_0 + delz1_0 * delz2_0;
   +    c0 /= r1_0 * r2_0;
   +    if (c0 > 1.0) c0 = 1.0;
   +    if (c0 < -1.0) c0 = -1.0;
   +    double s0 = sqrt(1.0 - c0 * c0);
   +    if (s0 < SMALL) s0 = SMALL;
   +    s0 = 1.0 / s0;
   +    const double dtheta0 = acos(c0) - coeff_theta0;
   +    const double tk0 = coeff_k * dtheta0;
   +    const double a0 = -2.0 * tk0 * s0;
   +    const double a11_0 = a0 * c0 / rsq1_0;
   +    const double a12_0 = -a0 / (r1_0 * r2_0);
   +    const double a22_0 = a0 * c0 / rsq2_0;
   +    const double f1x0 = a11_0 * delx1_0 + a12_0 * delx2_0;
   +    const double f1y0 = a11_0 * dely1_0 + a12_0 * dely2_0;
   +    const double f1z0 = a11_0 * delz1_0 + a12_0 * delz2_0;
   +    const double f3x0 = a22_0 * delx2_0 + a12_0 * delx1_0;
   +    const double f3y0 = a22_0 * dely2_0 + a12_0 * dely1_0;
   +    const double f3z0 = a22_0 * delz2_0 + a12_0 * delz1_0;
   +    f[i1_0].x += f1x0;
   +    f[i1_0].y += f1y0;
   +    f[i1_0].z += f1z0;
   +    f[i2_0].x -= f1x0 + f3x0;
   +    f[i2_0].y -= f1y0 + f3y0;
   +    f[i2_0].z -= f1z0 + f3z0;
   +    f[i3_0].x += f3x0;
   +    f[i3_0].y += f3y0;
   +    f[i3_0].z += f3z0;
   +
   +    double c1 = delx1_1 * delx2_1 + dely1_1 * dely2_1 + delz1_1 * delz2_1;
   +    c1 /= r1_1 * r2_1;
   +    if (c1 > 1.0) c1 = 1.0;
   +    if (c1 < -1.0) c1 = -1.0;
   +    double s1 = sqrt(1.0 - c1 * c1);
   +    if (s1 < SMALL) s1 = SMALL;
   +    s1 = 1.0 / s1;
   +    const double dtheta1 = acos(c1) - coeff_theta0;
   +    const double tk1 = coeff_k * dtheta1;
   +    const double a1 = -2.0 * tk1 * s1;
   +    const double a11_1 = a1 * c1 / rsq1_1;
   +    const double a12_1 = -a1 / (r1_1 * r2_1);
   +    const double a22_1 = a1 * c1 / rsq2_1;
   +    const double f1x1 = a11_1 * delx1_1 + a12_1 * delx2_1;
   +    const double f1y1 = a11_1 * dely1_1 + a12_1 * dely2_1;
   +    const double f1z1 = a11_1 * delz1_1 + a12_1 * delz2_1;
   +    const double f3x1 = a22_1 * delx2_1 + a12_1 * delx1_1;
   +    const double f3y1 = a22_1 * dely2_1 + a12_1 * dely1_1;
   +    const double f3z1 = a22_1 * delz2_1 + a12_1 * delz1_1;
   +    f[i1_1].x += f1x1;
   +    f[i1_1].y += f1y1;
   +    f[i1_1].z += f1z1;
   +    f[i2_1].x -= f1x1 + f3x1;
   +    f[i2_1].y -= f1y1 + f3y1;
   +    f[i2_1].z -= f1z1 + f3z1;
   +    f[i3_1].x += f3x1;
   +    f[i3_1].y += f3y1;
   +    f[i3_1].z += f3z1;
   +  }
   +
   +  for (; n < nanglelist; n++) {
   +    const int i1 = anglelist[n].a;
   +    const int i2 = anglelist[n].b;
   +    const int i3 = anglelist[n].c;
   +
   +    const double delx1 = x[i1].x - x[i2].x;
   +    const double dely1 = x[i1].y - x[i2].y;
   +    const double delz1 = x[i1].z - x[i2].z;
   +    const double delx2 = x[i3].x - x[i2].x;
   +    const double dely2 = x[i3].y - x[i2].y;
   +    const double delz2 = x[i3].z - x[i2].z;
   +
   +    const double rsq1 = delx1 * delx1 + dely1 * dely1 + delz1 * delz1;
   +    const double rsq2 = delx2 * delx2 + dely2 * dely2 + delz2 * delz2;
   +    const double r1 = sqrt(rsq1);
   +    const double r2 = sqrt(rsq2);
   +
   +    double c = delx1 * delx2 + dely1 * dely2 + delz1 * delz2;
   +    c /= r1 * r2;
   +    if (c > 1.0) c = 1.0;
   +    if (c < -1.0) c = -1.0;
   +
   +    double s = sqrt(1.0 - c * c);
   +    if (s < SMALL) s = SMALL;
   +    s = 1.0 / s;
   +
   +    const double dtheta = acos(c) - coeff_theta0;
   +    const double tk = coeff_k * dtheta;
   +    const double a = -2.0 * tk * s;
   +    const double a11 = a * c / rsq1;
   +    const double a12 = -a / (r1 * r2);
   +    const double a22 = a * c / rsq2;
   +
   +    const double f1x = a11 * delx1 + a12 * delx2;
   +    const double f1y = a11 * dely1 + a12 * dely2;
   +    const double f1z = a11 * delz1 + a12 * delz2;
   +    const double f3x = a22 * delx2 + a12 * delx1;
   +    const double f3y = a22 * dely2 + a12 * dely1;
   +    const double f3z = a22 * delz2 + a12 * delz1;
   +
   +    f[i1].x += f1x;
   +    f[i1].y += f1y;
   +    f[i1].z += f1z;
   +    f[i2].x -= f1x + f3x;
   +    f[i2].y -= f1y + f3y;
   +    f[i2].z -= f1z + f3z;
   +    f[i3].x += f3x;
   +    f[i3].y += f3y;
   +    f[i3].z += f3z;
   +  }
   +}
   +
    template <int EVFLAG, int EFLAG, int NEWTON_BOND>
    void AngleHarmonic::eval(int nanglelist)
    {
   diff --git a/src/MOLECULE/angle_harmonic.h b/src/MOLECULE/angle_harmonic.h
   index 104b1913ae..7e0e1a95af 100644
   --- a/src/MOLECULE/angle_harmonic.h
   +++ b/src/MOLECULE/angle_harmonic.h
   @@ -42,6 +42,7 @@ class AngleHarmonic : public Angle {
      double *k, *theta0;
    
      virtual void allocate();
   +  void eval_one_type_hot(int);
      template <int EVFLAG, int EFLAG, int NEWTON_BOND> void eval(int);
      template <int EVFLAG, int EFLAG, int NEWTON_BOND> void eval_one_type(int);
    };


.. _iter-0011-page-head-672cbe901cae: https://github.com/skilled-scipkg/lammps/commit/672cbe901cae
.. _iter-0011-guardrails-incumbent-cf4ba8e7204d7484a30eb354b07f05005552dd43: https://github.com/skilled-scipkg/lammps/commit/cf4ba8e7204d7484a30eb354b07f05005552dd43
.. _iter-0011-guardrails-candidate-672cbe901cae5b8e7041fb6f4c3ab883fd8e6576: https://github.com/skilled-scipkg/lammps/commit/672cbe901cae5b8e7041fb6f4c3ab883fd8e6576