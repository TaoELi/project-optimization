Iteration 0003 — 4a79c588ea08 (accepted)
========================================


GitHub commit: `4a79c588ea08 <iter-0003-page-head-4a79c588ea08_>`_
Published branch: `fermilink-optimize/lammps-bond <https://github.com/skilled-scipkg/lammps/tree/fermilink-optimize%2Flammps-bond>`_

Change summary
--------------


Add a `bond_class2` single-bond-type fast path that dispatches to a dedicated `eval_one_type` helper and hoists constant coefficients out of the inner loop.

Acceptance rationale
--------------------


Correctness passed and the single-bond-type `bond_class2` fast path improved weighted_median_bond_seconds from 0.28156 to 0.27301 (-3.04% vs incumbent).

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
| incumbent commit | `911f06a50100 <iter-0003-guardrails-incumbent-911f06a501008c67741eb0ff406d6554eccb81d7_>`_ |
+------------------+--------------------------------------------------------------------------------------------+
| candidate commit | `4a79c588ea08 <iter-0003-guardrails-candidate-4a79c588ea0869e5216d32c826251941d28332ac_>`_ |
+------------------+--------------------------------------------------------------------------------------------+
| incumbent metric | 0.28156                                                                                    |
+------------------+--------------------------------------------------------------------------------------------+
| candidate metric | 0.27301                                                                                    |
+------------------+--------------------------------------------------------------------------------------------+
| baseline metric  | 0.32627                                                                                    |
+------------------+--------------------------------------------------------------------------------------------+
| Δ vs incumbent   | +3.037% (lower-is-better sign)                                                             |
+------------------+--------------------------------------------------------------------------------------------+
| changed files    | src/CLASS2/bond_class2.cpp, src/CLASS2/bond_class2.h                                       |
+------------------+--------------------------------------------------------------------------------------------+


Diffstat
--------


.. code-block:: text

    src/CLASS2/bond_class2.cpp | 94 ++++++++++++++++++++++++++++++++++++++++++----
    src/CLASS2/bond_class2.h   |  1 +
    2 files changed, 87 insertions(+), 8 deletions(-)

Diff
----


:download:`download full diff <_diffs/iter_0003_4a79c588ea08.diff>`

.. code-block:: diff

   diff --git a/src/CLASS2/bond_class2.cpp b/src/CLASS2/bond_class2.cpp
   index f87f69fba7..4eb7c0fb39 100644
   --- a/src/CLASS2/bond_class2.cpp
   +++ b/src/CLASS2/bond_class2.cpp
   @@ -66,22 +66,38 @@ void BondClass2::compute(int eflag, int vflag)
    {
      int nbondlist = neighbor->nbondlist;
      int newton_bond = force->newton_bond;
   +  const int one_type = atom->nbondtypes == 1;
    
      ev_init(eflag,vflag);
    
      if (nbondlist <= 0) return;
    
   -  if (evflag) {
   -    if (eflag) {
   -      if (newton_bond) eval<1,1,1>(nbondlist);
   -      else eval<1,1,0>(nbondlist);
   +  if (one_type) {
   +    if (evflag) {
   +      if (eflag) {
   +        if (newton_bond) eval_one_type<1,1,1>(nbondlist);
   +        else eval_one_type<1,1,0>(nbondlist);
   +      } else {
   +        if (newton_bond) eval_one_type<1,0,1>(nbondlist);
   +        else eval_one_type<1,0,0>(nbondlist);
   +      }
        } else {
   -      if (newton_bond) eval<1,0,1>(nbondlist);
   -      else eval<1,0,0>(nbondlist);
   +      if (newton_bond) eval_one_type<0,0,1>(nbondlist);
   +      else eval_one_type<0,0,0>(nbondlist);
        }
      } else {
   -    if (newton_bond) eval<0,0,1>(nbondlist);
   -    else eval<0,0,0>(nbondlist);
   +    if (evflag) {
   +      if (eflag) {
   +        if (newton_bond) eval<1,1,1>(nbondlist);
   +        else eval<1,1,0>(nbondlist);
   +      } else {
   +        if (newton_bond) eval<1,0,1>(nbondlist);
   +        else eval<1,0,0>(nbondlist);
   +      }
   +    } else {
   +      if (newton_bond) eval<0,0,1>(nbondlist);
   +      else eval<0,0,0>(nbondlist);
   +    }
      }
    }
    
   @@ -145,6 +161,68 @@ void BondClass2::eval(int nbondlist)
      }
    }
    
   +template <int EVFLAG, int EFLAG, int NEWTON_BOND>
   +void BondClass2::eval_one_type(int nbondlist)
   +{
   +  int i1, i2;
   +  double delx, dely, delz, ebond, fbond;
   +  double rsq, r, dr, dr2, dr3, dr4, de_bond;
   +
   +  const auto *_noalias const x = (dbl3_t *) atom->x[0];
   +  auto *_noalias const f = (dbl3_t *) atom->f[0];    // NOLINT
   +  const int3_t *_noalias const bondlist = (int3_t *) neighbor->bondlist[0];
   +  const int nlocal = atom->nlocal;
   +  const double coeff_r0 = r0[1];
   +  const double coeff_k2 = k2[1];
   +  const double coeff_k3 = k3[1];
   +  const double coeff_k4 = k4[1];
   +  const double coeff_2k2 = 2.0 * coeff_k2;
   +  const double coeff_3k3 = 3.0 * coeff_k3;
   +  const double coeff_4k4 = 4.0 * coeff_k4;
   +
   +  ebond = 0.0;
   +
   +  for (int n = 0; n < nbondlist; n++) {
   +    i1 = bondlist[n].a;
   +    i2 = bondlist[n].b;
   +
   +    delx = x[i1].x - x[i2].x;
   +    dely = x[i1].y - x[i2].y;
   +    delz = x[i1].z - x[i2].z;
   +
   +    rsq = delx*delx + dely*dely + delz*delz;
   +    r = sqrt(rsq);
   +    dr = r - coeff_r0;
   +    dr2 = dr*dr;
   +    dr3 = dr2*dr;
   +    dr4 = dr3*dr;
   +
   +    // force & energy
   +
   +    de_bond = coeff_2k2*dr + coeff_3k3*dr2 + coeff_4k4*dr3;
   +    if (r > 0.0) fbond = -de_bond/r;
   +    else fbond = 0.0;
   +
   +    if (EFLAG) ebond = coeff_k2*dr2 + coeff_k3*dr3 + coeff_k4*dr4;
   +
   +    // apply force to each of 2 atoms
   +
   +    if (NEWTON_BOND || i1 < nlocal) {
   +      f[i1].x += delx*fbond;
   +      f[i1].y += dely*fbond;
   +      f[i1].z += delz*fbond;
   +    }
   +
   +    if (NEWTON_BOND || i2 < nlocal) {
   +      f[i2].x -= delx*fbond;
   +      f[i2].y -= dely*fbond;
   +      f[i2].z -= delz*fbond;
   +    }
   +
   +    if (EVFLAG) ev_tally(i1,i2,nlocal,NEWTON_BOND,ebond,fbond,delx,dely,delz);
   +  }
   +}
   +
    /* ---------------------------------------------------------------------- */
    
    void BondClass2::allocate()
   diff --git a/src/CLASS2/bond_class2.h b/src/CLASS2/bond_class2.h
   index 6c5ca116b1..de72ebe03b 100644
   --- a/src/CLASS2/bond_class2.h
   +++ b/src/CLASS2/bond_class2.h
   @@ -43,6 +43,7 @@ class BondClass2 : public Bond {
    
      virtual void allocate();
      template <int EVFLAG, int EFLAG, int NEWTON_BOND> void eval(int);
   +  template <int EVFLAG, int EFLAG, int NEWTON_BOND> void eval_one_type(int);
    };
    
    }    // namespace LAMMPS_NS


.. _iter-0003-page-head-4a79c588ea08: https://github.com/skilled-scipkg/lammps/commit/4a79c588ea08
.. _iter-0003-guardrails-incumbent-911f06a501008c67741eb0ff406d6554eccb81d7: https://github.com/skilled-scipkg/lammps/commit/911f06a501008c67741eb0ff406d6554eccb81d7
.. _iter-0003-guardrails-candidate-4a79c588ea0869e5216d32c826251941d28332ac: https://github.com/skilled-scipkg/lammps/commit/4a79c588ea0869e5216d32c826251941d28332ac