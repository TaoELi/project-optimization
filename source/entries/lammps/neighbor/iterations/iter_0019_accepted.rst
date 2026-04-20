Iteration 0019 — 11db74893452 (accepted)
========================================


GitHub commit: `11db74893452 <iter-0019-page-head-11db74893452_>`_
Published branch: `fermilink-optimize/lammps-neighbor <https://github.com/skilled-scipkg/lammps/tree/fermilink-optimize%2Flammps-neighbor>`_

Change summary
--------------


Optimize Neighbor movement tracking by scanning contiguous x/xhold arrays in check_distance() and bulk-copying xhold in build() while keeping the incumbent npair_bin fast path unchanged

Acceptance rationale
--------------------


Correctness held and weighted_median_neigh_seconds improved from 0.29094 to 0.28487 (-2.09% vs incumbent), clearing the 2.0% promotion gate.

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
| incumbent commit | `08f5ab4a4e02 <iter-0019-guardrails-incumbent-08f5ab4a4e026b06a3d38c7e3ebd6ecca2a25635_>`_ |
+------------------+--------------------------------------------------------------------------------------------+
| candidate commit | `11db74893452 <iter-0019-guardrails-candidate-11db748934527306a5e9d2d8c0c51eb47b03158a_>`_ |
+------------------+--------------------------------------------------------------------------------------------+
| incumbent metric | 0.29094                                                                                    |
+------------------+--------------------------------------------------------------------------------------------+
| candidate metric | 0.28487                                                                                    |
+------------------+--------------------------------------------------------------------------------------------+
| baseline metric  | 0.36158                                                                                    |
+------------------+--------------------------------------------------------------------------------------------+
| Δ vs incumbent   | +2.086% (lower-is-better sign)                                                             |
+------------------+--------------------------------------------------------------------------------------------+
| changed files    | src/neighbor.cpp                                                                           |
+------------------+--------------------------------------------------------------------------------------------+


Diffstat
--------


.. code-block:: text

    src/neighbor.cpp | 45 +++++++++++++++++++++++++++++++++++++--------
    1 file changed, 37 insertions(+), 8 deletions(-)

Diff
----


:download:`download full diff <_diffs/iter_0019_11db74893452.diff>`

.. code-block:: diff

   diff --git a/src/neighbor.cpp b/src/neighbor.cpp
   index 6496a8c39a..391afb5899 100644
   --- a/src/neighbor.cpp
   +++ b/src/neighbor.cpp
   @@ -2424,12 +2424,42 @@ int Neighbor::check_distance()
        }
      } else deltasq = triggersq;
    
   -  double **x = atom->x;
   +  using dbl3_t = double[3];
   +
   +  const auto * _noalias const x = (const dbl3_t *) atom->x[0];
   +  const auto * _noalias const xhold = (const dbl3_t *) this->xhold[0];
      int nlocal = atom->nlocal;
      if (includegroup) nlocal = atom->nfirst;
    
      int flag = 0;
   -  for (int i = 0; i < nlocal; i++) {
   +  int i = 0;
   +  for (; i + 3 < nlocal; i += 4) {
   +    const double delx0 = x[i][0] - xhold[i][0];
   +    const double dely0 = x[i][1] - xhold[i][1];
   +    const double delz0 = x[i][2] - xhold[i][2];
   +    const double rsq0 = delx0*delx0 + dely0*dely0 + delz0*delz0;
   +
   +    const double delx1 = x[i+1][0] - xhold[i+1][0];
   +    const double dely1 = x[i+1][1] - xhold[i+1][1];
   +    const double delz1 = x[i+1][2] - xhold[i+1][2];
   +    const double rsq1 = delx1*delx1 + dely1*dely1 + delz1*delz1;
   +
   +    const double delx2 = x[i+2][0] - xhold[i+2][0];
   +    const double dely2 = x[i+2][1] - xhold[i+2][1];
   +    const double delz2 = x[i+2][2] - xhold[i+2][2];
   +    const double rsq2 = delx2*delx2 + dely2*dely2 + delz2*delz2;
   +
   +    const double delx3 = x[i+3][0] - xhold[i+3][0];
   +    const double dely3 = x[i+3][1] - xhold[i+3][1];
   +    const double delz3 = x[i+3][2] - xhold[i+3][2];
   +    const double rsq3 = delx3*delx3 + dely3*dely3 + delz3*delz3;
   +
   +    if ((rsq0 > deltasq) || (rsq1 > deltasq) || (rsq2 > deltasq) || (rsq3 > deltasq)) {
   +      flag = 1;
   +      break;
   +    }
   +  }
   +  for (; i < nlocal; i++) {
        delx = x[i][0] - xhold[i][0];
        dely = x[i][1] - xhold[i][1];
        delz = x[i][2] - xhold[i][2];
   @@ -2472,18 +2502,17 @@ void Neighbor::build(int topoflag)
      // store current atom positions and box size if needed
    
      if (dist_check) {
   -    double **x = atom->x;
   +    using dbl3_t = double[3];
   +
        if (includegroup) nlocal = atom->nfirst;
        if (atom->nmax > maxhold) {
          maxhold = atom->nmax;
          memory->destroy(xhold);
          memory->create(xhold,maxhold,3,"neigh:xhold");
        }
   -    for (i = 0; i < nlocal; i++) {
   -      xhold[i][0] = x[i][0];
   -      xhold[i][1] = x[i][1];
   -      xhold[i][2] = x[i][2];
   -    }
   +    const auto * _noalias const x = (const dbl3_t *) atom->x[0];
   +    auto * _noalias const xhold = (dbl3_t *) this->xhold[0];
   +    if (nlocal > 0) std::memcpy(xhold, x, sizeof(dbl3_t) * nlocal);
        if (boxcheck) {
          if (triclinic == 0) {
            boxlo_hold[0] = bboxlo[0];


.. _iter-0019-page-head-11db74893452: https://github.com/skilled-scipkg/lammps/commit/11db74893452
.. _iter-0019-guardrails-incumbent-08f5ab4a4e026b06a3d38c7e3ebd6ecca2a25635: https://github.com/skilled-scipkg/lammps/commit/08f5ab4a4e026b06a3d38c7e3ebd6ecca2a25635
.. _iter-0019-guardrails-candidate-11db748934527306a5e9d2d8c0c51eb47b03158a: https://github.com/skilled-scipkg/lammps/commit/11db748934527306a5e9d2d8c0c51eb47b03158a