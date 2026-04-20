Iteration 0011 — 8c29b124a2c6 (accepted)
========================================


GitHub commit: `8c29b124a2c6 <iter-0011-page-head-8c29b124a2c6_>`_
Published branch: `fermilink-optimize/lammps-neighbor <https://github.com/skilled-scipkg/lammps/tree/fermilink-optimize%2Flammps-neighbor>`_

Change summary
--------------


Use contiguous `_noalias` aliases for hot coordinate/type/tag/molecule arrays in the incumbent `NPairBin<1,1,0,0,0>::build()` path

Acceptance rationale
--------------------


Correctness held and the candidate cut weighted_median_neigh_seconds from 0.31969 to 0.31150 (-2.56% vs incumbent), clearing the 2.0% promotion gate.

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
| incumbent commit | `c4128d3970b2 <iter-0011-guardrails-incumbent-c4128d3970b2bed911c9977d594b37a3fa372a85_>`_ |
+------------------+--------------------------------------------------------------------------------------------+
| candidate commit | `8c29b124a2c6 <iter-0011-guardrails-candidate-8c29b124a2c68bbc7efea2fbd01f26b5ad40288a_>`_ |
+------------------+--------------------------------------------------------------------------------------------+
| incumbent metric | 0.31969                                                                                    |
+------------------+--------------------------------------------------------------------------------------------+
| candidate metric | 0.3115                                                                                     |
+------------------+--------------------------------------------------------------------------------------------+
| baseline metric  | 0.36158                                                                                    |
+------------------+--------------------------------------------------------------------------------------------+
| Δ vs incumbent   | +2.562% (lower-is-better sign)                                                             |
+------------------+--------------------------------------------------------------------------------------------+
| changed files    | src/npair_bin.cpp                                                                          |
+------------------+--------------------------------------------------------------------------------------------+


Diffstat
--------


.. code-block:: text

    src/npair_bin.cpp | 18 ++++++++++--------
    1 file changed, 10 insertions(+), 8 deletions(-)

Diff
----


:download:`download full diff <_diffs/iter_0011_8c29b124a2c6.diff>`

.. code-block:: diff

   diff --git a/src/npair_bin.cpp b/src/npair_bin.cpp
   index 42917b2dda..178d975ed0 100644
   --- a/src/npair_bin.cpp
   +++ b/src/npair_bin.cpp
   @@ -46,13 +46,15 @@ void NPairBin<1, 1, 0, 0, 0>::build(NeighList *list)
      int *neighptr;
      double xtmp, ytmp, ztmp, delx, dely, delz, rsq;
    
   -  double **x = atom->x;
   -  int *type = atom->type;
   -  int *mask = atom->mask;
   -  tagint *tag = atom->tag;
   -  tagint *molecule = atom->molecule;
   -  tagint **special = atom->special;
   -  int **nspecial = atom->nspecial;
   +  using dbl3_t = double[3];
   +
   +  const auto * _noalias const x = (const dbl3_t *) atom->x[0];
   +  const int * _noalias const type = atom->type;
   +  int * _noalias const mask = atom->mask;
   +  const tagint * _noalias const tag = atom->tag;
   +  tagint * _noalias const molecule = atom->molecule;
   +  tagint * const * _noalias const special = atom->special;
   +  int * const * _noalias const nspecial = atom->nspecial;
      int nlocal = atom->nlocal;
      if (includegroup) nlocal = atom->nfirst;
    
   @@ -78,7 +80,7 @@ void NPairBin<1, 1, 0, 0, 0>::build(NeighList *list)
          neighptr = ipage->vget();
    
          itype = type[i];
   -      const double *const cutneighsq_i = cutneighsq[itype];
   +      const double * _noalias const cutneighsq_i = cutneighsq[itype];
          const tagint molecule_i = molecule[i];
          const tagint *const special_i = special[i];
          const int *const nspecial_i = nspecial[i];


.. _iter-0011-page-head-8c29b124a2c6: https://github.com/skilled-scipkg/lammps/commit/8c29b124a2c6
.. _iter-0011-guardrails-incumbent-c4128d3970b2bed911c9977d594b37a3fa372a85: https://github.com/skilled-scipkg/lammps/commit/c4128d3970b2bed911c9977d594b37a3fa372a85
.. _iter-0011-guardrails-candidate-8c29b124a2c68bbc7efea2fbd01f26b5ad40288a: https://github.com/skilled-scipkg/lammps/commit/8c29b124a2c68bbc7efea2fbd01f26b5ad40288a