Iteration 0018 — 08f5ab4a4e02 (accepted)
========================================


GitHub commit: `08f5ab4a4e02 <iter-0018-page-head-08f5ab4a4e02_>`_
Published branch: `fermilink-optimize/lammps-neighbor <https://github.com/skilled-scipkg/lammps/tree/fermilink-optimize%2Flammps-neighbor>`_

Change summary
--------------


Use explicit `jnext` linked-list traversal with bin/list `_noalias` aliases and delayed generic-only setup in the incumbent `NPairBin<1,1,0,0,0>` fast path

Acceptance rationale
--------------------


Correctness passed and weighted_median_neigh_seconds dropped from 0.31150 to 0.29094 (-6.60% vs incumbent), comfortably clearing the 2.0% promotion gate.

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
| incumbent commit | `8c29b124a2c6 <iter-0018-guardrails-incumbent-8c29b124a2c68bbc7efea2fbd01f26b5ad40288a_>`_ |
+------------------+--------------------------------------------------------------------------------------------+
| candidate commit | `08f5ab4a4e02 <iter-0018-guardrails-candidate-08f5ab4a4e026b06a3d38c7e3ebd6ecca2a25635_>`_ |
+------------------+--------------------------------------------------------------------------------------------+
| incumbent metric | 0.3115                                                                                     |
+------------------+--------------------------------------------------------------------------------------------+
| candidate metric | 0.29094                                                                                    |
+------------------+--------------------------------------------------------------------------------------------+
| baseline metric  | 0.36158                                                                                    |
+------------------+--------------------------------------------------------------------------------------------+
| Δ vs incumbent   | +6.600% (lower-is-better sign)                                                             |
+------------------+--------------------------------------------------------------------------------------------+
| changed files    | src/npair_bin.cpp                                                                          |
+------------------+--------------------------------------------------------------------------------------------+


Diffstat
--------


.. code-block:: text

    src/npair_bin.cpp | 60 +++++++++++++++++++++++++++++++++++++++----------------
    1 file changed, 43 insertions(+), 17 deletions(-)

Diff
----


:download:`download full diff <_diffs/iter_0018_08f5ab4a4e02.diff>`

.. code-block:: diff

   diff --git a/src/npair_bin.cpp b/src/npair_bin.cpp
   index 178d975ed0..dbf6866d61 100644
   --- a/src/npair_bin.cpp
   +++ b/src/npair_bin.cpp
   @@ -50,22 +50,19 @@ void NPairBin<1, 1, 0, 0, 0>::build(NeighList *list)
    
      const auto * _noalias const x = (const dbl3_t *) atom->x[0];
      const int * _noalias const type = atom->type;
   -  int * _noalias const mask = atom->mask;
      const tagint * _noalias const tag = atom->tag;
      tagint * _noalias const molecule = atom->molecule;
      tagint * const * _noalias const special = atom->special;
      int * const * _noalias const nspecial = atom->nspecial;
   +  const int * _noalias const atom2bin = this->atom2bin;
   +  const int * _noalias const bins = this->bins;
   +  const int * _noalias const binhead = this->binhead;
   +  const int * _noalias const stencil = this->stencil;
      int nlocal = atom->nlocal;
      if (includegroup) nlocal = atom->nfirst;
   -
   -  const bool moltemplate = (molecular == Atom::TEMPLATE);
   -  int *molindex = atom->molindex;
   -  int *molatom = atom->molatom;
   -  Molecule **onemols = atom->avec->onemols;
   -
   -  int *ilist = list->ilist;
   -  int *numneigh = list->numneigh;
   -  int **firstneigh = list->firstneigh;
   +  int * _noalias const ilist = list->ilist;
   +  int * _noalias const numneigh = list->numneigh;
   +  int ** _noalias const firstneigh = list->firstneigh;
      MyPage<int> *ipage = list->ipage;
    
      int inum = 0;
   @@ -90,13 +87,23 @@ void NPairBin<1, 1, 0, 0, 0>::build(NeighList *list)
    
          ibin = atom2bin[i];
    
   -      for (j = bins[i]; j >= 0; j = bins[j]) {
   +      for (j = bins[i]; j >= 0; ) {
   +        const int jnext = bins[j];
            const double *const xj = x[j];
            if (j >= nlocal) {
   -          if (xj[2] < ztmp) continue;
   +          if (xj[2] < ztmp) {
   +            j = jnext;
   +            continue;
   +          }
              if (xj[2] == ztmp) {
   -            if (xj[1] < ytmp) continue;
   -            if ((xj[1] == ytmp) && (xj[0] < xtmp)) continue;
   +            if (xj[1] < ytmp) {
   +              j = jnext;
   +              continue;
   +            }
   +            if ((xj[1] == ytmp) && (xj[0] < xtmp)) {
   +              j = jnext;
   +              continue;
   +            }
              }
            }
    
   @@ -107,9 +114,13 @@ void NPairBin<1, 1, 0, 0, 0>::build(NeighList *list)
            delz = ztmp - xj[2];
            rsq = delx * delx + dely * dely + delz * delz;
    
   -        if (rsq > cutneighsq_i[jtype]) continue;
   +        if (rsq > cutneighsq_i[jtype]) {
   +          j = jnext;
   +          continue;
   +        }
            if (molecule_i != molecule[j]) {
              neighptr[n++] = j;
   +          j = jnext;
              continue;
            }
    
   @@ -122,10 +133,13 @@ void NPairBin<1, 1, 0, 0, 0>::build(NeighList *list)
              neighptr[n++] = j;
            else if (which > 0)
              neighptr[n++] = j ^ (which << SBBITS);
   +
   +        j = jnext;
          }
    
          for (k = 1; k < nstencil; k++) {
   -        for (j = binhead[ibin + stencil[k]]; j >= 0; j = bins[j]) {
   +        for (j = binhead[ibin + stencil[k]]; j >= 0; ) {
   +          const int jnext = bins[j];
              const double *const xj = x[j];
              jtype = type[j];
    
   @@ -134,9 +148,13 @@ void NPairBin<1, 1, 0, 0, 0>::build(NeighList *list)
              delz = ztmp - xj[2];
              rsq = delx * delx + dely * dely + delz * delz;
    
   -          if (rsq > cutneighsq_i[jtype]) continue;
   +          if (rsq > cutneighsq_i[jtype]) {
   +            j = jnext;
   +            continue;
   +          }
              if (molecule_i != molecule[j]) {
                neighptr[n++] = j;
   +            j = jnext;
                continue;
              }
    
   @@ -149,6 +167,8 @@ void NPairBin<1, 1, 0, 0, 0>::build(NeighList *list)
                neighptr[n++] = j;
              else if (which > 0)
                neighptr[n++] = j ^ (which << SBBITS);
   +
   +          j = jnext;
            }
          }
    
   @@ -163,6 +183,12 @@ void NPairBin<1, 1, 0, 0, 0>::build(NeighList *list)
        return;
      }
    
   +  int * _noalias const mask = atom->mask;
   +  const bool moltemplate = (molecular == Atom::TEMPLATE);
   +  int * _noalias const molindex = atom->molindex;
   +  int * _noalias const molatom = atom->molatom;
   +  Molecule **onemols = atom->avec->onemols;
   +
      for (i = 0; i < nlocal; i++) {
        n = 0;
        neighptr = ipage->vget();


.. _iter-0018-page-head-08f5ab4a4e02: https://github.com/skilled-scipkg/lammps/commit/08f5ab4a4e02
.. _iter-0018-guardrails-incumbent-8c29b124a2c68bbc7efea2fbd01f26b5ad40288a: https://github.com/skilled-scipkg/lammps/commit/8c29b124a2c68bbc7efea2fbd01f26b5ad40288a
.. _iter-0018-guardrails-candidate-08f5ab4a4e026b06a3d38c7e3ebd6ecca2a25635: https://github.com/skilled-scipkg/lammps/commit/08f5ab4a4e026b06a3d38c7e3ebd6ecca2a25635