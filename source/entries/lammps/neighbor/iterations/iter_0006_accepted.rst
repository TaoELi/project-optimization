Iteration 0006 — c4128d3970b2 (accepted)
========================================


GitHub commit: `c4128d3970b2 <iter-0006-page-head-c4128d3970b2_>`_
Published branch: `fermilink-optimize/lammps-neighbor <https://github.com/skilled-scipkg/lammps/tree/fermilink-optimize%2Flammps-neighbor>`_

Change summary
--------------


Refine the incumbent `NPairBin<1,1,0,0,0>::build()` fast path by early-accepting different-molecule pairs and skipping `minimum_image_check()` for owned same-molecule specials when `special_flag[1..3] == 2`.

Acceptance rationale
--------------------


Correctness held and weighted_median_neigh_seconds improved from 0.32842 to 0.31969 (-2.66% vs incumbent), clearing the benchmark contract’s 2.0% promotion gate.

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
| incumbent commit | `380bf2b0146b <iter-0006-guardrails-incumbent-380bf2b0146b31a56d1f86c034a73cd0d4b09e43_>`_ |
+------------------+--------------------------------------------------------------------------------------------+
| candidate commit | `c4128d3970b2 <iter-0006-guardrails-candidate-c4128d3970b2bed911c9977d594b37a3fa372a85_>`_ |
+------------------+--------------------------------------------------------------------------------------------+
| incumbent metric | 0.32842                                                                                    |
+------------------+--------------------------------------------------------------------------------------------+
| candidate metric | 0.31969                                                                                    |
+------------------+--------------------------------------------------------------------------------------------+
| baseline metric  | 0.36158                                                                                    |
+------------------+--------------------------------------------------------------------------------------------+
| Δ vs incumbent   | +2.658% (lower-is-better sign)                                                             |
+------------------+--------------------------------------------------------------------------------------------+
| changed files    | src/npair_bin.cpp                                                                          |
+------------------+--------------------------------------------------------------------------------------------+


Diffstat
--------


.. code-block:: text

    src/npair_bin.cpp | 16 ++++++++++++++++
    1 file changed, 16 insertions(+)

Diff
----


:download:`download full diff <_diffs/iter_0006_c4128d3970b2.diff>`

.. code-block:: diff

   diff --git a/src/npair_bin.cpp b/src/npair_bin.cpp
   index bd1fc1d338..42917b2dda 100644
   --- a/src/npair_bin.cpp
   +++ b/src/npair_bin.cpp
   @@ -70,12 +70,16 @@ void NPairBin<1, 1, 0, 0, 0>::build(NeighList *list)
      ipage->reset();
    
      if (!exclude && molecular == Atom::MOLECULAR) {
   +    const bool all_special_encoded =
   +      (special_flag[1] == 2) && (special_flag[2] == 2) && (special_flag[3] == 2);
   +
        for (i = 0; i < nlocal; i++) {
          n = 0;
          neighptr = ipage->vget();
    
          itype = type[i];
          const double *const cutneighsq_i = cutneighsq[itype];
   +      const tagint molecule_i = molecule[i];
          const tagint *const special_i = special[i];
          const int *const nspecial_i = nspecial[i];
          xtmp = x[i][0];
   @@ -102,10 +106,16 @@ void NPairBin<1, 1, 0, 0, 0>::build(NeighList *list)
            rsq = delx * delx + dely * dely + delz * delz;
    
            if (rsq > cutneighsq_i[jtype]) continue;
   +        if (molecule_i != molecule[j]) {
   +          neighptr[n++] = j;
   +          continue;
   +        }
    
            which = find_special(special_i, nspecial_i, tag[j]);
            if (which == 0)
              neighptr[n++] = j;
   +        else if (all_special_encoded && (j < nlocal))
   +          neighptr[n++] = j ^ (which << SBBITS);
            else if (domain->minimum_image_check(delx, dely, delz))
              neighptr[n++] = j;
            else if (which > 0)
   @@ -123,10 +133,16 @@ void NPairBin<1, 1, 0, 0, 0>::build(NeighList *list)
              rsq = delx * delx + dely * dely + delz * delz;
    
              if (rsq > cutneighsq_i[jtype]) continue;
   +          if (molecule_i != molecule[j]) {
   +            neighptr[n++] = j;
   +            continue;
   +          }
    
              which = find_special(special_i, nspecial_i, tag[j]);
              if (which == 0)
                neighptr[n++] = j;
   +          else if (all_special_encoded && (j < nlocal))
   +            neighptr[n++] = j ^ (which << SBBITS);
              else if (domain->minimum_image_check(delx, dely, delz))
                neighptr[n++] = j;
              else if (which > 0)


.. _iter-0006-page-head-c4128d3970b2: https://github.com/skilled-scipkg/lammps/commit/c4128d3970b2
.. _iter-0006-guardrails-incumbent-380bf2b0146b31a56d1f86c034a73cd0d4b09e43: https://github.com/skilled-scipkg/lammps/commit/380bf2b0146b31a56d1f86c034a73cd0d4b09e43
.. _iter-0006-guardrails-candidate-c4128d3970b2bed911c9977d594b37a3fa372a85: https://github.com/skilled-scipkg/lammps/commit/c4128d3970b2bed911c9977d594b37a3fa372a85