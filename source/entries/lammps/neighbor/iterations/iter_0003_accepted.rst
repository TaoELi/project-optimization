Iteration 0003 — 380bf2b0146b (accepted)
========================================


GitHub commit: `380bf2b0146b <iter-0003-page-head-380bf2b0146b_>`_
Published branch: `fermilink-optimize/lammps-neighbor <https://github.com/skilled-scipkg/lammps/tree/fermilink-optimize%2Flammps-neighbor>`_

Change summary
--------------


Add a benchmark-specific fast path in `NPairBin<1,1,0,0,0>::build()` for standard molecular systems without neighbor exclusions, removing repeated `special_i` setup and `exclude` checks from the hot orthogonal half/bin/newton loops.

Acceptance rationale
--------------------


Candidate improves weighted_median_neigh_seconds from 0.33535 to 0.32842 with correctness intact, so it becomes the new incumbent.

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
| incumbent commit | `32b495fdb66b <iter-0003-guardrails-incumbent-32b495fdb66bfedb1ecf9eda9c0a1ec6d7bb6bb4_>`_ |
+------------------+--------------------------------------------------------------------------------------------+
| candidate commit | `380bf2b0146b <iter-0003-guardrails-candidate-380bf2b0146b31a56d1f86c034a73cd0d4b09e43_>`_ |
+------------------+--------------------------------------------------------------------------------------------+
| incumbent metric | 0.33535                                                                                    |
+------------------+--------------------------------------------------------------------------------------------+
| candidate metric | 0.32842                                                                                    |
+------------------+--------------------------------------------------------------------------------------------+
| baseline metric  | 0.36158                                                                                    |
+------------------+--------------------------------------------------------------------------------------------+
| Δ vs incumbent   | +2.066% (lower-is-better sign)                                                             |
+------------------+--------------------------------------------------------------------------------------------+
| changed files    | src/npair_bin.cpp                                                                          |
+------------------+--------------------------------------------------------------------------------------------+


Diffstat
--------


.. code-block:: text

    src/npair_bin.cpp | 76 +++++++++++++++++++++++++++++++++++++++++++++++++++++++
    1 file changed, 76 insertions(+)

Diff
----


:download:`download full diff <_diffs/iter_0003_380bf2b0146b.diff>`

.. code-block:: diff

   diff --git a/src/npair_bin.cpp b/src/npair_bin.cpp
   index 490bb45169..bd1fc1d338 100644
   --- a/src/npair_bin.cpp
   +++ b/src/npair_bin.cpp
   @@ -69,6 +69,82 @@ void NPairBin<1, 1, 0, 0, 0>::build(NeighList *list)
      int inum = 0;
      ipage->reset();
    
   +  if (!exclude && molecular == Atom::MOLECULAR) {
   +    for (i = 0; i < nlocal; i++) {
   +      n = 0;
   +      neighptr = ipage->vget();
   +
   +      itype = type[i];
   +      const double *const cutneighsq_i = cutneighsq[itype];
   +      const tagint *const special_i = special[i];
   +      const int *const nspecial_i = nspecial[i];
   +      xtmp = x[i][0];
   +      ytmp = x[i][1];
   +      ztmp = x[i][2];
   +
   +      ibin = atom2bin[i];
   +
   +      for (j = bins[i]; j >= 0; j = bins[j]) {
   +        const double *const xj = x[j];
   +        if (j >= nlocal) {
   +          if (xj[2] < ztmp) continue;
   +          if (xj[2] == ztmp) {
   +            if (xj[1] < ytmp) continue;
   +            if ((xj[1] == ytmp) && (xj[0] < xtmp)) continue;
   +          }
   +        }
   +
   +        jtype = type[j];
   +
   +        delx = xtmp - xj[0];
   +        dely = ytmp - xj[1];
   +        delz = ztmp - xj[2];
   +        rsq = delx * delx + dely * dely + delz * delz;
   +
   +        if (rsq > cutneighsq_i[jtype]) continue;
   +
   +        which = find_special(special_i, nspecial_i, tag[j]);
   +        if (which == 0)
   +          neighptr[n++] = j;
   +        else if (domain->minimum_image_check(delx, dely, delz))
   +          neighptr[n++] = j;
   +        else if (which > 0)
   +          neighptr[n++] = j ^ (which << SBBITS);
   +      }
   +
   +      for (k = 1; k < nstencil; k++) {
   +        for (j = binhead[ibin + stencil[k]]; j >= 0; j = bins[j]) {
   +          const double *const xj = x[j];
   +          jtype = type[j];
   +
   +          delx = xtmp - xj[0];
   +          dely = ytmp - xj[1];
   +          delz = ztmp - xj[2];
   +          rsq = delx * delx + dely * dely + delz * delz;
   +
   +          if (rsq > cutneighsq_i[jtype]) continue;
   +
   +          which = find_special(special_i, nspecial_i, tag[j]);
   +          if (which == 0)
   +            neighptr[n++] = j;
   +          else if (domain->minimum_image_check(delx, dely, delz))
   +            neighptr[n++] = j;
   +          else if (which > 0)
   +            neighptr[n++] = j ^ (which << SBBITS);
   +        }
   +      }
   +
   +      ilist[inum++] = i;
   +      firstneigh[i] = neighptr;
   +      numneigh[i] = n;
   +      ipage->vgot(n);
   +      if (ipage->status()) error->one(FLERR, Error::NOLASTLINE, "Neighbor list overflow, boost neigh_modify one" + utils::errorurl(36));
   +    }
   +
   +    list->inum = inum;
   +    return;
   +  }
   +
      for (i = 0; i < nlocal; i++) {
        n = 0;
        neighptr = ipage->vget();


.. _iter-0003-page-head-380bf2b0146b: https://github.com/skilled-scipkg/lammps/commit/380bf2b0146b
.. _iter-0003-guardrails-incumbent-32b495fdb66bfedb1ecf9eda9c0a1ec6d7bb6bb4: https://github.com/skilled-scipkg/lammps/commit/32b495fdb66bfedb1ecf9eda9c0a1ec6d7bb6bb4
.. _iter-0003-guardrails-candidate-380bf2b0146b31a56d1f86c034a73cd0d4b09e43: https://github.com/skilled-scipkg/lammps/commit/380bf2b0146b31a56d1f86c034a73cd0d4b09e43