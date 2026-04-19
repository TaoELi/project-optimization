Iteration 0015 — 01ed31155d08 (accepted)
========================================


Change summary
--------------


Tighten non-O Coulomb outer gating to `cut_coulsq` and reuse precomputed H-coordinate sums for O-path virial accumulation in `pair_lj_cut_tip4p_long::eval()` to cut unnecessary checks and arithmetic.

Acceptance rationale
--------------------


Candidate improved weighted_median_wall_seconds_per_100_steps by +2.5689% vs incumbent (0.271962072991 -> 0.264975531201) with correctness and guardrails passing and no hard reject.

Guardrails & metrics
--------------------


+------------------+---------------------------------------+
| field            | value                                 |
+==================+=======================================+
| decision         | ACCEPTED                              |
+------------------+---------------------------------------+
| correctness      | ok                                    |
+------------------+---------------------------------------+
| correctness mode | field_tolerances                      |
+------------------+---------------------------------------+
| hard reject      | no                                    |
+------------------+---------------------------------------+
| guardrail errors | 0                                     |
+------------------+---------------------------------------+
| incumbent commit | ecf63fa05a3a                          |
+------------------+---------------------------------------+
| candidate commit | 01ed31155d08                          |
+------------------+---------------------------------------+
| incumbent metric | 0.271962                              |
+------------------+---------------------------------------+
| candidate metric | 0.264976                              |
+------------------+---------------------------------------+
| baseline metric  | 0.34835                               |
+------------------+---------------------------------------+
| Δ vs incumbent   | +2.569% (lower-is-better sign)        |
+------------------+---------------------------------------+
| changed files    | src/KSPACE/pair_lj_cut_tip4p_long.cpp |
+------------------+---------------------------------------+


Diffstat
--------


.. code-block:: text

    src/KSPACE/pair_lj_cut_tip4p_long.cpp | 36 ++++++++++++++++++++++-------------
    1 file changed, 23 insertions(+), 13 deletions(-)

Diff
----


:download:`download full diff <_diffs/iter_0015_01ed31155d08.diff>`

.. code-block:: diff

   diff --git a/src/KSPACE/pair_lj_cut_tip4p_long.cpp b/src/KSPACE/pair_lj_cut_tip4p_long.cpp
   index 34aa06f95c..a436d3aafd 100644
   --- a/src/KSPACE/pair_lj_cut_tip4p_long.cpp
   +++ b/src/KSPACE/pair_lj_cut_tip4p_long.cpp
   @@ -180,6 +180,7 @@ void PairLJCutTIP4PLong::eval()
      double fiHxtmp,fiHytmp,fiHztmp;
      double xiH1x = 0.0, xiH1y = 0.0, xiH1z = 0.0;
      double xiH2x = 0.0, xiH2y = 0.0, xiH2z = 0.0;
   +  double xiHsumx = 0.0, xiHsumy = 0.0, xiHsumz = 0.0;
    
      inum = list->inum;
      ilist = list->ilist;
   @@ -222,6 +223,9 @@ void PairLJCutTIP4PLong::eval()
            xiH2x = xiH2[0];
            xiH2y = xiH2[1];
            xiH2z = xiH2[2];
   +        xiHsumx = xiH1x + xiH2x;
   +        xiHsumy = xiH1y + xiH2y;
   +        xiHsumz = xiH1z + xiH2z;
          }
        } else x1 = x[i];
    
   @@ -272,7 +276,10 @@ void PairLJCutTIP4PLong::eval()
          // adjust rsq and delxyz for off-site O charge(s) if necessary
          // but only if they are within reach
    
   -      if (rsq < cut_coulsqplus) {
   +      const double cut_coulsq_gate = (i_is_O || j_is_O) ? cut_coulsqplus :
   +        cut_coulsq;
   +
   +      if (rsq < cut_coulsq_gate) {
            if (i_is_O) {
    
              // if atom J = water O, set x2 = offset charge site
   @@ -387,12 +394,12 @@ void PairLJCutTIP4PLong::eval()
                fiHztmp += fHz;
    
                if (VFLAG) {
   -              v[0] = xtmp*fOx + xiH1x*fHx + xiH2x*fHx;
   -              v[1] = ytmp*fOy + xiH1y*fHy + xiH2y*fHy;
   -              v[2] = ztmp*fOz + xiH1z*fHz + xiH2z*fHz;
   -              v[3] = xtmp*fOy + xiH1x*fHy + xiH2x*fHy;
   -              v[4] = xtmp*fOz + xiH1x*fHz + xiH2x*fHz;
   -              v[5] = ytmp*fOz + xiH1y*fHz + xiH2y*fHz;
   +              v[0] = xtmp*fOx + xiHsumx*fHx;
   +              v[1] = ytmp*fOy + xiHsumy*fHy;
   +              v[2] = ztmp*fOz + xiHsumz*fHz;
   +              v[3] = xtmp*fOy + xiHsumx*fHy;
   +              v[4] = xtmp*fOz + xiHsumx*fHz;
   +              v[5] = ytmp*fOz + xiHsumy*fHz;
                }
                if (EVFLAG) {
                  vlist[n++] = i;
   @@ -446,12 +453,15 @@ void PairLJCutTIP4PLong::eval()
                if (VFLAG) {
                  xH1 = x[jH1];
                  xH2 = x[jH2];
   -              v[0] += xj[0]*fOx + xH1[0]*fHx + xH2[0]*fHx;
   -              v[1] += xj[1]*fOy + xH1[1]*fHy + xH2[1]*fHy;
   -              v[2] += xj[2]*fOz + xH1[2]*fHz + xH2[2]*fHz;
   -              v[3] += xj[0]*fOy + xH1[0]*fHy + xH2[0]*fHy;
   -              v[4] += xj[0]*fOz + xH1[0]*fHz + xH2[0]*fHz;
   -              v[5] += xj[1]*fOz + xH1[1]*fHz + xH2[1]*fHz;
   +              const double xHsumx = xH1[0] + xH2[0];
   +              const double xHsumy = xH1[1] + xH2[1];
   +              const double xHsumz = xH1[2] + xH2[2];
   +              v[0] += xj[0]*fOx + xHsumx*fHx;
   +              v[1] += xj[1]*fOy + xHsumy*fHy;
   +              v[2] += xj[2]*fOz + xHsumz*fHz;
   +              v[3] += xj[0]*fOy + xHsumx*fHy;
   +              v[4] += xj[0]*fOz + xHsumx*fHz;
   +              v[5] += xj[1]*fOz + xHsumy*fHz;
                }
                if (EVFLAG) {
                  vlist[n++] = j;
