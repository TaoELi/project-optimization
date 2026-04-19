Iteration 0027 — 5840dac0d9f1 (accepted)
========================================


Change summary
--------------


Inline TIP4P pair/PPPM cache-hit fast paths while keeping cache-miss reconstruction in dedicated out-of-line slow helpers

Acceptance rationale
--------------------


Correctness passed and 70.984 beats incumbent 72.919 by 2.65%, clearing the 2% acceptance bar.

Guardrails & metrics
--------------------


+------------------+--------------------------------------------------------------------------------------------------------------------------------+
| field            | value                                                                                                                          |
+==================+================================================================================================================================+
| decision         | ACCEPTED                                                                                                                       |
+------------------+--------------------------------------------------------------------------------------------------------------------------------+
| correctness      | ok                                                                                                                             |
+------------------+--------------------------------------------------------------------------------------------------------------------------------+
| correctness mode | field_tolerances                                                                                                               |
+------------------+--------------------------------------------------------------------------------------------------------------------------------+
| hard reject      | no                                                                                                                             |
+------------------+--------------------------------------------------------------------------------------------------------------------------------+
| guardrail errors | 0                                                                                                                              |
+------------------+--------------------------------------------------------------------------------------------------------------------------------+
| incumbent commit | 4cbc0fa85603                                                                                                                   |
+------------------+--------------------------------------------------------------------------------------------------------------------------------+
| candidate commit | 5840dac0d9f1                                                                                                                   |
+------------------+--------------------------------------------------------------------------------------------------------------------------------+
| incumbent metric | 72.919                                                                                                                         |
+------------------+--------------------------------------------------------------------------------------------------------------------------------+
| candidate metric | 70.984                                                                                                                         |
+------------------+--------------------------------------------------------------------------------------------------------------------------------+
| baseline metric  | 94.717                                                                                                                         |
+------------------+--------------------------------------------------------------------------------------------------------------------------------+
| Δ vs incumbent   | +2.654% (lower-is-better sign)                                                                                                 |
+------------------+--------------------------------------------------------------------------------------------------------------------------------+
| changed files    | src/KSPACE/pair_lj_cut_tip4p_long.cpp, src/KSPACE/pair_lj_cut_tip4p_long.h, src/KSPACE/pppm_tip4p.cpp, src/KSPACE/pppm_tip4p.h |
+------------------+--------------------------------------------------------------------------------------------------------------------------------+


Diffstat
--------


.. code-block:: text

    src/KSPACE/pair_lj_cut_tip4p_long.cpp | 44 ++++++++++++++---------------------
    src/KSPACE/pair_lj_cut_tip4p_long.h   | 17 +++++++++++++-
    src/KSPACE/pppm_tip4p.cpp             |  8 +++----
    src/KSPACE/pppm_tip4p.h               | 14 ++++++++++-
    4 files changed, 49 insertions(+), 34 deletions(-)

Diff
----


:download:`download full diff <_diffs/iter_0027_5840dac0d9f1.diff>`

.. code-block:: diff

   diff --git a/src/KSPACE/pair_lj_cut_tip4p_long.cpp b/src/KSPACE/pair_lj_cut_tip4p_long.cpp
   index 7d6d4a8fa2..c79b0ac40f 100644
   --- a/src/KSPACE/pair_lj_cut_tip4p_long.cpp
   +++ b/src/KSPACE/pair_lj_cut_tip4p_long.cpp
   @@ -514,36 +514,26 @@ void PairLJCutTIP4PLong::reset_tip4p_cache_stamps()
    
    /* ---------------------------------------------------------------------- */
    
   -double *PairLJCutTIP4PLong::tip4p_site(int i, int &iH1, int &iH2, double **x, tagint *tag,
   -                                       int *type)
   +double *PairLJCutTIP4PLong::tip4p_site_slow(int i, int &iH1, int &iH2, double **x,
   +                                            tagint *tag, int *type)
    {
      Tip4pSiteCache &cache = tip4p_cache[i];
    
   -  if (cache.neigh_stamp != neigh_stamp) {
   -    iH1 = atom->map(tag[i] + 1);
   -    iH2 = atom->map(tag[i] + 2);
   -    if (iH1 == -1 || iH2 == -1)
   -      error->one(FLERR,"TIP4P hydrogen is missing");
   -    if (type[iH1] != typeH || type[iH2] != typeH)
   -      error->one(FLERR,"TIP4P hydrogen has incorrect atom type");
   -
   -    iH1 = domain->closest_image(i,iH1);
   -    iH2 = domain->closest_image(i,iH2);
   -    cache.iH1 = iH1;
   -    cache.iH2 = iH2;
   -    cache.neigh_stamp = neigh_stamp;
   -    cache.site_stamp = 0;
   -
   -  } else {
   -    iH1 = cache.iH1;
   -    iH2 = cache.iH2;
   -  }
   -
   -  if (cache.site_stamp != site_stamp) {
   -    compute_newsite(x[i],x[iH1],x[iH2],cache.xM);
   -    cache.site_stamp = site_stamp;
   -  }
   -
   +  iH1 = atom->map(tag[i] + 1);
   +  iH2 = atom->map(tag[i] + 2);
   +  if (iH1 == -1 || iH2 == -1)
   +    error->one(FLERR,"TIP4P hydrogen is missing");
   +  if (type[iH1] != typeH || type[iH2] != typeH)
   +    error->one(FLERR,"TIP4P hydrogen has incorrect atom type");
   +
   +  iH1 = domain->closest_image(i,iH1);
   +  iH2 = domain->closest_image(i,iH2);
   +  cache.iH1 = iH1;
   +  cache.iH2 = iH2;
   +  cache.neigh_stamp = neigh_stamp;
   +
   +  compute_newsite(x[i],x[iH1],x[iH2],cache.xM);
   +  cache.site_stamp = site_stamp;
      return cache.xM;
    }
    
   diff --git a/src/KSPACE/pair_lj_cut_tip4p_long.h b/src/KSPACE/pair_lj_cut_tip4p_long.h
   index 233fbb0dcc..fc02aa2800 100644
   --- a/src/KSPACE/pair_lj_cut_tip4p_long.h
   +++ b/src/KSPACE/pair_lj_cut_tip4p_long.h
   @@ -65,7 +65,22 @@ class PairLJCutTIP4PLong : public PairLJCutCoulLong {
                  double, double);
      void grow_tip4p_cache();
      void reset_tip4p_cache_stamps();
   -  double *tip4p_site(int, int &, int &, double **, tagint *, int *);
   +  double *tip4p_site(int i, int &iH1, int &iH2, double **x, tagint *tag, int *type)
   +  {
   +    Tip4pSiteCache &cache = tip4p_cache[i];
   +
   +    if (cache.neigh_stamp != neigh_stamp) return tip4p_site_slow(i,iH1,iH2,x,tag,type);
   +
   +    iH1 = cache.iH1;
   +    iH2 = cache.iH2;
   +    if (cache.site_stamp != site_stamp) {
   +      compute_newsite(x[i],x[iH1],x[iH2],cache.xM);
   +      cache.site_stamp = site_stamp;
   +    }
   +
   +    return cache.xM;
   +  }
   +  double *tip4p_site_slow(int, int &, int &, double **, tagint *, int *);
      void compute_newsite(double *, double *, double *, double *);
    };
    
   diff --git a/src/KSPACE/pppm_tip4p.cpp b/src/KSPACE/pppm_tip4p.cpp
   index 4ddee644a8..b8bf5f351f 100644
   --- a/src/KSPACE/pppm_tip4p.cpp
   +++ b/src/KSPACE/pppm_tip4p.cpp
   @@ -108,15 +108,13 @@ void PPPMTIP4P::reset_tip4p_cache_stamps()
    
    /* ---------------------------------------------------------------------- */
    
   -void PPPMTIP4P::find_M_cached(int i, int &iH1, int &iH2, double *&xM)
   +void PPPMTIP4P::find_M_cached_slow(int i, int &iH1, int &iH2, double *&xM)
    {
      if (atom->nmax > tip4p_cache_nmax) grow_tip4p_cache();
    
      Tip4pCache &cache = tip4p_cache[i];
   -  if (cache.stamp != tip4p_cache_stamp) {
   -    find_M(i,cache.iH1,cache.iH2,cache.xM);
   -    cache.stamp = tip4p_cache_stamp;
   -  }
   +  find_M(i,cache.iH1,cache.iH2,cache.xM);
   +  cache.stamp = tip4p_cache_stamp;
    
      iH1 = cache.iH1;
      iH2 = cache.iH2;
   diff --git a/src/KSPACE/pppm_tip4p.h b/src/KSPACE/pppm_tip4p.h
   index d907f96dd3..5ee607541f 100644
   --- a/src/KSPACE/pppm_tip4p.h
   +++ b/src/KSPACE/pppm_tip4p.h
   @@ -56,7 +56,19 @@ class PPPMTIP4P : public PPPM {
    
      void grow_tip4p_cache();
      void reset_tip4p_cache_stamps();
   -  void find_M_cached(int, int &, int &, double *&);
   +  void find_M_cached(int i, int &iH1, int &iH2, double *&xM)
   +  {
   +    Tip4pCache &cache = tip4p_cache[i];
   +    if (cache.stamp != tip4p_cache_stamp) {
   +      find_M_cached_slow(i,iH1,iH2,xM);
   +      return;
   +    }
   +
   +    iH1 = cache.iH1;
   +    iH2 = cache.iH2;
   +    xM = cache.xM;
   +  }
   +  void find_M_cached_slow(int, int &, int &, double *&);
      void find_M(int, int &, int &, double *);
    };
    
