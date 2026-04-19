Iteration 0011 — 2a0b902a8b5a (accepted)
========================================


Change summary
--------------


Tighten exact oxygen-to-nonoxygen Coulomb prechecks and add a dedicated water-H `i` fast path in `pair_lj_cut_tip4p_long` that skips impossible LJ work while preserving the incumbent TIP4P force/tally order.

Acceptance rationale
--------------------


Correctness passed and 79.88 beats incumbent 81.885 by 2.45%, clearing the 2% acceptance bar.

Guardrails & metrics
--------------------


+------------------+----------------------------------------------------------------------------+
| field            | value                                                                      |
+==================+============================================================================+
| decision         | ACCEPTED                                                                   |
+------------------+----------------------------------------------------------------------------+
| correctness      | ok                                                                         |
+------------------+----------------------------------------------------------------------------+
| correctness mode | field_tolerances                                                           |
+------------------+----------------------------------------------------------------------------+
| hard reject      | no                                                                         |
+------------------+----------------------------------------------------------------------------+
| guardrail errors | 0                                                                          |
+------------------+----------------------------------------------------------------------------+
| incumbent commit | 3a45586c8afa                                                               |
+------------------+----------------------------------------------------------------------------+
| candidate commit | 2a0b902a8b5a                                                               |
+------------------+----------------------------------------------------------------------------+
| incumbent metric | 81.885                                                                     |
+------------------+----------------------------------------------------------------------------+
| candidate metric | 79.88                                                                      |
+------------------+----------------------------------------------------------------------------+
| baseline metric  | 94.717                                                                     |
+------------------+----------------------------------------------------------------------------+
| Δ vs incumbent   | +2.449% (lower-is-better sign)                                             |
+------------------+----------------------------------------------------------------------------+
| changed files    | src/KSPACE/pair_lj_cut_tip4p_long.cpp, src/KSPACE/pair_lj_cut_tip4p_long.h |
+------------------+----------------------------------------------------------------------------+


Diffstat
--------


.. code-block:: text

    src/KSPACE/pair_lj_cut_tip4p_long.cpp | 97 ++++++++++++++++++++++-------------
    src/KSPACE/pair_lj_cut_tip4p_long.h   |  4 +-
    2 files changed, 63 insertions(+), 38 deletions(-)

Diff
----


:download:`download full diff <_diffs/iter_0011_2a0b902a8b5a.diff>`

.. code-block:: diff

   diff --git a/src/KSPACE/pair_lj_cut_tip4p_long.cpp b/src/KSPACE/pair_lj_cut_tip4p_long.cpp
   index 5534584ad0..c7b794571f 100644
   --- a/src/KSPACE/pair_lj_cut_tip4p_long.cpp
   +++ b/src/KSPACE/pair_lj_cut_tip4p_long.cpp
   @@ -143,21 +143,31 @@ void PairLJCutTIP4PLong::eval()
        if (itype == typeO) {
          int iH1, iH2;
          double *x1 = tip4p_site(i,iH1,iH2,x,tag,type);
   -      eval_i<CTABLE,EVFLAG,EFLAG,VFLAG,1>(i,itype,qtmp,xtmp,ytmp,ztmp,iH1,iH2,x1,jlist,jnum,
   -                                          nlocal,cut_coulsqplus);
   +      const double site_delx = x1[0] - xtmp;
   +      const double site_dely = x1[1] - ytmp;
   +      const double site_delz = x1[2] - ztmp;
   +      const double cut_coul_plus_i_site =
   +        cut_coul + sqrt(site_delx*site_delx + site_dely*site_dely + site_delz*site_delz);
   +      const double cut_coulsq_i_h = cut_coul_plus_i_site * cut_coul_plus_i_site;
   +      eval_i<CTABLE,EVFLAG,EFLAG,VFLAG,1,0>(i,itype,qtmp,xtmp,ytmp,ztmp,iH1,iH2,x1,jlist,jnum,
   +                                            nlocal,cut_coulsq_i_h,cut_coulsqplus);
   +    } else if (itype == typeH) {
   +      eval_i<CTABLE,EVFLAG,EFLAG,VFLAG,0,1>(i,itype,qtmp,xtmp,ytmp,ztmp,-1,-1,x[i],jlist,jnum,
   +                                            nlocal,cut_coulsq,cut_coulsqplus);
        } else {
   -      eval_i<CTABLE,EVFLAG,EFLAG,VFLAG,0>(i,itype,qtmp,xtmp,ytmp,ztmp,-1,-1,x[i],jlist,jnum,
   -                                          nlocal,cut_coulsqplus);
   +      eval_i<CTABLE,EVFLAG,EFLAG,VFLAG,0,0>(i,itype,qtmp,xtmp,ytmp,ztmp,-1,-1,x[i],jlist,jnum,
   +                                            nlocal,cut_coulsq,cut_coulsqplus);
        }
      }
    }
    
    /* ---------------------------------------------------------------------- */
    
   -template <int CTABLE, int EVFLAG, int EFLAG, int VFLAG, int I_WATER>
   +template <int CTABLE, int EVFLAG, int EFLAG, int VFLAG, int I_WATER, int I_HYDROGEN>
    void PairLJCutTIP4PLong::eval_i(int i, int itype, double qtmp, double xtmp, double ytmp,
                                    double ztmp, int iH1, int iH2, double *x1, int *jlist,
   -                                int jnum, int nlocal, double cut_coulsqplus)
   +                                int jnum, int nlocal, double cut_coulsq_i_h,
   +                                double cut_coulsqplus)
    {
      int j, jj, jtype, itable, key = 0, n = 0, jH1 = -1, jH2 = -1;
      int vlist[6];
   @@ -178,12 +188,21 @@ void PairLJCutTIP4PLong::eval_i(int i, int itype, double qtmp, double xtmp, doub
      const double qqrd2e = force->qqrd2e;
      const double half_alpha = 0.5 * alpha;
      const double one_minus_alpha = 1.0 - alpha;
   -  double *cut_ljsq_i = cut_ljsq[itype];
   -  double *lj1_i = lj1[itype];
   -  double *lj2_i = lj2[itype];
   -  double *lj3_i = lj3[itype];
   -  double *lj4_i = lj4[itype];
   -  double *offset_i = offset[itype];
   +  double *cut_ljsq_i = nullptr;
   +  double *lj1_i = nullptr;
   +  double *lj2_i = nullptr;
   +  double *lj3_i = nullptr;
   +  double *lj4_i = nullptr;
   +  double *offset_i = nullptr;
   +
   +  if constexpr (!I_HYDROGEN) {
   +    cut_ljsq_i = cut_ljsq[itype];
   +    lj1_i = lj1[itype];
   +    lj2_i = lj2[itype];
   +    lj3_i = lj3[itype];
   +    lj4_i = lj4[itype];
   +    offset_i = offset[itype];
   +  }
    
      for (jj = 0; jj < jnum; jj++) {
        j = jlist[jj];
   @@ -196,39 +215,45 @@ void PairLJCutTIP4PLong::eval_i(int i, int itype, double qtmp, double xtmp, doub
        delz = ztmp - x[j][2];
        rsq = delx*delx + dely*dely + delz*delz;
        jtype = type[j];
   +    const bool j_is_oxygen = (jtype == typeO);
   +
   +    if constexpr (!I_HYDROGEN) {
   +      if (rsq < cut_ljsq_i[jtype]) {
   +        r2inv = 1.0/rsq;
   +        r6inv = r2inv*r2inv*r2inv;
   +        forcelj = r6inv * (lj1_i[jtype]*r6inv - lj2_i[jtype]);
   +        forcelj *= factor_lj * r2inv;
   +
   +        f[i][0] += delx*forcelj;
   +        f[i][1] += dely*forcelj;
   +        f[i][2] += delz*forcelj;
   +        f[j][0] -= delx*forcelj;
   +        f[j][1] -= dely*forcelj;
   +        f[j][2] -= delz*forcelj;
   +
   +        if constexpr (EFLAG) {
   +          evdwl = r6inv*(lj3_i[jtype]*r6inv-lj4_i[jtype]) - offset_i[jtype];
   +          evdwl *= factor_lj;
   +        } else evdwl = 0.0;
    
   -    if (rsq < cut_ljsq_i[jtype]) {
   -      r2inv = 1.0/rsq;
   -      r6inv = r2inv*r2inv*r2inv;
   -      forcelj = r6inv * (lj1_i[jtype]*r6inv - lj2_i[jtype]);
   -      forcelj *= factor_lj * r2inv;
   -
   -      f[i][0] += delx*forcelj;
   -      f[i][1] += dely*forcelj;
   -      f[i][2] += delz*forcelj;
   -      f[j][0] -= delx*forcelj;
   -      f[j][1] -= dely*forcelj;
   -      f[j][2] -= delz*forcelj;
   -
   -      if constexpr (EFLAG) {
   -        evdwl = r6inv*(lj3_i[jtype]*r6inv-lj4_i[jtype]) - offset_i[jtype];
   -        evdwl *= factor_lj;
   -      } else evdwl = 0.0;
   -
   -      if constexpr (EVFLAG)
   -        ev_tally(i,j,nlocal,1,evdwl,0.0,forcelj,delx,dely,delz);
   +        if constexpr (EVFLAG)
   +          ev_tally(i,j,nlocal,1,evdwl,0.0,forcelj,delx,dely,delz);
   +      }
        }
    
   -    if (rsq < cut_coulsqplus) {
   +    double cut_coulsq_precheck = cut_coulsqplus;
   +    if (!j_is_oxygen) cut_coulsq_precheck = cut_coulsq_i_h;
   +
   +    if (rsq < cut_coulsq_precheck) {
          if constexpr (I_WATER) {
   -        if (jtype == typeO) x2 = tip4p_site(j,jH1,jH2,x,tag,type);
   +        if (j_is_oxygen) x2 = tip4p_site(j,jH1,jH2,x,tag,type);
            else x2 = x[j];
    
            delx = x1[0] - x2[0];
            dely = x1[1] - x2[1];
            delz = x1[2] - x2[2];
            rsq = delx*delx + dely*dely + delz*delz;
   -      } else if (jtype == typeO) {
   +      } else if (j_is_oxygen) {
            x2 = tip4p_site(j,jH1,jH2,x,tag,type);
            delx = x1[0] - x2[0];
            dely = x1[1] - x2[1];
   @@ -328,7 +353,7 @@ void PairLJCutTIP4PLong::eval_i(int i, int itype, double qtmp, double xtmp, doub
              }
            }
    
   -        if (jtype != typeO) {
   +        if (!j_is_oxygen) {
              f[j][0] -= delx * cforce;
              f[j][1] -= dely * cforce;
              f[j][2] -= delz * cforce;
   diff --git a/src/KSPACE/pair_lj_cut_tip4p_long.h b/src/KSPACE/pair_lj_cut_tip4p_long.h
   index 52b7d39161..b97b0319cb 100644
   --- a/src/KSPACE/pair_lj_cut_tip4p_long.h
   +++ b/src/KSPACE/pair_lj_cut_tip4p_long.h
   @@ -50,9 +50,9 @@ class PairLJCutTIP4PLong : public PairLJCutCoulLong {
      double **newsite;    // locations of charge sites
    
      template <int, int, int, int> void eval();
   -  template <int, int, int, int, int>
   +  template <int, int, int, int, int, int>
      void eval_i(int, int, double, double, double, double, int, int, double *, int *, int, int,
   -              double);
   +              double, double);
      double *tip4p_site(int, int &, int &, double **, tagint *, int *);
      void compute_newsite(double *, double *, double *, double *);
    };
