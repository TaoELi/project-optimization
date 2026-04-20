Iteration 0026 — 4cbc0fa85603 (accepted)
========================================


GitHub commit: `4cbc0fa85603 <iter-0026-page-head-4cbc0fa85603_>`_
Published branch: `fermilink-optimize/lammps-tip4p <https://github.com/skilled-scipkg/lammps/tree/fermilink-optimize%2Flammps-tip4p>`_

Change summary
--------------


Contiguous generation-stamped TIP4P pair/PPPM caches plus a pure-water lj/cut/tip4p/long pair fast path that keeps exact oxygen-site cutoffs and incumbent force/tally ordering

Acceptance rationale
--------------------


Correctness passed and the candidate cuts weighted_median_pair_plus_kspace_seconds from 76.248 to 72.919, a 4.37% improvement over incumbent 2599e1f9f4c8 that clears the 2% acceptance bar.

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
| incumbent commit | `2599e1f9f4c8 <iter-0026-guardrails-incumbent-2599e1f9f4c8bb574886ac80f8e67052265a0ad4_>`_                                     |
+------------------+--------------------------------------------------------------------------------------------------------------------------------+
| candidate commit | `4cbc0fa85603 <iter-0026-guardrails-candidate-4cbc0fa856036c692930cd76775c671201e183cc_>`_                                     |
+------------------+--------------------------------------------------------------------------------------------------------------------------------+
| incumbent metric | 76.248                                                                                                                         |
+------------------+--------------------------------------------------------------------------------------------------------------------------------+
| candidate metric | 72.919                                                                                                                         |
+------------------+--------------------------------------------------------------------------------------------------------------------------------+
| baseline metric  | 94.717                                                                                                                         |
+------------------+--------------------------------------------------------------------------------------------------------------------------------+
| Δ vs incumbent   | +4.366% (lower-is-better sign)                                                                                                 |
+------------------+--------------------------------------------------------------------------------------------------------------------------------+
| changed files    | src/KSPACE/pair_lj_cut_tip4p_long.cpp, src/KSPACE/pair_lj_cut_tip4p_long.h, src/KSPACE/pppm_tip4p.cpp, src/KSPACE/pppm_tip4p.h |
+------------------+--------------------------------------------------------------------------------------------------------------------------------+


Diffstat
--------


.. code-block:: text

    src/KSPACE/pair_lj_cut_tip4p_long.cpp | 160 +++++++++++++++++++++++++---------
    src/KSPACE/pair_lj_cut_tip4p_long.h   |  20 ++++-
    src/KSPACE/pppm_tip4p.cpp             | 107 ++++++++++++-----------
    src/KSPACE/pppm_tip4p.h               |  18 +++-
    4 files changed, 207 insertions(+), 98 deletions(-)

Diff
----


:download:`download full diff <_diffs/iter_0026_4cbc0fa85603.diff>`

.. code-block:: diff

   diff --git a/src/KSPACE/pair_lj_cut_tip4p_long.cpp b/src/KSPACE/pair_lj_cut_tip4p_long.cpp
   index 7c4a239833..7d6d4a8fa2 100644
   --- a/src/KSPACE/pair_lj_cut_tip4p_long.cpp
   +++ b/src/KSPACE/pair_lj_cut_tip4p_long.cpp
   @@ -33,6 +33,7 @@
    
    #include <cmath>
    #include <cstring>
   +#include <limits>
    
    using namespace LAMMPS_NS;
    using namespace EwaldConst;
   @@ -56,6 +57,9 @@ PairLJCutTIP4PLong::PairLJCutTIP4PLong(LAMMPS *lmp) :
      nmax = 0;
      hneigh = nullptr;
      newsite = nullptr;
   +  neigh_stamp = 1;
   +  site_stamp = 1;
   +  tip4p_cache = nullptr;
    
      // TIP4P cannot compute virial as F dot r
      // due to finding bonded H atoms which are not near O atom
   @@ -69,6 +73,7 @@ PairLJCutTIP4PLong::~PairLJCutTIP4PLong()
    {
      memory->destroy(hneigh);
      memory->destroy(newsite);
   +  memory->destroy(tip4p_cache);
    }
    
    /* ---------------------------------------------------------------------- */
   @@ -77,23 +82,15 @@ void PairLJCutTIP4PLong::compute(int eflag, int vflag)
    {
      ev_init(eflag,vflag);
    
   -  // reallocate hneigh & newsite if necessary
   -  // initialize hneigh[0] to -1 on steps when reneighboring occurred
   -  // initialize hneigh[2] to 0 every step
   +  grow_tip4p_cache();
    
   -  const int nlocal = atom->nlocal;
   -  const int nall = nlocal + atom->nghost;
   -
   -  if (atom->nmax > nmax) {
   -    nmax = atom->nmax;
   -    memory->destroy(hneigh);
   -    memory->create(hneigh,nmax,3,"pair:hneigh");
   -    memory->destroy(newsite);
   -    memory->create(newsite,nmax,3,"pair:newsite");
   +  if (site_stamp == std::numeric_limits<int>::max()) reset_tip4p_cache_stamps();
   +  ++site_stamp;
   +
   +  if (neighbor->ago == 0) {
   +    if (neigh_stamp == std::numeric_limits<int>::max()) reset_tip4p_cache_stamps();
   +    ++neigh_stamp;
      }
   -  if (neighbor->ago == 0)
   -    for (int i = 0; i < nall; i++) hneigh[i][0] = -1;
   -  for (int i = 0; i < nall; i++) hneigh[i][2] = 0;
    
      if (!ncoultablebits) {
        if (evflag) {
   @@ -134,6 +131,7 @@ void PairLJCutTIP4PLong::eval()
      const int inum = list->inum;
      const int nlocal = atom->nlocal;
      const double cut_coulsqplus = (cut_coul+2.0*qdist) * (cut_coul+2.0*qdist);
   +  const bool pure_water = (atom->ntypes == 2);
    
      for (int ii = 0; ii < inum; ii++) {
        const int i = ilist[ii];
   @@ -154,21 +152,32 @@ void PairLJCutTIP4PLong::eval()
          const double cut_coul_plus_i_site =
            cut_coul + sqrt(site_delx*site_delx + site_dely*site_dely + site_delz*site_delz);
          const double cut_coulsq_i_h = cut_coul_plus_i_site * cut_coul_plus_i_site;
   -      eval_i<CTABLE,EVFLAG,EFLAG,VFLAG,1,0>(i,itype,qtmp,xtmp,ytmp,ztmp,iH1,iH2,x1,jlist,jnum,
   -                                            nlocal,cut_coulsq_i_h,cut_coulsqplus);
   +      if (pure_water)
   +        eval_i<CTABLE,EVFLAG,EFLAG,VFLAG,1,0,1>(i,itype,qtmp,xtmp,ytmp,ztmp,iH1,iH2,x1,jlist,
   +                                                jnum,nlocal,cut_coulsq_i_h,
   +                                                cut_coulsqplus);
   +      else
   +        eval_i<CTABLE,EVFLAG,EFLAG,VFLAG,1,0,0>(i,itype,qtmp,xtmp,ytmp,ztmp,iH1,iH2,x1,jlist,
   +                                                jnum,nlocal,cut_coulsq_i_h,
   +                                                cut_coulsqplus);
        } else if (itype == typeH) {
   -      eval_i<CTABLE,EVFLAG,EFLAG,VFLAG,0,1>(i,itype,qtmp,xtmp,ytmp,ztmp,-1,-1,x_d[i],jlist,jnum,
   -                                            nlocal,cut_coulsq,cut_coulsqplus);
   +      if (pure_water)
   +        eval_i<CTABLE,EVFLAG,EFLAG,VFLAG,0,1,1>(i,itype,qtmp,xtmp,ytmp,ztmp,-1,-1,x_d[i],jlist,
   +                                                jnum,nlocal,cut_coulsq,cut_coulsqplus);
   +      else
   +        eval_i<CTABLE,EVFLAG,EFLAG,VFLAG,0,1,0>(i,itype,qtmp,xtmp,ytmp,ztmp,-1,-1,x_d[i],jlist,
   +                                                jnum,nlocal,cut_coulsq,cut_coulsqplus);
        } else {
   -      eval_i<CTABLE,EVFLAG,EFLAG,VFLAG,0,0>(i,itype,qtmp,xtmp,ytmp,ztmp,-1,-1,x_d[i],jlist,jnum,
   -                                            nlocal,cut_coulsq,cut_coulsqplus);
   +      eval_i<CTABLE,EVFLAG,EFLAG,VFLAG,0,0,0>(i,itype,qtmp,xtmp,ytmp,ztmp,-1,-1,x_d[i],jlist,
   +                                              jnum,nlocal,cut_coulsq,cut_coulsqplus);
        }
      }
    }
    
    /* ---------------------------------------------------------------------- */
    
   -template <int CTABLE, int EVFLAG, int EFLAG, int VFLAG, int I_WATER, int I_HYDROGEN>
   +template <int CTABLE, int EVFLAG, int EFLAG, int VFLAG, int I_WATER, int I_HYDROGEN,
   +          int PURE_WATER>
    void PairLJCutTIP4PLong::eval_i(int i, int itype, double qtmp, double xtmp, double ytmp,
                                    double ztmp, int iH1, int iH2, double *x1, int *jlist,
                                    int jnum, int nlocal, double cut_coulsq_i_h,
   @@ -201,6 +210,12 @@ void PairLJCutTIP4PLong::eval_i(int i, int itype, double qtmp, double xtmp, doub
      double *lj3_i = nullptr;
      double *lj4_i = nullptr;
      double *offset_i = nullptr;
   +  double cut_ljsq_oo = 0.0;
   +  double lj1_oo = 0.0;
   +  double lj2_oo = 0.0;
   +  double lj3_oo = 0.0;
   +  double lj4_oo = 0.0;
   +  double offset_oo = 0.0;
    
      if constexpr (!I_HYDROGEN) {
        cut_ljsq_i = cut_ljsq[itype];
   @@ -209,12 +224,22 @@ void PairLJCutTIP4PLong::eval_i(int i, int itype, double qtmp, double xtmp, doub
        lj3_i = lj3[itype];
        lj4_i = lj4[itype];
        offset_i = offset[itype];
   +
   +    if constexpr (PURE_WATER) {
   +      cut_ljsq_oo = cut_ljsq_i[typeO];
   +      lj1_oo = lj1_i[typeO];
   +      lj2_oo = lj2_i[typeO];
   +      lj3_oo = lj3_i[typeO];
   +      lj4_oo = lj4_i[typeO];
   +      offset_oo = offset_i[typeO];
   +    }
      }
    
      for (jj = 0; jj < jnum; jj++) {
        j = jlist[jj];
        const int sb = sbmask(j);
   -    factor_coul = special_coul[sb];
   +    const bool special = (sb != 0);
   +    factor_coul = special ? special_coul[sb] : 1.0;
        j &= NEIGHMASK;
    
        delx = xtmp - x[j].x;
   @@ -223,11 +248,33 @@ void PairLJCutTIP4PLong::eval_i(int i, int itype, double qtmp, double xtmp, doub
        rsq = delx*delx + dely*dely + delz*delz;
        jtype = type[j];
        const bool j_is_oxygen = (jtype == typeO);
   -    const bool j_is_hydrogen = (jtype == typeH);
    
        if constexpr (!I_HYDROGEN) {
   -      if (!j_is_hydrogen && rsq < cut_ljsq_i[jtype]) {
   -        const double factor_lj = special_lj[sb];
   +      if constexpr (PURE_WATER) {
   +        if (j_is_oxygen && rsq < cut_ljsq_oo) {
   +          const double factor_lj = special ? special_lj[sb] : 1.0;
   +          r2inv = 1.0/rsq;
   +          r6inv = r2inv*r2inv*r2inv;
   +          forcelj = r6inv * (lj1_oo*r6inv - lj2_oo);
   +          forcelj *= factor_lj * r2inv;
   +
   +          f[i].x += delx*forcelj;
   +          f[i].y += dely*forcelj;
   +          f[i].z += delz*forcelj;
   +          f[j].x -= delx*forcelj;
   +          f[j].y -= dely*forcelj;
   +          f[j].z -= delz*forcelj;
   +
   +          if constexpr (EFLAG) {
   +            evdwl = r6inv*(lj3_oo*r6inv-lj4_oo) - offset_oo;
   +            evdwl *= factor_lj;
   +          } else evdwl = 0.0;
   +
   +          if constexpr (EVFLAG)
   +            ev_tally(i,j,nlocal,1,evdwl,0.0,forcelj,delx,dely,delz);
   +        }
   +      } else if (jtype != typeH && rsq < cut_ljsq_i[jtype]) {
   +        const double factor_lj = special ? special_lj[sb] : 1.0;
            r2inv = 1.0/rsq;
            r6inv = r2inv*r2inv*r2inv;
            forcelj = r6inv * (lj1_i[jtype]*r6inv - lj2_i[jtype]);
   @@ -250,8 +297,7 @@ void PairLJCutTIP4PLong::eval_i(int i, int itype, double qtmp, double xtmp, doub
          }
        }
    
   -    double cut_coulsq_precheck = cut_coulsqplus;
   -    if (!j_is_oxygen) cut_coulsq_precheck = cut_coulsq_i_h;
   +    const double cut_coulsq_precheck = j_is_oxygen ? cut_coulsqplus : cut_coulsq_i_h;
    
        if (rsq < cut_coulsq_precheck) {
          if constexpr (I_WATER) {
   @@ -439,10 +485,41 @@ void PairLJCutTIP4PLong::eval_i(int i, int itype, double qtmp, double xtmp, doub
    
    /* ---------------------------------------------------------------------- */
    
   +void PairLJCutTIP4PLong::grow_tip4p_cache()
   +{
   +  if (atom->nmax <= nmax) return;
   +
   +  nmax = atom->nmax;
   +  memory->destroy(tip4p_cache);
   +  memory->create(tip4p_cache,nmax,"pair:tip4p_cache");
   +
   +  for (int i = 0; i < nmax; i++) {
   +    tip4p_cache[i].neigh_stamp = 0;
   +    tip4p_cache[i].site_stamp = 0;
   +  }
   +}
   +
   +/* ---------------------------------------------------------------------- */
   +
   +void PairLJCutTIP4PLong::reset_tip4p_cache_stamps()
   +{
   +  for (int i = 0; i < nmax; i++) {
   +    tip4p_cache[i].neigh_stamp = 0;
   +    tip4p_cache[i].site_stamp = 0;
   +  }
   +
   +  neigh_stamp = 1;
   +  site_stamp = 1;
   +}
   +
   +/* ---------------------------------------------------------------------- */
   +
    double *PairLJCutTIP4PLong::tip4p_site(int i, int &iH1, int &iH2, double **x, tagint *tag,
                                           int *type)
    {
   -  if (hneigh[i][0] < 0) {
   +  Tip4pSiteCache &cache = tip4p_cache[i];
   +
   +  if (cache.neigh_stamp != neigh_stamp) {
        iH1 = atom->map(tag[i] + 1);
        iH2 = atom->map(tag[i] + 2);
        if (iH1 == -1 || iH2 == -1)
   @@ -452,21 +529,22 @@ double *PairLJCutTIP4PLong::tip4p_site(int i, int &iH1, int &iH2, double **x, ta
    
        iH1 = domain->closest_image(i,iH1);
        iH2 = domain->closest_image(i,iH2);
   -    compute_newsite(x[i],x[iH1],x[iH2],newsite[i]);
   -    hneigh[i][0] = iH1;
   -    hneigh[i][1] = iH2;
   -    hneigh[i][2] = 1;
   +    cache.iH1 = iH1;
   +    cache.iH2 = iH2;
   +    cache.neigh_stamp = neigh_stamp;
   +    cache.site_stamp = 0;
    
      } else {
   -    iH1 = hneigh[i][0];
   -    iH2 = hneigh[i][1];
   -    if (hneigh[i][2] == 0) {
   -      hneigh[i][2] = 1;
   -      compute_newsite(x[i],x[iH1],x[iH2],newsite[i]);
   -    }
   +    iH1 = cache.iH1;
   +    iH2 = cache.iH2;
   +  }
   +
   +  if (cache.site_stamp != site_stamp) {
   +    compute_newsite(x[i],x[iH1],x[iH2],cache.xM);
   +    cache.site_stamp = site_stamp;
      }
    
   -  return newsite[i];
   +  return cache.xM;
    }
    
    /* ----------------------------------------------------------------------
   @@ -679,6 +757,6 @@ double PairLJCutTIP4PLong::memory_usage()
    {
      double bytes = (double)maxeatom * sizeof(double);
      bytes += (double)maxvatom*6 * sizeof(double);
   -  bytes += (double)2 * nmax * sizeof(double);
   +  bytes += (double)nmax * sizeof(Tip4pSiteCache);
      return bytes;
    }
   diff --git a/src/KSPACE/pair_lj_cut_tip4p_long.h b/src/KSPACE/pair_lj_cut_tip4p_long.h
   index b97b0319cb..233fbb0dcc 100644
   --- a/src/KSPACE/pair_lj_cut_tip4p_long.h
   +++ b/src/KSPACE/pair_lj_cut_tip4p_long.h
   @@ -39,20 +39,32 @@ class PairLJCutTIP4PLong : public PairLJCutCoulLong {
      double memory_usage() override;
    
     protected:
   +  struct Tip4pSiteCache {
   +    int iH1;
   +    int iH2;
   +    int neigh_stamp;
   +    int site_stamp;
   +    double xM[3];
   +  };
   +
      std::string typeH_str, typeO_str, typeA_str, typeB_str;
      int typeH, typeO;    // atom types of TIP4P water H and O atoms
      int typeA, typeB;    // angle and bond types of TIP4P water
      double alpha;        // geometric constraint parameter for TIP4P
    
      int nmax;            // info on off-oxygen charge sites
   -  int **hneigh;        // 0,1 = indices of 2 H associated with O
   -                       // 2 = 0 if site loc not yet computed, 1 if yes
   -  double **newsite;    // locations of charge sites
   +  int **hneigh;        // legacy OPT cache storage
   +  double **newsite;    // legacy OPT cache storage
   +  int neigh_stamp;
   +  int site_stamp;
   +  Tip4pSiteCache *tip4p_cache;
    
      template <int, int, int, int> void eval();
   -  template <int, int, int, int, int, int>
   +  template <int, int, int, int, int, int, int>
      void eval_i(int, int, double, double, double, double, int, int, double *, int *, int, int,
                  double, double);
   +  void grow_tip4p_cache();
   +  void reset_tip4p_cache_stamps();
      double *tip4p_site(int, int &, int &, double **, tagint *, int *);
      void compute_newsite(double *, double *, double *, double *);
    };
   diff --git a/src/KSPACE/pppm_tip4p.cpp b/src/KSPACE/pppm_tip4p.cpp
   index ad04d3a89c..4ddee644a8 100644
   --- a/src/KSPACE/pppm_tip4p.cpp
   +++ b/src/KSPACE/pppm_tip4p.cpp
   @@ -27,6 +27,7 @@
    
    #include <cmath>
    #include <cstring>
   +#include <limits>
    
    using namespace LAMMPS_NS;
    using namespace MathConst;
   @@ -50,20 +51,21 @@ PPPMTIP4P::PPPMTIP4P(LAMMPS *lmp) : PPPM(lmp)
      tip4pflag = 1;
    
      tip4p_cache_nmax = 0;
   -  tip4p_h1 = nullptr;
   -  tip4p_h2 = nullptr;
   -  tip4p_cache_valid = nullptr;
   -  tip4p_xm = nullptr;
   +  tip4p_cache_stamp = 1;
   +  tip4p_oxygen_count = 0;
   +  tip4p_nonoxygen_count = 0;
   +  tip4p_oxygen = nullptr;
   +  tip4p_nonoxygen = nullptr;
   +  tip4p_cache = nullptr;
    }
    
    /* ---------------------------------------------------------------------- */
    
    PPPMTIP4P::~PPPMTIP4P()
    {
   -  memory->destroy(tip4p_h1);
   -  memory->destroy(tip4p_h2);
   -  memory->destroy(tip4p_cache_valid);
   -  memory->destroy(tip4p_xm);
   +  memory->destroy(tip4p_oxygen);
   +  memory->destroy(tip4p_nonoxygen);
   +  memory->destroy(tip4p_cache);
    }
    
    /* ---------------------------------------------------------------------- */
   @@ -86,14 +88,22 @@ void PPPMTIP4P::grow_tip4p_cache()
    
      tip4p_cache_nmax = atom->nmax;
    
   -  memory->destroy(tip4p_h1);
   -  memory->create(tip4p_h1,tip4p_cache_nmax,"pppm/tip4p:h1");
   -  memory->destroy(tip4p_h2);
   -  memory->create(tip4p_h2,tip4p_cache_nmax,"pppm/tip4p:h2");
   -  memory->destroy(tip4p_cache_valid);
   -  memory->create(tip4p_cache_valid,tip4p_cache_nmax,"pppm/tip4p:valid");
   -  memory->destroy(tip4p_xm);
   -  memory->create(tip4p_xm,tip4p_cache_nmax,3,"pppm/tip4p:xm");
   +  memory->destroy(tip4p_oxygen);
   +  memory->create(tip4p_oxygen,tip4p_cache_nmax,"pppm/tip4p:oxygen");
   +  memory->destroy(tip4p_nonoxygen);
   +  memory->create(tip4p_nonoxygen,tip4p_cache_nmax,"pppm/tip4p:nonoxygen");
   +  memory->destroy(tip4p_cache);
   +  memory->create(tip4p_cache,tip4p_cache_nmax,"pppm/tip4p:cache");
   +
   +  for (int i = 0; i < tip4p_cache_nmax; i++) tip4p_cache[i].stamp = 0;
   +}
   +
   +/* ---------------------------------------------------------------------- */
   +
   +void PPPMTIP4P::reset_tip4p_cache_stamps()
   +{
   +  for (int i = 0; i < tip4p_cache_nmax; i++) tip4p_cache[i].stamp = 0;
   +  tip4p_cache_stamp = 1;
    }
    
    /* ---------------------------------------------------------------------- */
   @@ -102,14 +112,15 @@ void PPPMTIP4P::find_M_cached(int i, int &iH1, int &iH2, double *&xM)
    {
      if (atom->nmax > tip4p_cache_nmax) grow_tip4p_cache();
    
   -  if (!tip4p_cache_valid[i]) {
   -    find_M(i,tip4p_h1[i],tip4p_h2[i],tip4p_xm[i]);
   -    tip4p_cache_valid[i] = 1;
   +  Tip4pCache &cache = tip4p_cache[i];
   +  if (cache.stamp != tip4p_cache_stamp) {
   +    find_M(i,cache.iH1,cache.iH2,cache.xM);
   +    cache.stamp = tip4p_cache_stamp;
      }
    
   -  iH1 = tip4p_h1[i];
   -  iH2 = tip4p_h2[i];
   -  xM = tip4p_xm[i];
   +  iH1 = cache.iH1;
   +  iH2 = cache.iH2;
   +  xM = cache.xM;
    }
    
    /* ----------------------------------------------------------------------
   @@ -131,10 +142,9 @@ void PPPMTIP4P::particle_map()
      const double boxloy = boxlo[1];
      const double boxloz = boxlo[2];
    
   -  if (nlocal > 0) {
   -    grow_tip4p_cache();
   -    std::memset(tip4p_cache_valid,0,nlocal*sizeof(int));
   -  }
   +  if (nlocal > 0) grow_tip4p_cache();
   +  if (tip4p_cache_stamp == std::numeric_limits<int>::max()) reset_tip4p_cache_stamps();
   +  ++tip4p_cache_stamp;
    
      if (!std::isfinite(boxlo[0]) || !std::isfinite(boxlo[1]) || !std::isfinite(boxlo[2]))
        error->one(FLERR,"Non-numeric box dimensions - simulation unstable" + utils::errorurl(6));
   @@ -186,9 +196,8 @@ void PPPMTIP4P::particle_map()
    
    void PPPMTIP4P::make_rho()
    {
   -  int l,m,n,nx,ny,nz,iH1,iH2;
   +  int i,l,m,n,nx,ny,nz;
      FFT_SCALAR dx,dy,dz,x0,y0,z0;
   -  double *xM;
    
      // clear 3d density array
    
   @@ -200,11 +209,9 @@ void PPPMTIP4P::make_rho()
      // (dx,dy,dz) = distance to "lower left" grid pt
      // (mx,my,mz) = global coords of moving stencil pt
    
   -  int *type = atom->type;
      double *q = atom->q;
      const auto * _noalias const x = (dbl3_t *) atom->x[0];
      const auto * _noalias const p2g = (int3_t *) part2grid[0];
   -  int nlocal = atom->nlocal;
      const int nx_out = nxhi_out - nxlo_out + 1;
      const int ny_out = nyhi_out - nylo_out + 1;
      const int nyz_out = nx_out * ny_out;
   @@ -212,10 +219,13 @@ void PPPMTIP4P::make_rho()
      const double boxloy = boxlo[1];
      const double boxloz = boxlo[2];
    
   -  for (int i = 0; i < nlocal; i++) {
   +  int *type = atom->type;
   +  const int nlocal = atom->nlocal;
   +
   +  for (i = 0; i < nlocal; i++) {
        double xix, xiy, xiz;
        if (type[i] == typeO) {
   -      find_M_cached(i,iH1,iH2,xM);
   +      const double *xM = tip4p_cache[i].xM;
          xix = xM[0];
          xiy = xM[1];
          xiz = xM[2];
   @@ -256,10 +266,9 @@ void PPPMTIP4P::make_rho()
    
    void PPPMTIP4P::fieldforce_ik()
    {
   -  int i,l,m,n,nx,ny,nz,iH1,iH2;
   +  int i,l,m,n,nx,ny,nz;
      FFT_SCALAR dx,dy,dz,x0,y0,z0;
      FFT_SCALAR ekx,eky,ekz;
   -  double *xM;
      double fx,fy,fz;
    
      // loop over my charges, interpolate electric field from nearby grid points
   @@ -278,8 +287,6 @@ void PPPMTIP4P::fieldforce_ik()
      FFT_SCALAR *vdy = &vdy_brick[nzlo_out][nylo_out][nxlo_out];
      FFT_SCALAR *vdz = &vdz_brick[nzlo_out][nylo_out][nxlo_out];
    
   -  int *type = atom->type;
   -  int nlocal = atom->nlocal;
      const int nx_out = nxhi_out - nxlo_out + 1;
      const int ny_out = nyhi_out - nylo_out + 1;
      const int nyz_out = nx_out * ny_out;
   @@ -287,13 +294,17 @@ void PPPMTIP4P::fieldforce_ik()
      const double boxloy = boxlo[1];
      const double boxloz = boxlo[2];
    
   +  int *type = atom->type;
   +  const int nlocal = atom->nlocal;
   +
      for (i = 0; i < nlocal; i++) {
   +    const Tip4pCache *cache = nullptr;
        double xix, xiy, xiz;
        if (type[i] == typeO) {
   -      find_M_cached(i,iH1,iH2,xM);
   -      xix = xM[0];
   -      xiy = xM[1];
   -      xiz = xM[2];
   +      cache = &tip4p_cache[i];
   +      xix = cache->xM[0];
   +      xiy = cache->xM[1];
   +      xiz = cache->xM[2];
        } else {
          xix = x[i].x;
          xiy = x[i].y;
   @@ -326,10 +337,8 @@ void PPPMTIP4P::fieldforce_ik()
          }
        }
    
   -    // convert E-field to force
   -
        const double qfactor = qqrd2e * scale * q[i];
   -    if (type[i] != typeO) {
   +    if (cache == nullptr) {
          f[i].x += qfactor*ekx;
          f[i].y += qfactor*eky;
          if (zforce_active) f[i].z += qfactor*ekz;
   @@ -343,13 +352,13 @@ void PPPMTIP4P::fieldforce_ik()
          f[i].y += fy * one_minus_alpha;
          if (zforce_active) f[i].z += fz * one_minus_alpha;
    
   -      f[iH1].x += half_alpha * fx;
   -      f[iH1].y += half_alpha * fy;
   -      if (zforce_active) f[iH1].z += half_alpha * fz;
   +      f[cache->iH1].x += half_alpha * fx;
   +      f[cache->iH1].y += half_alpha * fy;
   +      if (zforce_active) f[cache->iH1].z += half_alpha * fz;
    
   -      f[iH2].x += half_alpha * fx;
   -      f[iH2].y += half_alpha * fy;
   -      if (zforce_active) f[iH2].z += half_alpha * fz;
   +      f[cache->iH2].x += half_alpha * fx;
   +      f[cache->iH2].y += half_alpha * fy;
   +      if (zforce_active) f[cache->iH2].z += half_alpha * fz;
        }
      }
    }
   diff --git a/src/KSPACE/pppm_tip4p.h b/src/KSPACE/pppm_tip4p.h
   index abe3418959..d907f96dd3 100644
   --- a/src/KSPACE/pppm_tip4p.h
   +++ b/src/KSPACE/pppm_tip4p.h
   @@ -39,13 +39,23 @@ class PPPMTIP4P : public PPPM {
      void slabcorr() override;
    
     private:
   +  struct Tip4pCache {
   +    int iH1;
   +    int iH2;
   +    int stamp;
   +    double xM[3];
   +  };
   +
      int tip4p_cache_nmax;
   -  int *tip4p_h1;
   -  int *tip4p_h2;
   -  int *tip4p_cache_valid;
   -  double **tip4p_xm;
   +  int tip4p_cache_stamp;
   +  int tip4p_oxygen_count;
   +  int tip4p_nonoxygen_count;
   +  int *tip4p_oxygen;
   +  int *tip4p_nonoxygen;
   +  Tip4pCache *tip4p_cache;
    
      void grow_tip4p_cache();
   +  void reset_tip4p_cache_stamps();
      void find_M_cached(int, int &, int &, double *&);
      void find_M(int, int &, int &, double *);
    };


.. _iter-0026-page-head-4cbc0fa85603: https://github.com/skilled-scipkg/lammps/commit/4cbc0fa85603
.. _iter-0026-guardrails-incumbent-2599e1f9f4c8bb574886ac80f8e67052265a0ad4: https://github.com/skilled-scipkg/lammps/commit/2599e1f9f4c8bb574886ac80f8e67052265a0ad4
.. _iter-0026-guardrails-candidate-4cbc0fa856036c692930cd76775c671201e183cc: https://github.com/skilled-scipkg/lammps/commit/4cbc0fa856036c692930cd76775c671201e183cc