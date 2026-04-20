Iteration 0016 — ab55ea3e8508 (accepted)
========================================


GitHub commit: `ab55ea3e8508 <iter-0016-page-head-ab55ea3e8508_>`_
Published branch: `fermilink-optimize/lammps-tip4p <https://github.com/skilled-scipkg/lammps/tree/fermilink-optimize%2Flammps-tip4p>`_

Change summary
--------------


Restore safe hydrogen-neighbor LJ pruning in `pair_lj_cut_tip4p_long` and reuse cached TIP4P oxygen M/H state across PPPM `particle_map()`/`make_rho()`/`fieldforce_ik()` with `memset` density zeroing and linear IK brick indexing

Acceptance rationale
--------------------


Correctness passed and 78.12 beats incumbent 79.88 by 2.20%, clearing the 2% acceptance bar.

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
| incumbent commit | `2a0b902a8b5a <iter-0016-guardrails-incumbent-2a0b902a8b5a1a49c54d4236e27fa46ac09fb6ad_>`_ |
+------------------+--------------------------------------------------------------------------------------------+
| candidate commit | `ab55ea3e8508 <iter-0016-guardrails-candidate-ab55ea3e8508820bbebff62c9b733c8990e9a801_>`_ |
+------------------+--------------------------------------------------------------------------------------------+
| incumbent metric | 79.88                                                                                      |
+------------------+--------------------------------------------------------------------------------------------+
| candidate metric | 78.12                                                                                      |
+------------------+--------------------------------------------------------------------------------------------+
| baseline metric  | 94.717                                                                                     |
+------------------+--------------------------------------------------------------------------------------------+
| Δ vs incumbent   | +2.203% (lower-is-better sign)                                                             |
+------------------+--------------------------------------------------------------------------------------------+
| changed files    | src/KSPACE/pair_lj_cut_tip4p_long.cpp, src/KSPACE/pppm_tip4p.cpp, src/KSPACE/pppm_tip4p.h  |
+------------------+--------------------------------------------------------------------------------------------+


Diffstat
--------


.. code-block:: text

    src/KSPACE/pair_lj_cut_tip4p_long.cpp |  10 +--
    src/KSPACE/pppm_tip4p.cpp             | 131 +++++++++++++++++++++++++---------
    src/KSPACE/pppm_tip4p.h               |   9 +++
    3 files changed, 114 insertions(+), 36 deletions(-)

Diff
----


:download:`download full diff <_diffs/iter_0016_ab55ea3e8508.diff>`

.. code-block:: diff

   diff --git a/src/KSPACE/pair_lj_cut_tip4p_long.cpp b/src/KSPACE/pair_lj_cut_tip4p_long.cpp
   index c7b794571f..a070ca005d 100644
   --- a/src/KSPACE/pair_lj_cut_tip4p_long.cpp
   +++ b/src/KSPACE/pair_lj_cut_tip4p_long.cpp
   @@ -173,7 +173,7 @@ void PairLJCutTIP4PLong::eval_i(int i, int itype, double qtmp, double xtmp, doub
      int vlist[6];
      double delx, dely, delz, evdwl, ecoul, fraction, table;
      double r, rsq, r2inv, r6inv, forcecoul, forcelj, cforce;
   -  double factor_coul, factor_lj;
   +  double factor_coul;
      double grij, expm2, prefactor, t, erfc;
      double fO[3], fH[3], fd[3], v[6];
      double *x2, *xH1, *xH2;
   @@ -206,8 +206,8 @@ void PairLJCutTIP4PLong::eval_i(int i, int itype, double qtmp, double xtmp, doub
    
      for (jj = 0; jj < jnum; jj++) {
        j = jlist[jj];
   -    factor_lj = special_lj[sbmask(j)];
   -    factor_coul = special_coul[sbmask(j)];
   +    const int sb = sbmask(j);
   +    factor_coul = special_coul[sb];
        j &= NEIGHMASK;
    
        delx = xtmp - x[j][0];
   @@ -216,9 +216,11 @@ void PairLJCutTIP4PLong::eval_i(int i, int itype, double qtmp, double xtmp, doub
        rsq = delx*delx + dely*dely + delz*delz;
        jtype = type[j];
        const bool j_is_oxygen = (jtype == typeO);
   +    const bool j_is_hydrogen = (jtype == typeH);
    
        if constexpr (!I_HYDROGEN) {
   -      if (rsq < cut_ljsq_i[jtype]) {
   +      if (!j_is_hydrogen && rsq < cut_ljsq_i[jtype]) {
   +        const double factor_lj = special_lj[sb];
            r2inv = 1.0/rsq;
            r6inv = r2inv*r2inv*r2inv;
            forcelj = r6inv * (lj1_i[jtype]*r6inv - lj2_i[jtype]);
   diff --git a/src/KSPACE/pppm_tip4p.cpp b/src/KSPACE/pppm_tip4p.cpp
   index 631eebc65c..dae2a023f0 100644
   --- a/src/KSPACE/pppm_tip4p.cpp
   +++ b/src/KSPACE/pppm_tip4p.cpp
   @@ -23,8 +23,10 @@
    #include "force.h"
    #include "error.h"
    #include "math_const.h"
   +#include "memory.h"
    
    #include <cmath>
   +#include <cstring>
    
    using namespace LAMMPS_NS;
    using namespace MathConst;
   @@ -38,6 +40,22 @@ PPPMTIP4P::PPPMTIP4P(LAMMPS *lmp) : PPPM(lmp)
    {
      triclinic_support = 1;
      tip4pflag = 1;
   +
   +  tip4p_cache_nmax = 0;
   +  tip4p_h1 = nullptr;
   +  tip4p_h2 = nullptr;
   +  tip4p_cache_valid = nullptr;
   +  tip4p_xm = nullptr;
   +}
   +
   +/* ---------------------------------------------------------------------- */
   +
   +PPPMTIP4P::~PPPMTIP4P()
   +{
   +  memory->destroy(tip4p_h1);
   +  memory->destroy(tip4p_h2);
   +  memory->destroy(tip4p_cache_valid);
   +  memory->destroy(tip4p_xm);
    }
    
    /* ---------------------------------------------------------------------- */
   @@ -52,6 +70,40 @@ void PPPMTIP4P::init()
      PPPM::init();
    }
    
   +/* ---------------------------------------------------------------------- */
   +
   +void PPPMTIP4P::grow_tip4p_cache()
   +{
   +  if (atom->nmax <= tip4p_cache_nmax) return;
   +
   +  tip4p_cache_nmax = atom->nmax;
   +
   +  memory->destroy(tip4p_h1);
   +  memory->create(tip4p_h1,tip4p_cache_nmax,"pppm/tip4p:h1");
   +  memory->destroy(tip4p_h2);
   +  memory->create(tip4p_h2,tip4p_cache_nmax,"pppm/tip4p:h2");
   +  memory->destroy(tip4p_cache_valid);
   +  memory->create(tip4p_cache_valid,tip4p_cache_nmax,"pppm/tip4p:valid");
   +  memory->destroy(tip4p_xm);
   +  memory->create(tip4p_xm,tip4p_cache_nmax,3,"pppm/tip4p:xm");
   +}
   +
   +/* ---------------------------------------------------------------------- */
   +
   +void PPPMTIP4P::find_M_cached(int i, int &iH1, int &iH2, double *&xM)
   +{
   +  if (atom->nmax > tip4p_cache_nmax) grow_tip4p_cache();
   +
   +  if (!tip4p_cache_valid[i]) {
   +    find_M(i,tip4p_h1[i],tip4p_h2[i],tip4p_xm[i]);
   +    tip4p_cache_valid[i] = 1;
   +  }
   +
   +  iH1 = tip4p_h1[i];
   +  iH2 = tip4p_h2[i];
   +  xM = tip4p_xm[i];
   +}
   +
    /* ----------------------------------------------------------------------
       find center grid pt for each of my particles
       check that full stencil for the particle will fit in my 3d brick
   @@ -61,19 +113,24 @@ void PPPMTIP4P::init()
    void PPPMTIP4P::particle_map()
    {
      int nx,ny,nz,iH1,iH2;
   -  double *xi,xM[3];
   +  double *xi,*xM;
    
      int *type = atom->type;
      double **x = atom->x;
      int nlocal = atom->nlocal;
    
   +  if (nlocal > 0) {
   +    grow_tip4p_cache();
   +    std::memset(tip4p_cache_valid,0,nlocal*sizeof(int));
   +  }
   +
      if (!std::isfinite(boxlo[0]) || !std::isfinite(boxlo[1]) || !std::isfinite(boxlo[2]))
        error->one(FLERR,"Non-numeric box dimensions - simulation unstable" + utils::errorurl(6));
    
      int flag = 0;
      for (int i = 0; i < nlocal; i++) {
        if (type[i] == typeO) {
   -      find_M(i,iH1,iH2,xM);
   +      find_M_cached(i,iH1,iH2,xM);
          xi = xM;
        } else xi = x[i];
    
   @@ -110,14 +167,14 @@ void PPPMTIP4P::particle_map()
    
    void PPPMTIP4P::make_rho()
    {
   -  int i,l,m,n,nx,ny,nz,mx,my,mz,iH1,iH2;
   +  int l,m,n,nx,ny,nz,iH1,iH2;
      FFT_SCALAR dx,dy,dz,x0,y0,z0;
   -  double *xi,xM[3];
   +  double *xi,*xM;
    
      // clear 3d density array
    
   -  FFT_SCALAR *vec = &density_brick[nzlo_out][nylo_out][nxlo_out];
   -  for (i = 0; i < ngrid; i++) vec[i] = ZEROF;
   +  FFT_SCALAR *density = &density_brick[nzlo_out][nylo_out][nxlo_out];
   +  std::memset(density,0,ngrid*sizeof(FFT_SCALAR));
    
      // loop over my charges, add their contribution to nearby grid points
      // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
   @@ -128,10 +185,13 @@ void PPPMTIP4P::make_rho()
      double *q = atom->q;
      double **x = atom->x;
      int nlocal = atom->nlocal;
   +  const int nx_out = nxhi_out - nxlo_out + 1;
   +  const int ny_out = nyhi_out - nylo_out + 1;
   +  const int nyz_out = nx_out * ny_out;
    
      for (int i = 0; i < nlocal; i++) {
        if (type[i] == typeO) {
   -      find_M(i,iH1,iH2,xM);
   +      find_M_cached(i,iH1,iH2,xM);
          xi = xM;
        } else xi = x[i];
    
   @@ -146,15 +206,15 @@ void PPPMTIP4P::make_rho()
    
        z0 = delvolinv * q[i];
        for (n = nlower; n <= nupper; n++) {
   -      mz = n+nz;
   +      int index = ((nz+n-nzlo_out) * nyz_out) + (ny+nlower-nylo_out) * nx_out +
   +        (nx+nlower-nxlo_out);
          y0 = z0*rho1d[2][n];
          for (m = nlower; m <= nupper; m++) {
   -        my = m+ny;
            x0 = y0*rho1d[1][m];
            for (l = nlower; l <= nupper; l++) {
   -          mx = l+nx;
   -          density_brick[mz][my][mx] += x0*rho1d[0][l];
   +          density[index+l-nlower] += x0*rho1d[0][l];
            }
   +        index += nx_out;
          }
        }
      }
   @@ -166,12 +226,11 @@ void PPPMTIP4P::make_rho()
    
    void PPPMTIP4P::fieldforce_ik()
    {
   -  int i,l,m,n,nx,ny,nz,mx,my,mz;
   +  int i,l,m,n,nx,ny,nz,iH1,iH2;
      FFT_SCALAR dx,dy,dz,x0,y0,z0;
      FFT_SCALAR ekx,eky,ekz;
      double *xi;
   -  int iH1,iH2;
   -  double xM[3];
   +  double *xM;
      double fx,fy,fz;
    
      // loop over my charges, interpolate electric field from nearby grid points
   @@ -182,13 +241,22 @@ void PPPMTIP4P::fieldforce_ik()
      double *q = atom->q;
      double **x = atom->x;
      double **f = atom->f;
   +  const double half_alpha = 0.5 * alpha;
   +  const double one_minus_alpha = 1.0 - alpha;
   +  const int zforce_active = (slabflag != 2);
   +  FFT_SCALAR *vdx = &vdx_brick[nzlo_out][nylo_out][nxlo_out];
   +  FFT_SCALAR *vdy = &vdy_brick[nzlo_out][nylo_out][nxlo_out];
   +  FFT_SCALAR *vdz = &vdz_brick[nzlo_out][nylo_out][nxlo_out];
    
      int *type = atom->type;
      int nlocal = atom->nlocal;
   +  const int nx_out = nxhi_out - nxlo_out + 1;
   +  const int ny_out = nyhi_out - nylo_out + 1;
   +  const int nyz_out = nx_out * ny_out;
    
      for (i = 0; i < nlocal; i++) {
        if (type[i] == typeO) {
   -      find_M(i,iH1,iH2,xM);
   +      find_M_cached(i,iH1,iH2,xM);
          xi = xM;
        } else xi = x[i];
    
   @@ -203,18 +271,18 @@ void PPPMTIP4P::fieldforce_ik()
    
        ekx = eky = ekz = ZEROF;
        for (n = nlower; n <= nupper; n++) {
   -      mz = n+nz;
   +      int index = ((nz+n-nzlo_out) * nyz_out) + (ny+nlower-nylo_out) * nx_out +
   +        (nx+nlower-nxlo_out);
          z0 = rho1d[2][n];
          for (m = nlower; m <= nupper; m++) {
   -        my = m+ny;
            y0 = z0*rho1d[1][m];
            for (l = nlower; l <= nupper; l++) {
   -          mx = l+nx;
              x0 = y0*rho1d[0][l];
   -          ekx -= x0*vdx_brick[mz][my][mx];
   -          eky -= x0*vdy_brick[mz][my][mx];
   -          ekz -= x0*vdz_brick[mz][my][mx];
   +          ekx -= x0*vdx[index+l-nlower];
   +          eky -= x0*vdy[index+l-nlower];
   +          ekz -= x0*vdz[index+l-nlower];
            }
   +        index += nx_out;
          }
        }
    
   @@ -224,25 +292,24 @@ void PPPMTIP4P::fieldforce_ik()
        if (type[i] != typeO) {
          f[i][0] += qfactor*ekx;
          f[i][1] += qfactor*eky;
   -      if (slabflag != 2) f[i][2] += qfactor*ekz;
   +      if (zforce_active) f[i][2] += qfactor*ekz;
    
        } else {
          fx = qfactor * ekx;
          fy = qfactor * eky;
          fz = qfactor * ekz;
   -      find_M(i,iH1,iH2,xM);
    
   -      f[i][0] += fx*(1 - alpha);
   -      f[i][1] += fy*(1 - alpha);
   -      if (slabflag != 2) f[i][2] += fz*(1 - alpha);
   +      f[i][0] += fx * one_minus_alpha;
   +      f[i][1] += fy * one_minus_alpha;
   +      if (zforce_active) f[i][2] += fz * one_minus_alpha;
    
   -      f[iH1][0] += 0.5*alpha*fx;
   -      f[iH1][1] += 0.5*alpha*fy;
   -      if (slabflag != 2) f[iH1][2] += 0.5*alpha*fz;
   +      f[iH1][0] += half_alpha * fx;
   +      f[iH1][1] += half_alpha * fy;
   +      if (zforce_active) f[iH1][2] += half_alpha * fz;
    
   -      f[iH2][0] += 0.5*alpha*fx;
   -      f[iH2][1] += 0.5*alpha*fy;
   -      if (slabflag != 2) f[iH2][2] += 0.5*alpha*fz;
   +      f[iH2][0] += half_alpha * fx;
   +      f[iH2][1] += half_alpha * fy;
   +      if (zforce_active) f[iH2][2] += half_alpha * fz;
        }
      }
    }
   diff --git a/src/KSPACE/pppm_tip4p.h b/src/KSPACE/pppm_tip4p.h
   index 6af532c409..abe3418959 100644
   --- a/src/KSPACE/pppm_tip4p.h
   +++ b/src/KSPACE/pppm_tip4p.h
   @@ -27,6 +27,7 @@ namespace LAMMPS_NS {
    class PPPMTIP4P : public PPPM {
     public:
      PPPMTIP4P(class LAMMPS *);
   +  ~PPPMTIP4P() override;
      void init() override;
    
     protected:
   @@ -38,6 +39,14 @@ class PPPMTIP4P : public PPPM {
      void slabcorr() override;
    
     private:
   +  int tip4p_cache_nmax;
   +  int *tip4p_h1;
   +  int *tip4p_h2;
   +  int *tip4p_cache_valid;
   +  double **tip4p_xm;
   +
   +  void grow_tip4p_cache();
   +  void find_M_cached(int, int &, int &, double *&);
      void find_M(int, int &, int &, double *);
    };
    


.. _iter-0016-page-head-ab55ea3e8508: https://github.com/skilled-scipkg/lammps/commit/ab55ea3e8508
.. _iter-0016-guardrails-incumbent-2a0b902a8b5a1a49c54d4236e27fa46ac09fb6ad: https://github.com/skilled-scipkg/lammps/commit/2a0b902a8b5a1a49c54d4236e27fa46ac09fb6ad
.. _iter-0016-guardrails-candidate-ab55ea3e8508820bbebff62c9b733c8990e9a801: https://github.com/skilled-scipkg/lammps/commit/ab55ea3e8508820bbebff62c9b733c8990e9a801