Iteration 0021 — 2599e1f9f4c8 (accepted)
========================================


GitHub commit: `2599e1f9f4c8 <iter-0021-page-head-2599e1f9f4c8_>`_
Published branch: `fermilink-optimize/lammps-tip4p <https://github.com/skilled-scipkg/lammps/tree/fermilink-optimize%2Flammps-tip4p>`_

Change summary
--------------


Use contiguous alias-friendly x/f/part2grid views in pair_lj_cut_tip4p_long and TIP4P PPPM particle_map/make_rho/fieldforce_ik without changing force, tally, or cache semantics

Acceptance rationale
--------------------


Correctness passed and candidate 2599e1f9f4c8 reduces weighted_median_pair_plus_kspace_seconds from 78.12 to 76.248, a 2.40% improvement over incumbent ab55ea3e8508 that clears the 2% acceptance bar.

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
| incumbent commit | `ab55ea3e8508 <iter-0021-guardrails-incumbent-ab55ea3e8508820bbebff62c9b733c8990e9a801_>`_ |
+------------------+--------------------------------------------------------------------------------------------+
| candidate commit | `2599e1f9f4c8 <iter-0021-guardrails-candidate-2599e1f9f4c8bb574886ac80f8e67052265a0ad4_>`_ |
+------------------+--------------------------------------------------------------------------------------------+
| incumbent metric | 78.12                                                                                      |
+------------------+--------------------------------------------------------------------------------------------+
| candidate metric | 76.248                                                                                     |
+------------------+--------------------------------------------------------------------------------------------+
| baseline metric  | 94.717                                                                                     |
+------------------+--------------------------------------------------------------------------------------------+
| Δ vs incumbent   | +2.396% (lower-is-better sign)                                                             |
+------------------+--------------------------------------------------------------------------------------------+
| changed files    | src/KSPACE/pair_lj_cut_tip4p_long.cpp, src/KSPACE/pppm_tip4p.cpp                           |
+------------------+--------------------------------------------------------------------------------------------+


Diffstat
--------


.. code-block:: text

    src/KSPACE/pair_lj_cut_tip4p_long.cpp | 143 ++++++++++++++++++----------------
    src/KSPACE/pppm_tip4p.cpp             | 128 +++++++++++++++++++-----------
    2 files changed, 159 insertions(+), 112 deletions(-)

Diff
----


:download:`download full diff <_diffs/iter_0021_2599e1f9f4c8.diff>`

.. code-block:: diff

   diff --git a/src/KSPACE/pair_lj_cut_tip4p_long.cpp b/src/KSPACE/pair_lj_cut_tip4p_long.cpp
   index a070ca005d..7c4a239833 100644
   --- a/src/KSPACE/pair_lj_cut_tip4p_long.cpp
   +++ b/src/KSPACE/pair_lj_cut_tip4p_long.cpp
   @@ -37,6 +37,10 @@
    using namespace LAMMPS_NS;
    using namespace EwaldConst;
    
   +using dbl3_t = struct {
   +  double x, y, z;
   +};
   +
    /* ---------------------------------------------------------------------- */
    
    PairLJCutTIP4PLong::PairLJCutTIP4PLong(LAMMPS *lmp) :
   @@ -119,7 +123,8 @@ void PairLJCutTIP4PLong::compute(int eflag, int vflag)
    template <int CTABLE, int EVFLAG, int EFLAG, int VFLAG>
    void PairLJCutTIP4PLong::eval()
    {
   -  double **x = atom->x;
   +  const auto * _noalias const x = (dbl3_t *) atom->x[0];
   +  double **x_d = atom->x;
      double *q = atom->q;
      tagint *tag = atom->tag;
      int *type = atom->type;
   @@ -133,16 +138,16 @@ void PairLJCutTIP4PLong::eval()
      for (int ii = 0; ii < inum; ii++) {
        const int i = ilist[ii];
        const double qtmp = q[i];
   -    const double xtmp = x[i][0];
   -    const double ytmp = x[i][1];
   -    const double ztmp = x[i][2];
   +    const double xtmp = x[i].x;
   +    const double ytmp = x[i].y;
   +    const double ztmp = x[i].z;
        const int itype = type[i];
        int *jlist = firstneigh[i];
        const int jnum = numneigh[i];
    
        if (itype == typeO) {
          int iH1, iH2;
   -      double *x1 = tip4p_site(i,iH1,iH2,x,tag,type);
   +      double *x1 = tip4p_site(i,iH1,iH2,x_d,tag,type);
          const double site_delx = x1[0] - xtmp;
          const double site_dely = x1[1] - ytmp;
          const double site_delz = x1[2] - ztmp;
   @@ -152,10 +157,10 @@ void PairLJCutTIP4PLong::eval()
          eval_i<CTABLE,EVFLAG,EFLAG,VFLAG,1,0>(i,itype,qtmp,xtmp,ytmp,ztmp,iH1,iH2,x1,jlist,jnum,
                                                nlocal,cut_coulsq_i_h,cut_coulsqplus);
        } else if (itype == typeH) {
   -      eval_i<CTABLE,EVFLAG,EFLAG,VFLAG,0,1>(i,itype,qtmp,xtmp,ytmp,ztmp,-1,-1,x[i],jlist,jnum,
   +      eval_i<CTABLE,EVFLAG,EFLAG,VFLAG,0,1>(i,itype,qtmp,xtmp,ytmp,ztmp,-1,-1,x_d[i],jlist,jnum,
                                                nlocal,cut_coulsq,cut_coulsqplus);
        } else {
   -      eval_i<CTABLE,EVFLAG,EFLAG,VFLAG,0,0>(i,itype,qtmp,xtmp,ytmp,ztmp,-1,-1,x[i],jlist,jnum,
   +      eval_i<CTABLE,EVFLAG,EFLAG,VFLAG,0,0>(i,itype,qtmp,xtmp,ytmp,ztmp,-1,-1,x_d[i],jlist,jnum,
                                                nlocal,cut_coulsq,cut_coulsqplus);
        }
      }
   @@ -176,10 +181,12 @@ void PairLJCutTIP4PLong::eval_i(int i, int itype, double qtmp, double xtmp, doub
      double factor_coul;
      double grij, expm2, prefactor, t, erfc;
      double fO[3], fH[3], fd[3], v[6];
   -  double *x2, *xH1, *xH2;
   +  double *x2;
   +  const dbl3_t *xH1, *xH2;
    
   -  double **f = atom->f;
   -  double **x = atom->x;
   +  auto * _noalias const f = (dbl3_t *) atom->f[0];
   +  const auto * _noalias const x = (dbl3_t *) atom->x[0];
   +  double **x_d = atom->x;
      double *q = atom->q;
      tagint *tag = atom->tag;
      int *type = atom->type;
   @@ -210,9 +217,9 @@ void PairLJCutTIP4PLong::eval_i(int i, int itype, double qtmp, double xtmp, doub
        factor_coul = special_coul[sb];
        j &= NEIGHMASK;
    
   -    delx = xtmp - x[j][0];
   -    dely = ytmp - x[j][1];
   -    delz = ztmp - x[j][2];
   +    delx = xtmp - x[j].x;
   +    dely = ytmp - x[j].y;
   +    delz = ztmp - x[j].z;
        rsq = delx*delx + dely*dely + delz*delz;
        jtype = type[j];
        const bool j_is_oxygen = (jtype == typeO);
   @@ -226,12 +233,12 @@ void PairLJCutTIP4PLong::eval_i(int i, int itype, double qtmp, double xtmp, doub
            forcelj = r6inv * (lj1_i[jtype]*r6inv - lj2_i[jtype]);
            forcelj *= factor_lj * r2inv;
    
   -        f[i][0] += delx*forcelj;
   -        f[i][1] += dely*forcelj;
   -        f[i][2] += delz*forcelj;
   -        f[j][0] -= delx*forcelj;
   -        f[j][1] -= dely*forcelj;
   -        f[j][2] -= delz*forcelj;
   +        f[i].x += delx*forcelj;
   +        f[i].y += dely*forcelj;
   +        f[i].z += delz*forcelj;
   +        f[j].x -= delx*forcelj;
   +        f[j].y -= dely*forcelj;
   +        f[j].z -= delz*forcelj;
    
            if constexpr (EFLAG) {
              evdwl = r6inv*(lj3_i[jtype]*r6inv-lj4_i[jtype]) - offset_i[jtype];
   @@ -248,15 +255,15 @@ void PairLJCutTIP4PLong::eval_i(int i, int itype, double qtmp, double xtmp, doub
    
        if (rsq < cut_coulsq_precheck) {
          if constexpr (I_WATER) {
   -        if (j_is_oxygen) x2 = tip4p_site(j,jH1,jH2,x,tag,type);
   -        else x2 = x[j];
   +        if (j_is_oxygen) x2 = tip4p_site(j,jH1,jH2,x_d,tag,type);
   +        else x2 = x_d[j];
    
            delx = x1[0] - x2[0];
            dely = x1[1] - x2[1];
            delz = x1[2] - x2[2];
            rsq = delx*delx + dely*dely + delz*delz;
          } else if (j_is_oxygen) {
   -        x2 = tip4p_site(j,jH1,jH2,x,tag,type);
   +        x2 = tip4p_site(j,jH1,jH2,x_d,tag,type);
            delx = x1[0] - x2[0];
            dely = x1[1] - x2[1];
            delz = x1[2] - x2[2];
   @@ -299,9 +306,9 @@ void PairLJCutTIP4PLong::eval_i(int i, int itype, double qtmp, double xtmp, doub
            }
    
            if constexpr (!I_WATER) {
   -          f[i][0] += delx * cforce;
   -          f[i][1] += dely * cforce;
   -          f[i][2] += delz * cforce;
   +          f[i].x += delx * cforce;
   +          f[i].y += dely * cforce;
   +          f[i].z += delz * cforce;
    
              if constexpr (VFLAG) {
                v[0] = xtmp * delx * cforce;
   @@ -326,27 +333,27 @@ void PairLJCutTIP4PLong::eval_i(int i, int itype, double qtmp, double xtmp, doub
              fH[1] = half_alpha * fd[1];
              fH[2] = half_alpha * fd[2];
    
   -          f[i][0] += fO[0];
   -          f[i][1] += fO[1];
   -          f[i][2] += fO[2];
   +          f[i].x += fO[0];
   +          f[i].y += fO[1];
   +          f[i].z += fO[2];
    
   -          f[iH1][0] += fH[0];
   -          f[iH1][1] += fH[1];
   -          f[iH1][2] += fH[2];
   +          f[iH1].x += fH[0];
   +          f[iH1].y += fH[1];
   +          f[iH1].z += fH[2];
    
   -          f[iH2][0] += fH[0];
   -          f[iH2][1] += fH[1];
   -          f[iH2][2] += fH[2];
   +          f[iH2].x += fH[0];
   +          f[iH2].y += fH[1];
   +          f[iH2].z += fH[2];
    
              if constexpr (VFLAG) {
   -            xH1 = x[iH1];
   -            xH2 = x[iH2];
   -            v[0] = xtmp*fO[0] + xH1[0]*fH[0] + xH2[0]*fH[0];
   -            v[1] = ytmp*fO[1] + xH1[1]*fH[1] + xH2[1]*fH[1];
   -            v[2] = ztmp*fO[2] + xH1[2]*fH[2] + xH2[2]*fH[2];
   -            v[3] = xtmp*fO[1] + xH1[0]*fH[1] + xH2[0]*fH[1];
   -            v[4] = xtmp*fO[2] + xH1[0]*fH[2] + xH2[0]*fH[2];
   -            v[5] = ytmp*fO[2] + xH1[1]*fH[2] + xH2[1]*fH[2];
   +            xH1 = &x[iH1];
   +            xH2 = &x[iH2];
   +            v[0] = xtmp*fO[0] + xH1->x*fH[0] + xH2->x*fH[0];
   +            v[1] = ytmp*fO[1] + xH1->y*fH[1] + xH2->y*fH[1];
   +            v[2] = ztmp*fO[2] + xH1->z*fH[2] + xH2->z*fH[2];
   +            v[3] = xtmp*fO[1] + xH1->x*fH[1] + xH2->x*fH[1];
   +            v[4] = xtmp*fO[2] + xH1->x*fH[2] + xH2->x*fH[2];
   +            v[5] = ytmp*fO[2] + xH1->y*fH[2] + xH2->y*fH[2];
              }
              if constexpr (EVFLAG) {
                vlist[n++] = i;
   @@ -356,17 +363,17 @@ void PairLJCutTIP4PLong::eval_i(int i, int itype, double qtmp, double xtmp, doub
            }
    
            if (!j_is_oxygen) {
   -          f[j][0] -= delx * cforce;
   -          f[j][1] -= dely * cforce;
   -          f[j][2] -= delz * cforce;
   +          f[j].x -= delx * cforce;
   +          f[j].y -= dely * cforce;
   +          f[j].z -= delz * cforce;
    
              if constexpr (VFLAG) {
   -            v[0] -= x[j][0] * delx * cforce;
   -            v[1] -= x[j][1] * dely * cforce;
   -            v[2] -= x[j][2] * delz * cforce;
   -            v[3] -= x[j][0] * dely * cforce;
   -            v[4] -= x[j][0] * delz * cforce;
   -            v[5] -= x[j][1] * delz * cforce;
   +            v[0] -= x[j].x * delx * cforce;
   +            v[1] -= x[j].y * dely * cforce;
   +            v[2] -= x[j].z * delz * cforce;
   +            v[3] -= x[j].x * dely * cforce;
   +            v[4] -= x[j].x * delz * cforce;
   +            v[5] -= x[j].y * delz * cforce;
              }
              if constexpr (EVFLAG) vlist[n++] = j;
    
   @@ -385,27 +392,27 @@ void PairLJCutTIP4PLong::eval_i(int i, int itype, double qtmp, double xtmp, doub
              fH[1] = half_alpha * fd[1];
              fH[2] = half_alpha * fd[2];
    
   -          f[j][0] += fO[0];
   -          f[j][1] += fO[1];
   -          f[j][2] += fO[2];
   +          f[j].x += fO[0];
   +          f[j].y += fO[1];
   +          f[j].z += fO[2];
    
   -          f[jH1][0] += fH[0];
   -          f[jH1][1] += fH[1];
   -          f[jH1][2] += fH[2];
   +          f[jH1].x += fH[0];
   +          f[jH1].y += fH[1];
   +          f[jH1].z += fH[2];
    
   -          f[jH2][0] += fH[0];
   -          f[jH2][1] += fH[1];
   -          f[jH2][2] += fH[2];
   +          f[jH2].x += fH[0];
   +          f[jH2].y += fH[1];
   +          f[jH2].z += fH[2];
    
              if constexpr (VFLAG) {
   -            xH1 = x[jH1];
   -            xH2 = x[jH2];
   -            v[0] += x[j][0]*fO[0] + xH1[0]*fH[0] + xH2[0]*fH[0];
   -            v[1] += x[j][1]*fO[1] + xH1[1]*fH[1] + xH2[1]*fH[1];
   -            v[2] += x[j][2]*fO[2] + xH1[2]*fH[2] + xH2[2]*fH[2];
   -            v[3] += x[j][0]*fO[1] + xH1[0]*fH[1] + xH2[0]*fH[1];
   -            v[4] += x[j][0]*fO[2] + xH1[0]*fH[2] + xH2[0]*fH[2];
   -            v[5] += x[j][1]*fO[2] + xH1[1]*fH[2] + xH2[1]*fH[2];
   +            xH1 = &x[jH1];
   +            xH2 = &x[jH2];
   +            v[0] += x[j].x*fO[0] + xH1->x*fH[0] + xH2->x*fH[0];
   +            v[1] += x[j].y*fO[1] + xH1->y*fH[1] + xH2->y*fH[1];
   +            v[2] += x[j].z*fO[2] + xH1->z*fH[2] + xH2->z*fH[2];
   +            v[3] += x[j].x*fO[1] + xH1->x*fH[1] + xH2->x*fH[1];
   +            v[4] += x[j].x*fO[2] + xH1->x*fH[2] + xH2->x*fH[2];
   +            v[5] += x[j].y*fO[2] + xH1->y*fH[2] + xH2->y*fH[2];
              }
              if constexpr (EVFLAG) {
                vlist[n++] = j;
   diff --git a/src/KSPACE/pppm_tip4p.cpp b/src/KSPACE/pppm_tip4p.cpp
   index dae2a023f0..ad04d3a89c 100644
   --- a/src/KSPACE/pppm_tip4p.cpp
   +++ b/src/KSPACE/pppm_tip4p.cpp
   @@ -34,6 +34,14 @@ using namespace MathConst;
    static constexpr FFT_SCALAR ZEROF = 0.0;
    static constexpr int OFFSET = 16384;
    
   +using dbl3_t = struct {
   +  double x, y, z;
   +};
   +
   +using int3_t = struct {
   +  int a, b, t;
   +};
   +
    /* ---------------------------------------------------------------------- */
    
    PPPMTIP4P::PPPMTIP4P(LAMMPS *lmp) : PPPM(lmp)
   @@ -113,11 +121,15 @@ void PPPMTIP4P::find_M_cached(int i, int &iH1, int &iH2, double *&xM)
    void PPPMTIP4P::particle_map()
    {
      int nx,ny,nz,iH1,iH2;
   -  double *xi,*xM;
   +  double *xM;
    
      int *type = atom->type;
   -  double **x = atom->x;
   +  const auto * _noalias const x = (dbl3_t *) atom->x[0];
   +  auto * _noalias const p2g = (int3_t *) part2grid[0];
      int nlocal = atom->nlocal;
   +  const double boxlox = boxlo[0];
   +  const double boxloy = boxlo[1];
   +  const double boxloz = boxlo[2];
    
      if (nlocal > 0) {
        grow_tip4p_cache();
   @@ -129,22 +141,29 @@ void PPPMTIP4P::particle_map()
    
      int flag = 0;
      for (int i = 0; i < nlocal; i++) {
   +    double xix, xiy, xiz;
        if (type[i] == typeO) {
          find_M_cached(i,iH1,iH2,xM);
   -      xi = xM;
   -    } else xi = x[i];
   +      xix = xM[0];
   +      xiy = xM[1];
   +      xiz = xM[2];
   +    } else {
   +      xix = x[i].x;
   +      xiy = x[i].y;
   +      xiz = x[i].z;
   +    }
    
        // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
        // current particle coord can be outside global and local box
        // add/subtract OFFSET to avoid int(-0.75) = 0 when want it to be -1
    
   -    nx = static_cast<int> ((xi[0]-boxlo[0])*delxinv+shift) - OFFSET;
   -    ny = static_cast<int> ((xi[1]-boxlo[1])*delyinv+shift) - OFFSET;
   -    nz = static_cast<int> ((xi[2]-boxlo[2])*delzinv+shift) - OFFSET;
   +    nx = static_cast<int> ((xix-boxlox)*delxinv+shift) - OFFSET;
   +    ny = static_cast<int> ((xiy-boxloy)*delyinv+shift) - OFFSET;
   +    nz = static_cast<int> ((xiz-boxloz)*delzinv+shift) - OFFSET;
    
   -    part2grid[i][0] = nx;
   -    part2grid[i][1] = ny;
   -    part2grid[i][2] = nz;
   +    p2g[i].a = nx;
   +    p2g[i].b = ny;
   +    p2g[i].t = nz;
    
        // check that entire stencil around nx,ny,nz will fit in my 3d brick
    
   @@ -169,11 +188,11 @@ void PPPMTIP4P::make_rho()
    {
      int l,m,n,nx,ny,nz,iH1,iH2;
      FFT_SCALAR dx,dy,dz,x0,y0,z0;
   -  double *xi,*xM;
   +  double *xM;
    
      // clear 3d density array
    
   -  FFT_SCALAR *density = &density_brick[nzlo_out][nylo_out][nxlo_out];
   +  FFT_SCALAR * _noalias const density = &density_brick[nzlo_out][nylo_out][nxlo_out];
      std::memset(density,0,ngrid*sizeof(FFT_SCALAR));
    
      // loop over my charges, add their contribution to nearby grid points
   @@ -183,24 +202,35 @@ void PPPMTIP4P::make_rho()
    
      int *type = atom->type;
      double *q = atom->q;
   -  double **x = atom->x;
   +  const auto * _noalias const x = (dbl3_t *) atom->x[0];
   +  const auto * _noalias const p2g = (int3_t *) part2grid[0];
      int nlocal = atom->nlocal;
      const int nx_out = nxhi_out - nxlo_out + 1;
      const int ny_out = nyhi_out - nylo_out + 1;
      const int nyz_out = nx_out * ny_out;
   +  const double boxlox = boxlo[0];
   +  const double boxloy = boxlo[1];
   +  const double boxloz = boxlo[2];
    
      for (int i = 0; i < nlocal; i++) {
   +    double xix, xiy, xiz;
        if (type[i] == typeO) {
          find_M_cached(i,iH1,iH2,xM);
   -      xi = xM;
   -    } else xi = x[i];
   +      xix = xM[0];
   +      xiy = xM[1];
   +      xiz = xM[2];
   +    } else {
   +      xix = x[i].x;
   +      xiy = x[i].y;
   +      xiz = x[i].z;
   +    }
    
   -    nx = part2grid[i][0];
   -    ny = part2grid[i][1];
   -    nz = part2grid[i][2];
   -    dx = nx+shiftone - (xi[0]-boxlo[0])*delxinv;
   -    dy = ny+shiftone - (xi[1]-boxlo[1])*delyinv;
   -    dz = nz+shiftone - (xi[2]-boxlo[2])*delzinv;
   +    nx = p2g[i].a;
   +    ny = p2g[i].b;
   +    nz = p2g[i].t;
   +    dx = nx+shiftone - (xix-boxlox)*delxinv;
   +    dy = ny+shiftone - (xiy-boxloy)*delyinv;
   +    dz = nz+shiftone - (xiz-boxloz)*delzinv;
    
        compute_rho1d(dx,dy,dz);
    
   @@ -229,7 +259,6 @@ void PPPMTIP4P::fieldforce_ik()
      int i,l,m,n,nx,ny,nz,iH1,iH2;
      FFT_SCALAR dx,dy,dz,x0,y0,z0;
      FFT_SCALAR ekx,eky,ekz;
   -  double *xi;
      double *xM;
      double fx,fy,fz;
    
   @@ -239,8 +268,9 @@ void PPPMTIP4P::fieldforce_ik()
      // (mx,my,mz) = global coords of moving stencil pt
      // ek = 3 components of E-field on particle
      double *q = atom->q;
   -  double **x = atom->x;
   -  double **f = atom->f;
   +  const auto * _noalias const x = (dbl3_t *) atom->x[0];
   +  auto * _noalias const f = (dbl3_t *) atom->f[0];
   +  const auto * _noalias const p2g = (int3_t *) part2grid[0];
      const double half_alpha = 0.5 * alpha;
      const double one_minus_alpha = 1.0 - alpha;
      const int zforce_active = (slabflag != 2);
   @@ -253,19 +283,29 @@ void PPPMTIP4P::fieldforce_ik()
      const int nx_out = nxhi_out - nxlo_out + 1;
      const int ny_out = nyhi_out - nylo_out + 1;
      const int nyz_out = nx_out * ny_out;
   +  const double boxlox = boxlo[0];
   +  const double boxloy = boxlo[1];
   +  const double boxloz = boxlo[2];
    
      for (i = 0; i < nlocal; i++) {
   +    double xix, xiy, xiz;
        if (type[i] == typeO) {
          find_M_cached(i,iH1,iH2,xM);
   -      xi = xM;
   -    } else xi = x[i];
   +      xix = xM[0];
   +      xiy = xM[1];
   +      xiz = xM[2];
   +    } else {
   +      xix = x[i].x;
   +      xiy = x[i].y;
   +      xiz = x[i].z;
   +    }
    
   -    nx = part2grid[i][0];
   -    ny = part2grid[i][1];
   -    nz = part2grid[i][2];
   -    dx = nx+shiftone - (xi[0]-boxlo[0])*delxinv;
   -    dy = ny+shiftone - (xi[1]-boxlo[1])*delyinv;
   -    dz = nz+shiftone - (xi[2]-boxlo[2])*delzinv;
   +    nx = p2g[i].a;
   +    ny = p2g[i].b;
   +    nz = p2g[i].t;
   +    dx = nx+shiftone - (xix-boxlox)*delxinv;
   +    dy = ny+shiftone - (xiy-boxloy)*delyinv;
   +    dz = nz+shiftone - (xiz-boxloz)*delzinv;
    
        compute_rho1d(dx,dy,dz);
    
   @@ -290,26 +330,26 @@ void PPPMTIP4P::fieldforce_ik()
    
        const double qfactor = qqrd2e * scale * q[i];
        if (type[i] != typeO) {
   -      f[i][0] += qfactor*ekx;
   -      f[i][1] += qfactor*eky;
   -      if (zforce_active) f[i][2] += qfactor*ekz;
   +      f[i].x += qfactor*ekx;
   +      f[i].y += qfactor*eky;
   +      if (zforce_active) f[i].z += qfactor*ekz;
    
        } else {
          fx = qfactor * ekx;
          fy = qfactor * eky;
          fz = qfactor * ekz;
    
   -      f[i][0] += fx * one_minus_alpha;
   -      f[i][1] += fy * one_minus_alpha;
   -      if (zforce_active) f[i][2] += fz * one_minus_alpha;
   +      f[i].x += fx * one_minus_alpha;
   +      f[i].y += fy * one_minus_alpha;
   +      if (zforce_active) f[i].z += fz * one_minus_alpha;
    
   -      f[iH1][0] += half_alpha * fx;
   -      f[iH1][1] += half_alpha * fy;
   -      if (zforce_active) f[iH1][2] += half_alpha * fz;
   +      f[iH1].x += half_alpha * fx;
   +      f[iH1].y += half_alpha * fy;
   +      if (zforce_active) f[iH1].z += half_alpha * fz;
    
   -      f[iH2][0] += half_alpha * fx;
   -      f[iH2][1] += half_alpha * fy;
   -      if (zforce_active) f[iH2][2] += half_alpha * fz;
   +      f[iH2].x += half_alpha * fx;
   +      f[iH2].y += half_alpha * fy;
   +      if (zforce_active) f[iH2].z += half_alpha * fz;
        }
      }
    }


.. _iter-0021-page-head-2599e1f9f4c8: https://github.com/skilled-scipkg/lammps/commit/2599e1f9f4c8
.. _iter-0021-guardrails-incumbent-ab55ea3e8508820bbebff62c9b733c8990e9a801: https://github.com/skilled-scipkg/lammps/commit/ab55ea3e8508820bbebff62c9b733c8990e9a801
.. _iter-0021-guardrails-candidate-2599e1f9f4c8bb574886ac80f8e67052265a0ad4: https://github.com/skilled-scipkg/lammps/commit/2599e1f9f4c8bb574886ac80f8e67052265a0ad4