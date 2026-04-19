Iteration 0003 — 3a45586c8afa (accepted)
========================================


Change summary
--------------


Specialize `pair_lj_cut_tip4p_long` with compile-time energy/virial/table dispatch, a shared TIP4P site-cache helper, and separate oxygen/non-oxygen `i` hot paths while preserving the original force writes and `ev_tally_tip4p()` semantics.

Acceptance rationale
--------------------


Correctness passed with no guardrail errors, and the candidate improves the primary metric from 94.717 to 81.885, a 13.55% gain that clears the 2% acceptance bar.

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
| incumbent commit | e7c0ed95a333                                                               |
+------------------+----------------------------------------------------------------------------+
| candidate commit | 3a45586c8afa                                                               |
+------------------+----------------------------------------------------------------------------+
| incumbent metric | 94.717                                                                     |
+------------------+----------------------------------------------------------------------------+
| candidate metric | 81.885                                                                     |
+------------------+----------------------------------------------------------------------------+
| baseline metric  | 94.717                                                                     |
+------------------+----------------------------------------------------------------------------+
| Δ vs incumbent   | +13.548% (lower-is-better sign)                                            |
+------------------+----------------------------------------------------------------------------+
| changed files    | src/KSPACE/pair_lj_cut_tip4p_long.cpp, src/KSPACE/pair_lj_cut_tip4p_long.h |
+------------------+----------------------------------------------------------------------------+


Diffstat
--------


.. code-block:: text

    src/KSPACE/pair_lj_cut_tip4p_long.cpp | 616 +++++++++++++++++-----------------
    src/KSPACE/pair_lj_cut_tip4p_long.h   |   5 +
    2 files changed, 322 insertions(+), 299 deletions(-)

Diff
----


:download:`download full diff <_diffs/iter_0003_3a45586c8afa.diff>`

.. code-block:: diff

   diff --git a/src/KSPACE/pair_lj_cut_tip4p_long.cpp b/src/KSPACE/pair_lj_cut_tip4p_long.cpp
   index d9f8daf751..5534584ad0 100644
   --- a/src/KSPACE/pair_lj_cut_tip4p_long.cpp
   +++ b/src/KSPACE/pair_lj_cut_tip4p_long.cpp
   @@ -71,28 +71,14 @@ PairLJCutTIP4PLong::~PairLJCutTIP4PLong()
    
    void PairLJCutTIP4PLong::compute(int eflag, int vflag)
    {
   -  int i,j,ii,jj,inum,jnum,itype,jtype,itable,key;
   -  int n,vlist[6];
   -  int iH1,iH2,jH1,jH2;
   -  double qtmp,xtmp,ytmp,ztmp,delx,dely,delz,evdwl,ecoul;
   -  double fraction,table;
   -  double r,r2inv,r6inv,forcecoul,forcelj,cforce;
   -  double factor_coul,factor_lj;
   -  double grij,expm2,prefactor,t,erfc;
   -  double fO[3],fH[3],fd[3],v[6];
   -  double *x1,*x2,*xH1,*xH2;
   -  int *ilist,*jlist,*numneigh,**firstneigh;
   -  double rsq;
   -
   -  evdwl = ecoul = 0.0;
      ev_init(eflag,vflag);
    
      // reallocate hneigh & newsite if necessary
      // initialize hneigh[0] to -1 on steps when reneighboring occurred
      // initialize hneigh[2] to 0 every step
    
   -  int nlocal = atom->nlocal;
   -  int nall = nlocal + atom->nghost;
   +  const int nlocal = atom->nlocal;
   +  const int nall = nlocal + atom->nghost;
    
      if (atom->nmax > nmax) {
        nmax = atom->nmax;
   @@ -102,321 +88,353 @@ void PairLJCutTIP4PLong::compute(int eflag, int vflag)
        memory->create(newsite,nmax,3,"pair:newsite");
      }
      if (neighbor->ago == 0)
   -    for (i = 0; i < nall; i++) hneigh[i][0] = -1;
   -  for (i = 0; i < nall; i++) hneigh[i][2] = 0;
   +    for (int i = 0; i < nall; i++) hneigh[i][0] = -1;
   +  for (int i = 0; i < nall; i++) hneigh[i][2] = 0;
   +
   +  if (!ncoultablebits) {
   +    if (evflag) {
   +      if (eflag) {
   +        if (vflag) return eval<1,1,1,1>();
   +        else return eval<1,1,1,0>();
   +      } else {
   +        if (vflag) return eval<1,1,0,1>();
   +        else return eval<1,1,0,0>();
   +      }
   +    } else return eval<1,0,0,0>();
   +  } else {
   +    if (evflag) {
   +      if (eflag) {
   +        if (vflag) return eval<0,1,1,1>();
   +        else return eval<0,1,1,0>();
   +      } else {
   +        if (vflag) return eval<0,1,0,1>();
   +        else return eval<0,1,0,0>();
   +      }
   +    } else return eval<0,0,0,0>();
   +  }
   +}
    
   -  double **f = atom->f;
   +/* ---------------------------------------------------------------------- */
   +
   +template <int CTABLE, int EVFLAG, int EFLAG, int VFLAG>
   +void PairLJCutTIP4PLong::eval()
   +{
      double **x = atom->x;
      double *q = atom->q;
      tagint *tag = atom->tag;
      int *type = atom->type;
   -  double *special_coul = force->special_coul;
   -  double *special_lj = force->special_lj;
   -  int newton_pair = force->newton_pair;
   -  double qqrd2e = force->qqrd2e;
   -  double cut_coulsqplus = (cut_coul+2.0*qdist) * (cut_coul+2.0*qdist);
   +  int *ilist = list->ilist;
   +  int *numneigh = list->numneigh;
   +  int **firstneigh = list->firstneigh;
   +  const int inum = list->inum;
   +  const int nlocal = atom->nlocal;
   +  const double cut_coulsqplus = (cut_coul+2.0*qdist) * (cut_coul+2.0*qdist);
   +
   +  for (int ii = 0; ii < inum; ii++) {
   +    const int i = ilist[ii];
   +    const double qtmp = q[i];
   +    const double xtmp = x[i][0];
   +    const double ytmp = x[i][1];
   +    const double ztmp = x[i][2];
   +    const int itype = type[i];
   +    int *jlist = firstneigh[i];
   +    const int jnum = numneigh[i];
    
   -  inum = list->inum;
   -  ilist = list->ilist;
   -  numneigh = list->numneigh;
   -  firstneigh = list->firstneigh;
   -
   -  // loop over neighbors of my atoms
   +    if (itype == typeO) {
   +      int iH1, iH2;
   +      double *x1 = tip4p_site(i,iH1,iH2,x,tag,type);
   +      eval_i<CTABLE,EVFLAG,EFLAG,VFLAG,1>(i,itype,qtmp,xtmp,ytmp,ztmp,iH1,iH2,x1,jlist,jnum,
   +                                          nlocal,cut_coulsqplus);
   +    } else {
   +      eval_i<CTABLE,EVFLAG,EFLAG,VFLAG,0>(i,itype,qtmp,xtmp,ytmp,ztmp,-1,-1,x[i],jlist,jnum,
   +                                          nlocal,cut_coulsqplus);
   +    }
   +  }
   +}
    
   -  for (ii = 0; ii < inum; ii++) {
   -    i = ilist[ii];
   -    qtmp = q[i];
   -    xtmp = x[i][0];
   -    ytmp = x[i][1];
   -    ztmp = x[i][2];
   -    itype = type[i];
   +/* ---------------------------------------------------------------------- */
    
   -    // if atom I = water O, set x1 = offset charge site
   -    // else x1 = x of atom I
   +template <int CTABLE, int EVFLAG, int EFLAG, int VFLAG, int I_WATER>
   +void PairLJCutTIP4PLong::eval_i(int i, int itype, double qtmp, double xtmp, double ytmp,
   +                                double ztmp, int iH1, int iH2, double *x1, int *jlist,
   +                                int jnum, int nlocal, double cut_coulsqplus)
   +{
   +  int j, jj, jtype, itable, key = 0, n = 0, jH1 = -1, jH2 = -1;
   +  int vlist[6];
   +  double delx, dely, delz, evdwl, ecoul, fraction, table;
   +  double r, rsq, r2inv, r6inv, forcecoul, forcelj, cforce;
   +  double factor_coul, factor_lj;
   +  double grij, expm2, prefactor, t, erfc;
   +  double fO[3], fH[3], fd[3], v[6];
   +  double *x2, *xH1, *xH2;
    
   -    if (itype == typeO) {
   -      if (hneigh[i][0] < 0) {
   -        iH1 = atom->map(tag[i] + 1);
   -        iH2 = atom->map(tag[i] + 2);
   -        if (iH1 == -1 || iH2 == -1)
   -          error->one(FLERR,"TIP4P hydrogen is missing");
   -        if (atom->type[iH1] != typeH || atom->type[iH2] != typeH)
   -          error->one(FLERR,"TIP4P hydrogen has incorrect atom type");
   -        // set iH1,iH2 to closest image to O
   -        iH1 = domain->closest_image(i,iH1);
   -        iH2 = domain->closest_image(i,iH2);
   -        compute_newsite(x[i],x[iH1],x[iH2],newsite[i]);
   -        hneigh[i][0] = iH1;
   -        hneigh[i][1] = iH2;
   -        hneigh[i][2] = 1;
   +  double **f = atom->f;
   +  double **x = atom->x;
   +  double *q = atom->q;
   +  tagint *tag = atom->tag;
   +  int *type = atom->type;
   +  double *special_coul = force->special_coul;
   +  double *special_lj = force->special_lj;
   +  const double qqrd2e = force->qqrd2e;
   +  const double half_alpha = 0.5 * alpha;
   +  const double one_minus_alpha = 1.0 - alpha;
   +  double *cut_ljsq_i = cut_ljsq[itype];
   +  double *lj1_i = lj1[itype];
   +  double *lj2_i = lj2[itype];
   +  double *lj3_i = lj3[itype];
   +  double *lj4_i = lj4[itype];
   +  double *offset_i = offset[itype];
   +
   +  for (jj = 0; jj < jnum; jj++) {
   +    j = jlist[jj];
   +    factor_lj = special_lj[sbmask(j)];
   +    factor_coul = special_coul[sbmask(j)];
   +    j &= NEIGHMASK;
   +
   +    delx = xtmp - x[j][0];
   +    dely = ytmp - x[j][1];
   +    delz = ztmp - x[j][2];
   +    rsq = delx*delx + dely*dely + delz*delz;
   +    jtype = type[j];
   +
   +    if (rsq < cut_ljsq_i[jtype]) {
   +      r2inv = 1.0/rsq;
   +      r6inv = r2inv*r2inv*r2inv;
   +      forcelj = r6inv * (lj1_i[jtype]*r6inv - lj2_i[jtype]);
   +      forcelj *= factor_lj * r2inv;
   +
   +      f[i][0] += delx*forcelj;
   +      f[i][1] += dely*forcelj;
   +      f[i][2] += delz*forcelj;
   +      f[j][0] -= delx*forcelj;
   +      f[j][1] -= dely*forcelj;
   +      f[j][2] -= delz*forcelj;
   +
   +      if constexpr (EFLAG) {
   +        evdwl = r6inv*(lj3_i[jtype]*r6inv-lj4_i[jtype]) - offset_i[jtype];
   +        evdwl *= factor_lj;
   +      } else evdwl = 0.0;
   +
   +      if constexpr (EVFLAG)
   +        ev_tally(i,j,nlocal,1,evdwl,0.0,forcelj,delx,dely,delz);
   +    }
    
   -      } else {
   -        iH1 = hneigh[i][0];
   -        iH2 = hneigh[i][1];
   -        if (hneigh[i][2] == 0) {
   -          hneigh[i][2] = 1;
   -          compute_newsite(x[i],x[iH1],x[iH2],newsite[i]);
   -        }
   -      }
   -      x1 = newsite[i];
   -    } else x1 = x[i];
   -
   -    jlist = firstneigh[i];
   -    jnum = numneigh[i];
   -
   -    for (jj = 0; jj < jnum; jj++) {
   -      j = jlist[jj];
   -      factor_lj = special_lj[sbmask(j)];
   -      factor_coul = special_coul[sbmask(j)];
   -      j &= NEIGHMASK;
   -
   -      delx = xtmp - x[j][0];
   -      dely = ytmp - x[j][1];
   -      delz = ztmp - x[j][2];
   -      rsq = delx*delx + dely*dely + delz*delz;
   -      jtype = type[j];
   -
   -      // LJ interaction based on true rsq
   -
   -      if (rsq < cut_ljsq[itype][jtype]) {
   -        r2inv = 1.0/rsq;
   -        r6inv = r2inv*r2inv*r2inv;
   -        forcelj = r6inv * (lj1[itype][jtype]*r6inv - lj2[itype][jtype]);
   -        forcelj *= factor_lj * r2inv;
   -
   -        f[i][0] += delx*forcelj;
   -        f[i][1] += dely*forcelj;
   -        f[i][2] += delz*forcelj;
   -        f[j][0] -= delx*forcelj;
   -        f[j][1] -= dely*forcelj;
   -        f[j][2] -= delz*forcelj;
   -
   -        if (eflag) {
   -          evdwl = r6inv*(lj3[itype][jtype]*r6inv-lj4[itype][jtype]) -
   -            offset[itype][jtype];
   -          evdwl *= factor_lj;
   -        } else evdwl = 0.0;
   -
   -        if (evflag) ev_tally(i,j,nlocal,newton_pair,
   -                             evdwl,0.0,forcelj,delx,dely,delz);
   +    if (rsq < cut_coulsqplus) {
   +      if constexpr (I_WATER) {
   +        if (jtype == typeO) x2 = tip4p_site(j,jH1,jH2,x,tag,type);
   +        else x2 = x[j];
   +
   +        delx = x1[0] - x2[0];
   +        dely = x1[1] - x2[1];
   +        delz = x1[2] - x2[2];
   +        rsq = delx*delx + dely*dely + delz*delz;
   +      } else if (jtype == typeO) {
   +        x2 = tip4p_site(j,jH1,jH2,x,tag,type);
   +        delx = x1[0] - x2[0];
   +        dely = x1[1] - x2[1];
   +        delz = x1[2] - x2[2];
   +        rsq = delx*delx + dely*dely + delz*delz;
          }
    
   -      // adjust rsq and delxyz for off-site O charge(s) if necessary
   -      // but only if they are within reach
   -
   -      if (rsq < cut_coulsqplus) {
   -        if (itype == typeO || jtype == typeO) {
   -
   -          // if atom J = water O, set x2 = offset charge site
   -          // else x2 = x of atom J
   -
   -          if (jtype == typeO) {
   -            if (hneigh[j][0] < 0) {
   -              jH1 = atom->map(tag[j] + 1);
   -              jH2 = atom->map(tag[j] + 2);
   -              if (jH1 == -1 || jH2 == -1)
   -                error->one(FLERR,"TIP4P hydrogen is missing");
   -              if (atom->type[jH1] != typeH || atom->type[jH2] != typeH)
   -                error->one(FLERR,"TIP4P hydrogen has incorrect atom type");
   -              // set jH1,jH2 to closest image to O
   -              jH1 = domain->closest_image(j,jH1);
   -              jH2 = domain->closest_image(j,jH2);
   -              compute_newsite(x[j],x[jH1],x[jH2],newsite[j]);
   -              hneigh[j][0] = jH1;
   -              hneigh[j][1] = jH2;
   -              hneigh[j][2] = 1;
   -
   -            } else {
   -              jH1 = hneigh[j][0];
   -              jH2 = hneigh[j][1];
   -              if (hneigh[j][2] == 0) {
   -                hneigh[j][2] = 1;
   -                compute_newsite(x[j],x[jH1],x[jH2],newsite[j]);
   -              }
   -            }
   -            x2 = newsite[j];
   -          } else x2 = x[j];
   -
   -          delx = x1[0] - x2[0];
   -          dely = x1[1] - x2[1];
   -          delz = x1[2] - x2[2];
   -          rsq = delx*delx + dely*dely + delz*delz;
   -        }
   -
   -        // Coulombic interaction based on modified rsq
   -
   -        if (rsq < cut_coulsq) {
   -          r2inv = 1 / rsq;
   -          if (!ncoultablebits || rsq <= tabinnersq) {
   -            r = sqrt(rsq);
   -            grij = g_ewald * r;
   -            expm2 = exp(-grij*grij);
   -            t = 1.0 / (1.0 + EWALD_P*grij);
   -            erfc = t * (A1+t*(A2+t*(A3+t*(A4+t*A5)))) * expm2;
   -            prefactor = qqrd2e * qtmp*q[j]/r;
   -            forcecoul = prefactor * (erfc + EWALD_F*grij*expm2);
   -            if (factor_coul < 1.0) {
   -              forcecoul -= (1.0-factor_coul)*prefactor;
   -            }
   -          } else {
   -            union_int_float_t rsq_lookup;
   -            rsq_lookup.f = rsq;
   -            itable = rsq_lookup.i & ncoulmask;
   -            itable >>= ncoulshiftbits;
   -            fraction = ((double) rsq_lookup.f - rtable[itable]) * drtable[itable];
   -            table = ftable[itable] + fraction*dftable[itable];
   -            forcecoul = qtmp*q[j] * table;
   -            if (factor_coul < 1.0) {
   -              table = ctable[itable] + fraction*dctable[itable];
   -              prefactor = qtmp*q[j] * table;
   -              forcecoul -= (1.0-factor_coul)*prefactor;
   -            }
   +      if (rsq < cut_coulsq) {
   +        r2inv = 1.0 / rsq;
   +        if (CTABLE || rsq <= tabinnersq) {
   +          r = sqrt(rsq);
   +          grij = g_ewald * r;
   +          expm2 = exp(-grij*grij);
   +          t = 1.0 / (1.0 + EWALD_P*grij);
   +          erfc = t * (A1+t*(A2+t*(A3+t*(A4+t*A5)))) * expm2;
   +          prefactor = qqrd2e * qtmp*q[j]/r;
   +          forcecoul = prefactor * (erfc + EWALD_F*grij*expm2);
   +          if (factor_coul < 1.0) {
   +            forcecoul -= (1.0-factor_coul)*prefactor;
              }
   +        } else {
   +          union_int_float_t rsq_lookup;
   +          rsq_lookup.f = rsq;
   +          itable = rsq_lookup.i & ncoulmask;
   +          itable >>= ncoulshiftbits;
   +          fraction = ((double) rsq_lookup.f - rtable[itable]) * drtable[itable];
   +          table = ftable[itable] + fraction*dftable[itable];
   +          forcecoul = qtmp*q[j] * table;
   +          if (factor_coul < 1.0) {
   +            table = ctable[itable] + fraction*dctable[itable];
   +            prefactor = qtmp*q[j] * table;
   +            forcecoul -= (1.0-factor_coul)*prefactor;
   +          }
   +        }
    
   -          cforce = forcecoul * r2inv;
   -
   -          // if i,j are not O atoms, force is applied directly
   -          // if i or j are O atoms, force is on fictitious atom & partitioned
   -          // force partitioning due to Feenstra, J Comp Chem, 20, 786 (1999)
   -          // f_f = fictitious force, fO = f_f (1 - 2 alpha), fH = alpha f_f
   -          // preserves total force and torque on water molecule
   -          // virial = sum(r x F) where each water's atoms are near xi and xj
   -          // vlist stores 2,4,6 atoms whose forces contribute to virial
   +        cforce = forcecoul * r2inv;
    
   +        if constexpr (EVFLAG) {
              n = 0;
   -          key = 0;
   -
   -          if (itype != typeO) {
   -            f[i][0] += delx * cforce;
   -            f[i][1] += dely * cforce;
   -            f[i][2] += delz * cforce;
   -
   -            if (vflag) {
   -              v[0] = x[i][0] * delx * cforce;
   -              v[1] = x[i][1] * dely * cforce;
   -              v[2] = x[i][2] * delz * cforce;
   -              v[3] = x[i][0] * dely * cforce;
   -              v[4] = x[i][0] * delz * cforce;
   -              v[5] = x[i][1] * delz * cforce;
   -            }
   -            vlist[n++] = i;
   +          key = I_WATER;
   +        }
    
   -          } else {
   -            key++;
   -
   -            fd[0] = delx*cforce;
   -            fd[1] = dely*cforce;
   -            fd[2] = delz*cforce;
   -
   -            fO[0] = fd[0]*(1 - alpha);
   -            fO[1] = fd[1]*(1 - alpha);
   -            fO[2] = fd[2]*(1 - alpha);
   -
   -            fH[0] = 0.5 * alpha * fd[0];
   -            fH[1] = 0.5 * alpha * fd[1];
   -            fH[2] = 0.5 * alpha * fd[2];
   -
   -            f[i][0] += fO[0];
   -            f[i][1] += fO[1];
   -            f[i][2] += fO[2];
   -
   -            f[iH1][0] += fH[0];
   -            f[iH1][1] += fH[1];
   -            f[iH1][2] += fH[2];
   -
   -            f[iH2][0] += fH[0];
   -            f[iH2][1] += fH[1];
   -            f[iH2][2] += fH[2];
   -
   -            if (vflag) {
   -              xH1 = x[iH1];
   -              xH2 = x[iH2];
   -              v[0] = x[i][0]*fO[0] + xH1[0]*fH[0] + xH2[0]*fH[0];
   -              v[1] = x[i][1]*fO[1] + xH1[1]*fH[1] + xH2[1]*fH[1];
   -              v[2] = x[i][2]*fO[2] + xH1[2]*fH[2] + xH2[2]*fH[2];
   -              v[3] = x[i][0]*fO[1] + xH1[0]*fH[1] + xH2[0]*fH[1];
   -              v[4] = x[i][0]*fO[2] + xH1[0]*fH[2] + xH2[0]*fH[2];
   -              v[5] = x[i][1]*fO[2] + xH1[1]*fH[2] + xH2[1]*fH[2];
   -            }
   +        if constexpr (!I_WATER) {
   +          f[i][0] += delx * cforce;
   +          f[i][1] += dely * cforce;
   +          f[i][2] += delz * cforce;
   +
   +          if constexpr (VFLAG) {
   +            v[0] = xtmp * delx * cforce;
   +            v[1] = ytmp * dely * cforce;
   +            v[2] = ztmp * delz * cforce;
   +            v[3] = xtmp * dely * cforce;
   +            v[4] = xtmp * delz * cforce;
   +            v[5] = ytmp * delz * cforce;
   +          }
   +          if constexpr (EVFLAG) vlist[n++] = i;
   +
   +        } else {
   +          fd[0] = delx*cforce;
   +          fd[1] = dely*cforce;
   +          fd[2] = delz*cforce;
   +
   +          fO[0] = fd[0] * one_minus_alpha;
   +          fO[1] = fd[1] * one_minus_alpha;
   +          fO[2] = fd[2] * one_minus_alpha;
   +
   +          fH[0] = half_alpha * fd[0];
   +          fH[1] = half_alpha * fd[1];
   +          fH[2] = half_alpha * fd[2];
   +
   +          f[i][0] += fO[0];
   +          f[i][1] += fO[1];
   +          f[i][2] += fO[2];
   +
   +          f[iH1][0] += fH[0];
   +          f[iH1][1] += fH[1];
   +          f[iH1][2] += fH[2];
   +
   +          f[iH2][0] += fH[0];
   +          f[iH2][1] += fH[1];
   +          f[iH2][2] += fH[2];
   +
   +          if constexpr (VFLAG) {
   +            xH1 = x[iH1];
   +            xH2 = x[iH2];
   +            v[0] = xtmp*fO[0] + xH1[0]*fH[0] + xH2[0]*fH[0];
   +            v[1] = ytmp*fO[1] + xH1[1]*fH[1] + xH2[1]*fH[1];
   +            v[2] = ztmp*fO[2] + xH1[2]*fH[2] + xH2[2]*fH[2];
   +            v[3] = xtmp*fO[1] + xH1[0]*fH[1] + xH2[0]*fH[1];
   +            v[4] = xtmp*fO[2] + xH1[0]*fH[2] + xH2[0]*fH[2];
   +            v[5] = ytmp*fO[2] + xH1[1]*fH[2] + xH2[1]*fH[2];
   +          }
   +          if constexpr (EVFLAG) {
                vlist[n++] = i;
                vlist[n++] = iH1;
                vlist[n++] = iH2;
              }
   +        }
    
   -          if (jtype != typeO) {
   -            f[j][0] -= delx * cforce;
   -            f[j][1] -= dely * cforce;
   -            f[j][2] -= delz * cforce;
   -
   -            if (vflag) {
   -              v[0] -= x[j][0] * delx * cforce;
   -              v[1] -= x[j][1] * dely * cforce;
   -              v[2] -= x[j][2] * delz * cforce;
   -              v[3] -= x[j][0] * dely * cforce;
   -              v[4] -= x[j][0] * delz * cforce;
   -              v[5] -= x[j][1] * delz * cforce;
   -            }
   -            vlist[n++] = j;
   -
   -          } else {
   -            key += 2;
   -
   -            fd[0] = -delx*cforce;
   -            fd[1] = -dely*cforce;
   -            fd[2] = -delz*cforce;
   -
   -            fO[0] = fd[0]*(1 - alpha);
   -            fO[1] = fd[1]*(1 - alpha);
   -            fO[2] = fd[2]*(1 - alpha);
   -
   -            fH[0] = 0.5 * alpha * fd[0];
   -            fH[1] = 0.5 * alpha * fd[1];
   -            fH[2] = 0.5 * alpha * fd[2];
   -
   -            f[j][0] += fO[0];
   -            f[j][1] += fO[1];
   -            f[j][2] += fO[2];
   -
   -            f[jH1][0] += fH[0];
   -            f[jH1][1] += fH[1];
   -            f[jH1][2] += fH[2];
   -
   -            f[jH2][0] += fH[0];
   -            f[jH2][1] += fH[1];
   -            f[jH2][2] += fH[2];
   -
   -            if (vflag) {
   -              xH1 = x[jH1];
   -              xH2 = x[jH2];
   -              v[0] += x[j][0]*fO[0] + xH1[0]*fH[0] + xH2[0]*fH[0];
   -              v[1] += x[j][1]*fO[1] + xH1[1]*fH[1] + xH2[1]*fH[1];
   -              v[2] += x[j][2]*fO[2] + xH1[2]*fH[2] + xH2[2]*fH[2];
   -              v[3] += x[j][0]*fO[1] + xH1[0]*fH[1] + xH2[0]*fH[1];
   -              v[4] += x[j][0]*fO[2] + xH1[0]*fH[2] + xH2[0]*fH[2];
   -              v[5] += x[j][1]*fO[2] + xH1[1]*fH[2] + xH2[1]*fH[2];
   -            }
   +        if (jtype != typeO) {
   +          f[j][0] -= delx * cforce;
   +          f[j][1] -= dely * cforce;
   +          f[j][2] -= delz * cforce;
   +
   +          if constexpr (VFLAG) {
   +            v[0] -= x[j][0] * delx * cforce;
   +            v[1] -= x[j][1] * dely * cforce;
   +            v[2] -= x[j][2] * delz * cforce;
   +            v[3] -= x[j][0] * dely * cforce;
   +            v[4] -= x[j][0] * delz * cforce;
   +            v[5] -= x[j][1] * delz * cforce;
   +          }
   +          if constexpr (EVFLAG) vlist[n++] = j;
   +
   +        } else {
   +          if constexpr (EVFLAG) key += 2;
   +
   +          fd[0] = -delx*cforce;
   +          fd[1] = -dely*cforce;
   +          fd[2] = -delz*cforce;
   +
   +          fO[0] = fd[0] * one_minus_alpha;
   +          fO[1] = fd[1] * one_minus_alpha;
   +          fO[2] = fd[2] * one_minus_alpha;
   +
   +          fH[0] = half_alpha * fd[0];
   +          fH[1] = half_alpha * fd[1];
   +          fH[2] = half_alpha * fd[2];
   +
   +          f[j][0] += fO[0];
   +          f[j][1] += fO[1];
   +          f[j][2] += fO[2];
   +
   +          f[jH1][0] += fH[0];
   +          f[jH1][1] += fH[1];
   +          f[jH1][2] += fH[2];
   +
   +          f[jH2][0] += fH[0];
   +          f[jH2][1] += fH[1];
   +          f[jH2][2] += fH[2];
   +
   +          if constexpr (VFLAG) {
   +            xH1 = x[jH1];
   +            xH2 = x[jH2];
   +            v[0] += x[j][0]*fO[0] + xH1[0]*fH[0] + xH2[0]*fH[0];
   +            v[1] += x[j][1]*fO[1] + xH1[1]*fH[1] + xH2[1]*fH[1];
   +            v[2] += x[j][2]*fO[2] + xH1[2]*fH[2] + xH2[2]*fH[2];
   +            v[3] += x[j][0]*fO[1] + xH1[0]*fH[1] + xH2[0]*fH[1];
   +            v[4] += x[j][0]*fO[2] + xH1[0]*fH[2] + xH2[0]*fH[2];
   ... (88 more lines; see download link above)
