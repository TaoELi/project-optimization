Iteration 0016 — 9f1f066ed3fb (accepted)
========================================


Change summary
--------------


Refactor `pair_lj_cut_tip4p_long::eval()` Coulomb off-site reach handling to early-continue `i_is_O`/`j_is_O` gating and cache `q[j]` once per Coulomb interaction, while preserving existing cutoff and TIP4P force/virial semantics.

Acceptance rationale
--------------------


Primary metric improved by +20.50% vs incumbent (0.264975531201 -> 0.210648385517) with correctness/guardrails passing and no hard reject.

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
| incumbent commit | 01ed31155d08                          |
+------------------+---------------------------------------+
| candidate commit | 9f1f066ed3fb                          |
+------------------+---------------------------------------+
| incumbent metric | 0.264976                              |
+------------------+---------------------------------------+
| candidate metric | 0.210648                              |
+------------------+---------------------------------------+
| baseline metric  | 0.34835                               |
+------------------+---------------------------------------+
| Δ vs incumbent   | +20.503% (lower-is-better sign)       |
+------------------+---------------------------------------+
| changed files    | src/KSPACE/pair_lj_cut_tip4p_long.cpp |
+------------------+---------------------------------------+


Diffstat
--------


.. code-block:: text

    src/KSPACE/pair_lj_cut_tip4p_long.cpp | 382 +++++++++++++++++-----------------
    1 file changed, 192 insertions(+), 190 deletions(-)

Diff
----


:download:`download full diff <_diffs/iter_0016_9f1f066ed3fb.diff>`

.. code-block:: diff

   diff --git a/src/KSPACE/pair_lj_cut_tip4p_long.cpp b/src/KSPACE/pair_lj_cut_tip4p_long.cpp
   index a436d3aafd..d191123e11 100644
   --- a/src/KSPACE/pair_lj_cut_tip4p_long.cpp
   +++ b/src/KSPACE/pair_lj_cut_tip4p_long.cpp
   @@ -276,211 +276,213 @@ void PairLJCutTIP4PLong::eval()
          // adjust rsq and delxyz for off-site O charge(s) if necessary
          // but only if they are within reach
    
   -      const double cut_coulsq_gate = (i_is_O || j_is_O) ? cut_coulsqplus :
   -        cut_coulsq;
   -
   -      if (rsq < cut_coulsq_gate) {
   -        if (i_is_O) {
   -
   -          // if atom J = water O, set x2 = offset charge site
   -          // else x2 = x of atom J
   -
   -          if (j_is_O) {
   -            jH1 = hneigh2d[j][0];
   -            jH2 = hneigh2d[j][1];
   -            x2 = newsite2d[j];
   -          } else x2 = x[j];
   -
   -          delx = x1[0] - x2[0];
   -          dely = x1[1] - x2[1];
   -          delz = x1[2] - x2[2];
   -          rsq = delx*delx + dely*dely + delz*delz;
   -        } else if (j_is_O) {
   +      if (i_is_O) {
   +        if (rsq >= cut_coulsqplus) continue;
   +
   +        // if atom J = water O, set x2 = offset charge site
   +        // else x2 = x of atom J
   +
   +        if (j_is_O) {
              jH1 = hneigh2d[j][0];
              jH2 = hneigh2d[j][1];
              x2 = newsite2d[j];
   +        } else x2 = xj;
   +
   +        delx = x1[0] - x2[0];
   +        dely = x1[1] - x2[1];
   +        delz = x1[2] - x2[2];
   +        rsq = delx*delx + dely*dely + delz*delz;
   +      } else if (j_is_O) {
   +        if (rsq >= cut_coulsqplus) continue;
   +
   +        jH1 = hneigh2d[j][0];
   +        jH2 = hneigh2d[j][1];
   +        x2 = newsite2d[j];
   +
   +        delx = x1[0] - x2[0];
   +        dely = x1[1] - x2[1];
   +        delz = x1[2] - x2[2];
   +        rsq = delx*delx + dely*dely + delz*delz;
   +      } else if (rsq >= cut_coulsq) {
   +        continue;
   +      }
    
   -          delx = x1[0] - x2[0];
   -          dely = x1[1] - x2[1];
   -          delz = x1[2] - x2[2];
   -          rsq = delx*delx + dely*dely + delz*delz;
   +      // Coulombic interaction based on modified rsq
   +
   +      if (rsq < cut_coulsq) {
   +        const double qj = q[j];
   +        double qiqj = 0.0;
   +        factor_coul = 1.0;
   +        r2inv = 1 / rsq;
   +        if (CTABLE || rsq <= tabinnersq) {
   +          r = sqrt(rsq);
   +          grij = g_ewald * r;
   +          expm2 = exp(-grij*grij);
   +          t = 1.0 / (1.0 + EWALD_P*grij);
   +          erfc = t * (A1+t*(A2+t*(A3+t*(A4+t*A5)))) * expm2;
   +          prefactor = qqrd2e_qi * qj / r;
   +          forcecoul = prefactor * (erfc + EWALD_F*grij*expm2);
   +          if (sb) factor_coul = special_coul[sb];
   +          if (factor_coul < 1.0) {
   +            forcecoul -= (1.0-factor_coul)*prefactor;
   +          }
   +        } else {
   +          union_int_float_t rsq_lookup;
   +          rsq_lookup.f = rsq;
   +          itable = rsq_lookup.i & ncoulmask;
   +          itable >>= ncoulshiftbits;
   +          fraction = ((double) rsq_lookup.f - rtable[itable]) * drtable[itable];
   +          table = ftable[itable] + fraction*dftable[itable];
   +          qiqj = qtmp * qj;
   +          forcecoul = qiqj * table;
   +          if (sb) factor_coul = special_coul[sb];
   +          if (factor_coul < 1.0) {
   +            table = ctable[itable] + fraction*dctable[itable];
   +            prefactor = qiqj * table;
   +            forcecoul -= (1.0-factor_coul)*prefactor;
   +          }
            }
    
   -        // Coulombic interaction based on modified rsq
   -
   -        if (rsq < cut_coulsq) {
   -          const double qiqj = qtmp * q[j];
   -          factor_coul = 1.0;
   -          r2inv = 1 / rsq;
   -          if (CTABLE || rsq <= tabinnersq) {
   -            r = sqrt(rsq);
   -            grij = g_ewald * r;
   -            expm2 = exp(-grij*grij);
   -            t = 1.0 / (1.0 + EWALD_P*grij);
   -            erfc = t * (A1+t*(A2+t*(A3+t*(A4+t*A5)))) * expm2;
   -            prefactor = qqrd2e_qi * q[j] / r;
   -            forcecoul = prefactor * (erfc + EWALD_F*grij*expm2);
   -            if (sb) factor_coul = special_coul[sb];
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
   -            forcecoul = qiqj * table;
   -            if (sb) factor_coul = special_coul[sb];
   -            if (factor_coul < 1.0) {
   -              table = ctable[itable] + fraction*dctable[itable];
   -              prefactor = qiqj * table;
   -              forcecoul -= (1.0-factor_coul)*prefactor;
   -            }
   -          }
   +        cforce = forcecoul * r2inv;
    
   -          cforce = forcecoul * r2inv;
   +        // if i,j are not O atoms, force is applied directly
   +        // if i or j are O atoms, force is on fictitious atom & partitioned
   +        // force partitioning due to Feenstra, J Comp Chem, 20, 786 (1999)
   +        // f_f = fictitious force, fO = f_f (1 - 2 alpha), fH = alpha f_f
   +        // preserves total force and torque on water molecule
   +        // virial = sum(r x F) where each water's atoms are near xi and xj
   +        // vlist stores 2,4,6 atoms whose forces contribute to virial
    
   -          // if i,j are not O atoms, force is applied directly
   -          // if i or j are O atoms, force is on fictitious atom & partitioned
   -          // force partitioning due to Feenstra, J Comp Chem, 20, 786 (1999)
   -          // f_f = fictitious force, fO = f_f (1 - 2 alpha), fH = alpha f_f
   -          // preserves total force and torque on water molecule
   -          // virial = sum(r x F) where each water's atoms are near xi and xj
   -          // vlist stores 2,4,6 atoms whose forces contribute to virial
   +        if (EVFLAG) {
   +          n = 0;
   +          key = 0;
   +        }
    
   +        if (!i_is_O) {
   +          fxtmp += delx * cforce;
   +          fytmp += dely * cforce;
   +          fztmp += delz * cforce;
   +
   +          if (VFLAG) {
   +            v[0] = xtmp * delx * cforce;
   +            v[1] = ytmp * dely * cforce;
   +            v[2] = ztmp * delz * cforce;
   +            v[3] = xtmp * dely * cforce;
   +            v[4] = xtmp * delz * cforce;
   +            v[5] = ytmp * delz * cforce;
   +          }
   +          if (EVFLAG) vlist[n++] = i;
   +
   +        } else {
   +          if (EVFLAG) key += 1;
   +
   +          fdx = delx*cforce;
   +          fdy = dely*cforce;
   +          fdz = delz*cforce;
   +
   +          fOx = fdx*alpha_O;
   +          fOy = fdy*alpha_O;
   +          fOz = fdz*alpha_O;
   +
   +          fHx = alpha_H * fdx;
   +          fHy = alpha_H * fdy;
   +          fHz = alpha_H * fdz;
   +
   +          fxtmp += fOx;
   +          fytmp += fOy;
   +          fztmp += fOz;
   +
   +          fiHxtmp += fHx;
   +          fiHytmp += fHy;
   +          fiHztmp += fHz;
   +
   +          if (VFLAG) {
   +            v[0] = xtmp*fOx + xiHsumx*fHx;
   +            v[1] = ytmp*fOy + xiHsumy*fHy;
   +            v[2] = ztmp*fOz + xiHsumz*fHz;
   +            v[3] = xtmp*fOy + xiHsumx*fHy;
   +            v[4] = xtmp*fOz + xiHsumx*fHz;
   +            v[5] = ytmp*fOz + xiHsumy*fHz;
   +          }
              if (EVFLAG) {
   -            n = 0;
   -            key = 0;
   +            vlist[n++] = i;
   +            vlist[n++] = iH1;
   +            vlist[n++] = iH2;
              }
   +        }
    
   -          if (!i_is_O) {
   -            fxtmp += delx * cforce;
   -            fytmp += dely * cforce;
   -            fztmp += delz * cforce;
   -
   -            if (VFLAG) {
   -              v[0] = xtmp * delx * cforce;
   -              v[1] = ytmp * dely * cforce;
   -              v[2] = ztmp * delz * cforce;
   -              v[3] = xtmp * dely * cforce;
   -              v[4] = xtmp * delz * cforce;
   -              v[5] = ytmp * delz * cforce;
   -            }
   -            if (EVFLAG) vlist[n++] = i;
   -
   -          } else {
   -            if (EVFLAG) key += 1;
   -
   -            fdx = delx*cforce;
   -            fdy = dely*cforce;
   -            fdz = delz*cforce;
   -
   -            fOx = fdx*alpha_O;
   -            fOy = fdy*alpha_O;
   -            fOz = fdz*alpha_O;
   -
   -            fHx = alpha_H * fdx;
   -            fHy = alpha_H * fdy;
   -            fHz = alpha_H * fdz;
   -
   -            fxtmp += fOx;
   -            fytmp += fOy;
   -            fztmp += fOz;
   -
   -            fiHxtmp += fHx;
   -            fiHytmp += fHy;
   -            fiHztmp += fHz;
   -
   -            if (VFLAG) {
   -              v[0] = xtmp*fOx + xiHsumx*fHx;
   -              v[1] = ytmp*fOy + xiHsumy*fHy;
   -              v[2] = ztmp*fOz + xiHsumz*fHz;
   -              v[3] = xtmp*fOy + xiHsumx*fHy;
   -              v[4] = xtmp*fOz + xiHsumx*fHz;
   -              v[5] = ytmp*fOz + xiHsumy*fHz;
   -            }
   -            if (EVFLAG) {
   -              vlist[n++] = i;
   -              vlist[n++] = iH1;
   -              vlist[n++] = iH2;
   -            }
   +        if (!j_is_O) {
   +          f[j][0] -= delx * cforce;
   +          f[j][1] -= dely * cforce;
   +          f[j][2] -= delz * cforce;
   +
   +          if (VFLAG) {
   +            v[0] -= xj[0] * delx * cforce;
   +            v[1] -= xj[1] * dely * cforce;
   +            v[2] -= xj[2] * delz * cforce;
   +            v[3] -= xj[0] * dely * cforce;
   +            v[4] -= xj[0] * delz * cforce;
   +            v[5] -= xj[1] * delz * cforce;
              }
   -
   -          if (!j_is_O) {
   -            f[j][0] -= delx * cforce;
   -            f[j][1] -= dely * cforce;
   -            f[j][2] -= delz * cforce;
   -
   -            if (VFLAG) {
   -              v[0] -= xj[0] * delx * cforce;
   -              v[1] -= xj[1] * dely * cforce;
   -              v[2] -= xj[2] * delz * cforce;
   -              v[3] -= xj[0] * dely * cforce;
   -              v[4] -= xj[0] * delz * cforce;
   -              v[5] -= xj[1] * delz * cforce;
   -            }
   -            if (EVFLAG) vlist[n++] = j;
   -
   -          } else {
   -            if (EVFLAG) key += 2;
   -
   -            fdx = -delx*cforce;
   -            fdy = -dely*cforce;
   -            fdz = -delz*cforce;
   -
   -            fOx = fdx*alpha_O;
   -            fOy = fdy*alpha_O;
   -            fOz = fdz*alpha_O;
   -
   -            fHx = alpha_H * fdx;
   -            fHy = alpha_H * fdy;
   -            fHz = alpha_H * fdz;
   -
   -            f[j][0] += fOx;
   -            f[j][1] += fOy;
   -            f[j][2] += fOz;
   -
   -            f[jH1][0] += fHx;
   -            f[jH1][1] += fHy;
   -            f[jH1][2] += fHz;
   -
   -            f[jH2][0] += fHx;
   -            f[jH2][1] += fHy;
   -            f[jH2][2] += fHz;
   -
   -            if (VFLAG) {
   -              xH1 = x[jH1];
   -              xH2 = x[jH2];
   -              const double xHsumx = xH1[0] + xH2[0];
   -              const double xHsumy = xH1[1] + xH2[1];
   -              const double xHsumz = xH1[2] + xH2[2];
   -              v[0] += xj[0]*fOx + xHsumx*fHx;
   -              v[1] += xj[1]*fOy + xHsumy*fHy;
   -              v[2] += xj[2]*fOz + xHsumz*fHz;
   -              v[3] += xj[0]*fOy + xHsumx*fHy;
   -              v[4] += xj[0]*fOz + xHsumx*fHz;
   -              v[5] += xj[1]*fOz + xHsumy*fHz;
   -            }
   -            if (EVFLAG) {
   -              vlist[n++] = j;
   -              vlist[n++] = jH1;
   -              vlist[n++] = jH2;
   -            }
   +          if (EVFLAG) vlist[n++] = j;
   +
   +        } else {
   +          if (EVFLAG) key += 2;
   +
   +          fdx = -delx*cforce;
   +          fdy = -dely*cforce;
   +          fdz = -delz*cforce;
   +
   +          fOx = fdx*alpha_O;
   +          fOy = fdy*alpha_O;
   +          fOz = fdz*alpha_O;
   +
   +          fHx = alpha_H * fdx;
   +          fHy = alpha_H * fdy;
   +          fHz = alpha_H * fdz;
   +
   +          f[j][0] += fOx;
   +          f[j][1] += fOy;
   +          f[j][2] += fOz;
   +
   +          f[jH1][0] += fHx;
   +          f[jH1][1] += fHy;
   +          f[jH1][2] += fHz;
   +
   +          f[jH2][0] += fHx;
   +          f[jH2][1] += fHy;
   +          f[jH2][2] += fHz;
   +
   +          if (VFLAG) {
   +            xH1 = x[jH1];
   +            xH2 = x[jH2];
   +            const double xHsumx = xH1[0] + xH2[0];
   +            const double xHsumy = xH1[1] + xH2[1];
   +            const double xHsumz = xH1[2] + xH2[2];
   +            v[0] += xj[0]*fOx + xHsumx*fHx;
   +            v[1] += xj[1]*fOy + xHsumy*fHy;
   +            v[2] += xj[2]*fOz + xHsumz*fHz;
   +            v[3] += xj[0]*fOy + xHsumx*fHy;
   +            v[4] += xj[0]*fOz + xHsumx*fHz;
   +            v[5] += xj[1]*fOz + xHsumy*fHz;
   +          }
   +          if (EVFLAG) {
   +            vlist[n++] = j;
   +            vlist[n++] = jH1;
   +            vlist[n++] = jH2;
              }
   -
   -          if (EFLAG) {
   -            if (CTABLE || rsq <= tabinnersq)
   -              ecoul = prefactor*erfc;
   -            else {
   -              table = etable[itable] + fraction*detable[itable];
   -              ecoul = qiqj * table;
   -            }
   -            if (factor_coul < 1.0) ecoul -= (1.0-factor_coul)*prefactor;
   -          } else ecoul = 0.0;
   -          if (EVFLAG) ev_tally_tip4p(key,vlist,v,ecoul,alpha);
            }
   +
   +        if (EFLAG) {
   +          if (CTABLE || rsq <= tabinnersq)
   +            ecoul = prefactor*erfc;
   +          else {
   +            table = etable[itable] + fraction*detable[itable];
   +            ecoul = qiqj * table;
   +          }
   +          if (factor_coul < 1.0) ecoul -= (1.0-factor_coul)*prefactor;
   +        } else ecoul = 0.0;
   +        if (EVFLAG) ev_tally_tip4p(key,vlist,v,ecoul,alpha);
          }
        }
    
