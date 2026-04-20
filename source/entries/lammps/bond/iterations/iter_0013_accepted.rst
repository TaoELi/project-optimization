Iteration 0013 — 30584dbe2ab0 (accepted)
========================================


GitHub commit: `30584dbe2ab0 <iter-0013-page-head-30584dbe2ab0_>`_
Published branch: `fermilink-optimize/lammps-bond <https://github.com/skilled-scipkg/lammps/tree/fermilink-optimize%2Flammps-bond>`_

Change summary
--------------


Replace the one-type bond-to-angle reuse with an angle-oriented cache built in `bond_class2`, then consume it from a dedicated branch-free, two-way-unrolled `angle_harmonic` hot path with direct fallback to the incumbent local-geometry path when bond/angle pairing does not align.

Acceptance rationale
--------------------


Correctness passed and weighted_median_bond_seconds improved 3.87% versus incumbent, clearing the 2% promotion threshold.

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
| incumbent commit | `0826ffd39444 <iter-0013-guardrails-incumbent-0826ffd3944439fee85546aa3c9671ac4478beb4_>`_ |
+------------------+--------------------------------------------------------------------------------------------+
| candidate commit | `30584dbe2ab0 <iter-0013-guardrails-candidate-30584dbe2ab0c1af8b2a7a566a91b31999295371_>`_ |
+------------------+--------------------------------------------------------------------------------------------+
| incumbent metric | 0.23436                                                                                    |
+------------------+--------------------------------------------------------------------------------------------+
| candidate metric | 0.22529                                                                                    |
+------------------+--------------------------------------------------------------------------------------------+
| baseline metric  | 0.32627                                                                                    |
+------------------+--------------------------------------------------------------------------------------------+
| Δ vs incumbent   | +3.870% (lower-is-better sign)                                                             |
+------------------+--------------------------------------------------------------------------------------------+
| changed files    | src/CLASS2/bond_class2.cpp, src/CLASS2/bond_class2.h, src/MOLECULE/angle_harmonic.cpp      |
+------------------+--------------------------------------------------------------------------------------------+


Diffstat
--------


.. code-block:: text

    src/CLASS2/bond_class2.cpp      | 101 ++++++++++++++----------
    src/CLASS2/bond_class2.h        |  12 +++
    src/MOLECULE/angle_harmonic.cpp | 171 ++++++++++++++++++++++++++--------------
    3 files changed, 180 insertions(+), 104 deletions(-)

Diff
----


:download:`download full diff <_diffs/iter_0013_30584dbe2ab0.diff>`

.. code-block:: diff

   diff --git a/src/CLASS2/bond_class2.cpp b/src/CLASS2/bond_class2.cpp
   index 358375a6dc..097547c488 100644
   --- a/src/CLASS2/bond_class2.cpp
   +++ b/src/CLASS2/bond_class2.cpp
   @@ -40,17 +40,15 @@ using int3_t = struct {
      int a, b, t;
    };
    
   +using int4_t = struct {
   +  int a, b, c, t;
   +};
   +
    namespace LAMMPS_NS {
   -std::vector<int> bond_class2_hot_i1;
   -std::vector<int> bond_class2_hot_i2;
   -std::vector<double> bond_class2_hot_delx;
   -std::vector<double> bond_class2_hot_dely;
   -std::vector<double> bond_class2_hot_delz;
   -std::vector<double> bond_class2_hot_rsq;
   -std::vector<double> bond_class2_hot_r;
   -bigint bond_class2_hot_timestep = -1;
   -int bond_class2_hot_nbondlist = 0;
   -bool bond_class2_hot_valid = false;
   +std::vector<BondClass2AngleHotEntry> bond_class2_hot_angle_cache;
   +bigint bond_class2_hot_angle_timestep = -1;
   +int bond_class2_hot_angle_nanglelist = 0;
   +bool bond_class2_hot_angle_valid = false;
    }    // namespace LAMMPS_NS
    
    /* ---------------------------------------------------------------------- */
   @@ -83,7 +81,7 @@ void BondClass2::compute(int eflag, int vflag)
      int newton_bond = force->newton_bond;
      const int one_type = atom->nbondtypes == 1;
    
   -  bond_class2_hot_valid = false;
   +  bond_class2_hot_angle_valid = false;
    
      ev_init(eflag,vflag);
    
   @@ -128,20 +126,22 @@ void BondClass2::eval_one_type_hot(int nbondlist)
      const auto *_noalias const x = (dbl3_t *) atom->x[0];
      auto *_noalias const f = (dbl3_t *) atom->f[0];    // NOLINT
      const int3_t *_noalias const bondlist = (int3_t *) neighbor->bondlist[0];
   +  const int nanglelist = neighbor->nanglelist;
   +  const bool maybe_angle_cache = nbondlist == 2 * nanglelist;
   +  const int4_t *_noalias const anglelist = maybe_angle_cache ? (int4_t *) neighbor->anglelist[0] : nullptr;
      const double coeff_r0 = r0[1];
      const double coeff_2k2 = 2.0 * k2[1];
      const double coeff_3k3 = 3.0 * k3[1];
      const double coeff_4k4 = 4.0 * k4[1];
    
   -  bond_class2_hot_i1.resize(nbondlist);
   -  bond_class2_hot_i2.resize(nbondlist);
   -  bond_class2_hot_delx.resize(nbondlist);
   -  bond_class2_hot_dely.resize(nbondlist);
   -  bond_class2_hot_delz.resize(nbondlist);
   -  bond_class2_hot_rsq.resize(nbondlist);
   -  bond_class2_hot_r.resize(nbondlist);
   -  bond_class2_hot_nbondlist = nbondlist;
   -  bond_class2_hot_timestep = update->ntimestep;
   +  bool angle_cache_valid = maybe_angle_cache;
   +  if (maybe_angle_cache) {
   +    bond_class2_hot_angle_cache.resize(nanglelist);
   +    bond_class2_hot_angle_nanglelist = nanglelist;
   +    bond_class2_hot_angle_timestep = update->ntimestep;
   +  } else {
   +    bond_class2_hot_angle_nanglelist = 0;
   +  }
    
      int n = 0;
    
   @@ -164,20 +164,43 @@ void BondClass2::eval_one_type_hot(int nbondlist)
        const double r0 = sqrt(rsq0);
        const double r1 = sqrt(rsq1);
    
   -    bond_class2_hot_i1[n] = i1_0;
   -    bond_class2_hot_i2[n] = i2_0;
   -    bond_class2_hot_delx[n] = delx0;
   -    bond_class2_hot_dely[n] = dely0;
   -    bond_class2_hot_delz[n] = delz0;
   -    bond_class2_hot_rsq[n] = rsq0;
   -    bond_class2_hot_r[n] = r0;
   -    bond_class2_hot_i1[n + 1] = i1_1;
   -    bond_class2_hot_i2[n + 1] = i2_1;
   -    bond_class2_hot_delx[n + 1] = delx1;
   -    bond_class2_hot_dely[n + 1] = dely1;
   -    bond_class2_hot_delz[n + 1] = delz1;
   -    bond_class2_hot_rsq[n + 1] = rsq1;
   -    bond_class2_hot_r[n + 1] = r1;
   +    if (angle_cache_valid) {
   +      const int angle_index = n >> 1;
   +      const int angle_i1 = anglelist[angle_index].a;
   +      const int angle_i2 = anglelist[angle_index].b;
   +      const int angle_i3 = anglelist[angle_index].c;
   +
   +      if (i1_0 == angle_i2 && i1_1 == angle_i2) {
   +        auto &cache = bond_class2_hot_angle_cache[angle_index];
   +        if (i2_0 == angle_i1 && i2_1 == angle_i3) {
   +          cache.delx1 = -delx0;
   +          cache.dely1 = -dely0;
   +          cache.delz1 = -delz0;
   +          cache.rsq1 = rsq0;
   +          cache.r1 = r0;
   +          cache.delx2 = -delx1;
   +          cache.dely2 = -dely1;
   +          cache.delz2 = -delz1;
   +          cache.rsq2 = rsq1;
   +          cache.r2 = r1;
   +        } else if (i2_0 == angle_i3 && i2_1 == angle_i1) {
   +          cache.delx1 = -delx1;
   +          cache.dely1 = -dely1;
   +          cache.delz1 = -delz1;
   +          cache.rsq1 = rsq1;
   +          cache.r1 = r1;
   +          cache.delx2 = -delx0;
   +          cache.dely2 = -dely0;
   +          cache.delz2 = -delz0;
   +          cache.rsq2 = rsq0;
   +          cache.r2 = r0;
   +        } else {
   +          angle_cache_valid = false;
   +        }
   +      } else {
   +        angle_cache_valid = false;
   +      }
   +    }
    
        const double dr0 = r0 - coeff_r0;
        const double dr1 = r1 - coeff_r0;
   @@ -226,14 +249,6 @@ void BondClass2::eval_one_type_hot(int nbondlist)
        const double rsq = delx * delx + dely * dely + delz * delz;
        const double r = sqrt(rsq);
    
   -    bond_class2_hot_i1[n] = i1;
   -    bond_class2_hot_i2[n] = i2;
   -    bond_class2_hot_delx[n] = delx;
   -    bond_class2_hot_dely[n] = dely;
   -    bond_class2_hot_delz[n] = delz;
   -    bond_class2_hot_rsq[n] = rsq;
   -    bond_class2_hot_r[n] = r;
   -
        const double dr = r - coeff_r0;
        const double dr2 = dr * dr;
        const double dr3 = dr2 * dr;
   @@ -254,7 +269,7 @@ void BondClass2::eval_one_type_hot(int nbondlist)
        f[i2].z -= fz;
      }
    
   -  bond_class2_hot_valid = true;
   +  bond_class2_hot_angle_valid = angle_cache_valid;
    }
    
    template <int EVFLAG, int EFLAG, int NEWTON_BOND>
   diff --git a/src/CLASS2/bond_class2.h b/src/CLASS2/bond_class2.h
   index d8ae006de7..6bcddc6e9e 100644
   --- a/src/CLASS2/bond_class2.h
   +++ b/src/CLASS2/bond_class2.h
   @@ -22,8 +22,20 @@ BondStyle(class2,BondClass2);
    
    #include "bond.h"
    
   +#include <vector>
   +
    namespace LAMMPS_NS {
    
   +struct BondClass2AngleHotEntry {
   +  double delx1, dely1, delz1, rsq1, r1;
   +  double delx2, dely2, delz2, rsq2, r2;
   +};
   +
   +extern std::vector<BondClass2AngleHotEntry> bond_class2_hot_angle_cache;
   +extern bigint bond_class2_hot_angle_timestep;
   +extern int bond_class2_hot_angle_nanglelist;
   +extern bool bond_class2_hot_angle_valid;
   +
    class BondClass2 : public Bond {
     public:
      BondClass2(class LAMMPS *);
   diff --git a/src/MOLECULE/angle_harmonic.cpp b/src/MOLECULE/angle_harmonic.cpp
   index af1e8d9451..0e42dc22eb 100644
   --- a/src/MOLECULE/angle_harmonic.cpp
   +++ b/src/MOLECULE/angle_harmonic.cpp
   @@ -13,6 +13,7 @@
    
    #include "angle_harmonic.h"
    
   +#include "../CLASS2/bond_class2.h"
    #include "atom.h"
    #include "comm.h"
    #include "domain.h"
   @@ -41,19 +42,6 @@ using int4_t = struct {
      int a, b, c, t;
    };
    
   -namespace LAMMPS_NS {
   -extern std::vector<int> bond_class2_hot_i1;
   -extern std::vector<int> bond_class2_hot_i2;
   -extern std::vector<double> bond_class2_hot_delx;
   -extern std::vector<double> bond_class2_hot_dely;
   -extern std::vector<double> bond_class2_hot_delz;
   -extern std::vector<double> bond_class2_hot_rsq;
   -extern std::vector<double> bond_class2_hot_r;
   -extern bigint bond_class2_hot_timestep;
   -extern int bond_class2_hot_nbondlist;
   -extern bool bond_class2_hot_valid;
   -}    // namespace LAMMPS_NS
   -
    /* ---------------------------------------------------------------------- */
    
    AngleHarmonic::AngleHarmonic(LAMMPS *_lmp) : Angle(_lmp)
   @@ -128,58 +116,119 @@ void AngleHarmonic::eval_one_type_hot(int nanglelist)
      const double coeff_k = k[1];
      const double coeff_theta0 = theta0[1];
    
   -  if (bond_class2_hot_valid && bond_class2_hot_timestep == update->ntimestep &&
   -      bond_class2_hot_nbondlist == 2 * nanglelist) {
   -    for (int n = 0; n < nanglelist; n++) {
   +  if (bond_class2_hot_angle_valid && bond_class2_hot_angle_timestep == update->ntimestep &&
   +      bond_class2_hot_angle_nanglelist == nanglelist) {
   +    int n = 0;
   +    for (; n + 1 < nanglelist; n += 2) {
   +      const int i1_0 = anglelist[n].a;
   +      const int i2_0 = anglelist[n].b;
   +      const int i3_0 = anglelist[n].c;
   +      const int i1_1 = anglelist[n + 1].a;
   +      const int i2_1 = anglelist[n + 1].b;
   +      const int i3_1 = anglelist[n + 1].c;
   +      const auto &_noalias cache0 = bond_class2_hot_angle_cache[n];
   +      const auto &_noalias cache1 = bond_class2_hot_angle_cache[n + 1];
   +
   +      const double delx1_0 = cache0.delx1;
   +      const double dely1_0 = cache0.dely1;
   +      const double delz1_0 = cache0.delz1;
   +      const double delx2_0 = cache0.delx2;
   +      const double dely2_0 = cache0.dely2;
   +      const double delz2_0 = cache0.delz2;
   +      const double rsq1_0 = cache0.rsq1;
   +      const double rsq2_0 = cache0.rsq2;
   +      const double r1_0 = cache0.r1;
   +      const double r2_0 = cache0.r2;
   +
   +      double c0 = delx1_0 * delx2_0 + dely1_0 * dely2_0 + delz1_0 * delz2_0;
   +      c0 /= r1_0 * r2_0;
   +      if (c0 > 1.0) c0 = 1.0;
   +      if (c0 < -1.0) c0 = -1.0;
   +      double s0 = sqrt(1.0 - c0 * c0);
   +      if (s0 < SMALL) s0 = SMALL;
   +      s0 = 1.0 / s0;
   +      const double dtheta0 = acos(c0) - coeff_theta0;
   +      const double tk0 = coeff_k * dtheta0;
   +      const double a0 = -2.0 * tk0 * s0;
   +      const double a11_0 = a0 * c0 / rsq1_0;
   +      const double a12_0 = -a0 / (r1_0 * r2_0);
   +      const double a22_0 = a0 * c0 / rsq2_0;
   +      const double f1x0 = a11_0 * delx1_0 + a12_0 * delx2_0;
   +      const double f1y0 = a11_0 * dely1_0 + a12_0 * dely2_0;
   +      const double f1z0 = a11_0 * delz1_0 + a12_0 * delz2_0;
   +      const double f3x0 = a22_0 * delx2_0 + a12_0 * delx1_0;
   +      const double f3y0 = a22_0 * dely2_0 + a12_0 * dely1_0;
   +      const double f3z0 = a22_0 * delz2_0 + a12_0 * delz1_0;
   +      f[i1_0].x += f1x0;
   +      f[i1_0].y += f1y0;
   +      f[i1_0].z += f1z0;
   +      f[i2_0].x -= f1x0 + f3x0;
   +      f[i2_0].y -= f1y0 + f3y0;
   +      f[i2_0].z -= f1z0 + f3z0;
   +      f[i3_0].x += f3x0;
   +      f[i3_0].y += f3y0;
   +      f[i3_0].z += f3z0;
   +
   +      const double delx1_1 = cache1.delx1;
   +      const double dely1_1 = cache1.dely1;
   +      const double delz1_1 = cache1.delz1;
   +      const double delx2_1 = cache1.delx2;
   +      const double dely2_1 = cache1.dely2;
   +      const double delz2_1 = cache1.delz2;
   +      const double rsq1_1 = cache1.rsq1;
   +      const double rsq2_1 = cache1.rsq2;
   +      const double r1_1 = cache1.r1;
   +      const double r2_1 = cache1.r2;
   +
   +      double c1 = delx1_1 * delx2_1 + dely1_1 * dely2_1 + delz1_1 * delz2_1;
   +      c1 /= r1_1 * r2_1;
   +      if (c1 > 1.0) c1 = 1.0;
   +      if (c1 < -1.0) c1 = -1.0;
   +      double s1 = sqrt(1.0 - c1 * c1);
   +      if (s1 < SMALL) s1 = SMALL;
   +      s1 = 1.0 / s1;
   +      const double dtheta1 = acos(c1) - coeff_theta0;
   +      const double tk1 = coeff_k * dtheta1;
   +      const double a1 = -2.0 * tk1 * s1;
   +      const double a11_1 = a1 * c1 / rsq1_1;
   +      const double a12_1 = -a1 / (r1_1 * r2_1);
   +      const double a22_1 = a1 * c1 / rsq2_1;
   +      const double f1x1 = a11_1 * delx1_1 + a12_1 * delx2_1;
   +      const double f1y1 = a11_1 * dely1_1 + a12_1 * dely2_1;
   +      const double f1z1 = a11_1 * delz1_1 + a12_1 * delz2_1;
   +      const double f3x1 = a22_1 * delx2_1 + a12_1 * delx1_1;
   +      const double f3y1 = a22_1 * dely2_1 + a12_1 * dely1_1;
   +      const double f3z1 = a22_1 * delz2_1 + a12_1 * delz1_1;
   +      f[i1_1].x += f1x1;
   +      f[i1_1].y += f1y1;
   +      f[i1_1].z += f1z1;
   +      f[i2_1].x -= f1x1 + f3x1;
   +      f[i2_1].y -= f1y1 + f3y1;
   +      f[i2_1].z -= f1z1 + f3z1;
   +      f[i3_1].x += f3x1;
   +      f[i3_1].y += f3y1;
   +      f[i3_1].z += f3z1;
   +    }
   +
   +    for (; n < nanglelist; n++) {
          const int i1 = anglelist[n].a;
          const int i2 = anglelist[n].b;
          const int i3 = anglelist[n].c;
   -      const int b0 = 2 * n;
   -      const int b1 = b0 + 1;
   -
   -      double delx1 = 0.0, dely1 = 0.0, delz1 = 0.0;
   -      double delx2 = 0.0, dely2 = 0.0, delz2 = 0.0;
   -      double rsq1 = 0.0, rsq2 = 0.0, r1 = 0.0, r2 = 0.0;
   -      if (bond_class2_hot_i1[b0] == i2 && bond_class2_hot_i2[b0] == i1 && bond_class2_hot_i1[b1] == i2 &&
   -          bond_class2_hot_i2[b1] == i3) {
   -        delx1 = -bond_class2_hot_delx[b0];
   -        dely1 = -bond_class2_hot_dely[b0];
   -        delz1 = -bond_class2_hot_delz[b0];
   -        rsq1 = bond_class2_hot_rsq[b0];
   -        r1 = bond_class2_hot_r[b0];
   -        delx2 = -bond_class2_hot_delx[b1];
   -        dely2 = -bond_class2_hot_dely[b1];
   -        delz2 = -bond_class2_hot_delz[b1];
   -        rsq2 = bond_class2_hot_rsq[b1];
   -        r2 = bond_class2_hot_r[b1];
   -      } else if (bond_class2_hot_i1[b0] == i2 && bond_class2_hot_i2[b0] == i3 &&
   -                 bond_class2_hot_i1[b1] == i2 && bond_class2_hot_i2[b1] == i1) {
   -        delx1 = -bond_class2_hot_delx[b1];
   -        dely1 = -bond_class2_hot_dely[b1];
   -        delz1 = -bond_class2_hot_delz[b1];
   -        rsq1 = bond_class2_hot_rsq[b1];
   -        r1 = bond_class2_hot_r[b1];
   -        delx2 = -bond_class2_hot_delx[b0];
   -        dely2 = -bond_class2_hot_dely[b0];
   -        delz2 = -bond_class2_hot_delz[b0];
   -        rsq2 = bond_class2_hot_rsq[b0];
   -        r2 = bond_class2_hot_r[b0];
   -      } else {
   -        delx1 = x[i1].x - x[i2].x;
   -        dely1 = x[i1].y - x[i2].y;
   -        delz1 = x[i1].z - x[i2].z;
   -        delx2 = x[i3].x - x[i2].x;
   -        dely2 = x[i3].y - x[i2].y;
   -        delz2 = x[i3].z - x[i2].z;
   -        rsq1 = delx1 * delx1 + dely1 * dely1 + delz1 * delz1;
   -        rsq2 = delx2 * delx2 + dely2 * dely2 + delz2 * delz2;
   -        r1 = sqrt(rsq1);
   -        r2 = sqrt(rsq2);
   -      }
   +      const auto &_noalias cache = bond_class2_hot_angle_cache[n];
   +
   +      const double delx1 = cache.delx1;
   +      const double dely1 = cache.dely1;
   +      const double delz1 = cache.delz1;
   +      const double delx2 = cache.delx2;
   +      const double dely2 = cache.dely2;
   +      const double delz2 = cache.delz2;
   +      const double rsq1 = cache.rsq1;
   +      const double rsq2 = cache.rsq2;
   +      const double r1 = cache.r1;
   +      const double r2 = cache.r2;
    
   -      const double r12 = r1 * r2;
          double c = delx1 * delx2 + dely1 * dely2 + delz1 * delz2;
   -      c /= r12;
   +      c /= r1 * r2;
          if (c > 1.0) c = 1.0;
          if (c < -1.0) c = -1.0;
    
   @@ -191,7 +240,7 @@ void AngleHarmonic::eval_one_type_hot(int nanglelist)
          const double tk = coeff_k * dtheta;
          const double a = -2.0 * tk * s;
          const double a11 = a * c / rsq1;
   -      const double a12 = -a / r12;
   +      const double a12 = -a / (r1 * r2);
          const double a22 = a * c / rsq2;
    
          const double f1x = a11 * delx1 + a12 * delx2;


.. _iter-0013-page-head-30584dbe2ab0: https://github.com/skilled-scipkg/lammps/commit/30584dbe2ab0
.. _iter-0013-guardrails-incumbent-0826ffd3944439fee85546aa3c9671ac4478beb4: https://github.com/skilled-scipkg/lammps/commit/0826ffd3944439fee85546aa3c9671ac4478beb4
.. _iter-0013-guardrails-candidate-30584dbe2ab0c1af8b2a7a566a91b31999295371: https://github.com/skilled-scipkg/lammps/commit/30584dbe2ab0c1af8b2a7a566a91b31999295371