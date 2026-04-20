Iteration 0025 — bf07f28c4118 (accepted)
========================================


GitHub commit: `bf07f28c4118 <iter-0025-page-head-bf07f28c4118_>`_
Published branch: `fermilink-optimize/lammps-bond <https://github.com/skilled-scipkg/lammps/tree/fermilink-optimize%2Flammps-bond>`_

Change summary
--------------


Compact the hot angle cache to exact r12 payloads and widen the one-type cached bond_class2 producer to four-way exact unrolling while preserving the incumbent angle_harmonic arithmetic order and fallback guard

Acceptance rationale
--------------------


Correctness passed and weighted_median_bond_seconds improved from 0.22529 to 0.21875, a 2.90% gain over the incumbent that clears the 2% promotion threshold.

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
| incumbent commit | `30584dbe2ab0 <iter-0025-guardrails-incumbent-30584dbe2ab0c1af8b2a7a566a91b31999295371_>`_ |
+------------------+--------------------------------------------------------------------------------------------+
| candidate commit | `bf07f28c4118 <iter-0025-guardrails-candidate-bf07f28c411869f7ae2ab685092e6fe0372178f5_>`_ |
+------------------+--------------------------------------------------------------------------------------------+
| incumbent metric | 0.22529                                                                                    |
+------------------+--------------------------------------------------------------------------------------------+
| candidate metric | 0.21875                                                                                    |
+------------------+--------------------------------------------------------------------------------------------+
| baseline metric  | 0.32627                                                                                    |
+------------------+--------------------------------------------------------------------------------------------+
| Δ vs incumbent   | +2.903% (lower-is-better sign)                                                             |
+------------------+--------------------------------------------------------------------------------------------+
| changed files    | src/CLASS2/bond_class2.cpp, src/CLASS2/bond_class2.h, src/MOLECULE/angle_harmonic.cpp      |
+------------------+--------------------------------------------------------------------------------------------+


Diffstat
--------


.. code-block:: text

    src/CLASS2/bond_class2.cpp      | 191 ++++++++++++++++++++++++++++++++--------
    src/CLASS2/bond_class2.h        |   4 +-
    src/MOLECULE/angle_harmonic.cpp |  21 ++---
    3 files changed, 167 insertions(+), 49 deletions(-)

Diff
----


:download:`download full diff <_diffs/iter_0025_bf07f28c4118.diff>`

.. code-block:: diff

   diff --git a/src/CLASS2/bond_class2.cpp b/src/CLASS2/bond_class2.cpp
   index 097547c488..4250318473 100644
   --- a/src/CLASS2/bond_class2.cpp
   +++ b/src/CLASS2/bond_class2.cpp
   @@ -44,6 +44,44 @@ using int4_t = struct {
      int a, b, c, t;
    };
    
   +static inline __attribute__((always_inline)) bool fill_hot_angle_cache_entry(
   +    BondClass2AngleHotEntry &cache, const int4_t &angle, const int i1_0, const int i2_0,
   +    const int i1_1, const int i2_1, const double delx0, const double dely0, const double delz0,
   +    const double rsq0, const double r0, const double delx1, const double dely1,
   +    const double delz1, const double rsq1, const double r1)
   +{
   +  if (i1_0 != angle.b || i1_1 != angle.b) return false;
   +
   +  const double r12 = r0 * r1;
   +  if (i2_0 == angle.a && i2_1 == angle.c) {
   +    cache.delx1 = -delx0;
   +    cache.dely1 = -dely0;
   +    cache.delz1 = -delz0;
   +    cache.rsq1 = rsq0;
   +    cache.delx2 = -delx1;
   +    cache.dely2 = -dely1;
   +    cache.delz2 = -delz1;
   +    cache.rsq2 = rsq1;
   +    cache.r12 = r12;
   +    return true;
   +  }
   +
   +  if (i2_0 == angle.c && i2_1 == angle.a) {
   +    cache.delx1 = -delx1;
   +    cache.dely1 = -dely1;
   +    cache.delz1 = -delz1;
   +    cache.rsq1 = rsq1;
   +    cache.delx2 = -delx0;
   +    cache.dely2 = -dely0;
   +    cache.delz2 = -delz0;
   +    cache.rsq2 = rsq0;
   +    cache.r12 = r12;
   +    return true;
   +  }
   +
   +  return false;
   +}
   +
    namespace LAMMPS_NS {
    std::vector<BondClass2AngleHotEntry> bond_class2_hot_angle_cache;
    bigint bond_class2_hot_angle_timestep = -1;
   @@ -145,7 +183,120 @@ void BondClass2::eval_one_type_hot(int nbondlist)
    
      int n = 0;
    
   -  // Two-way unrolling gives the compiler more independent exact work around the sqrt/div pair.
   +  // Four-way unrolling exposes more independent exact work around the sqrt/div pair.
   +  for (; n + 3 < nbondlist; n += 4) {
   +    const int i1_0 = bondlist[n].a;
   +    const int i2_0 = bondlist[n].b;
   +    const int i1_1 = bondlist[n + 1].a;
   +    const int i2_1 = bondlist[n + 1].b;
   +    const int i1_2 = bondlist[n + 2].a;
   +    const int i2_2 = bondlist[n + 2].b;
   +    const int i1_3 = bondlist[n + 3].a;
   +    const int i2_3 = bondlist[n + 3].b;
   +
   +    const double delx0 = x[i1_0].x - x[i2_0].x;
   +    const double dely0 = x[i1_0].y - x[i2_0].y;
   +    const double delz0 = x[i1_0].z - x[i2_0].z;
   +    const double delx1 = x[i1_1].x - x[i2_1].x;
   +    const double dely1 = x[i1_1].y - x[i2_1].y;
   +    const double delz1 = x[i1_1].z - x[i2_1].z;
   +    const double delx2 = x[i1_2].x - x[i2_2].x;
   +    const double dely2 = x[i1_2].y - x[i2_2].y;
   +    const double delz2 = x[i1_2].z - x[i2_2].z;
   +    const double delx3 = x[i1_3].x - x[i2_3].x;
   +    const double dely3 = x[i1_3].y - x[i2_3].y;
   +    const double delz3 = x[i1_3].z - x[i2_3].z;
   +
   +    const double rsq0 = delx0 * delx0 + dely0 * dely0 + delz0 * delz0;
   +    const double rsq1 = delx1 * delx1 + dely1 * dely1 + delz1 * delz1;
   +    const double rsq2 = delx2 * delx2 + dely2 * dely2 + delz2 * delz2;
   +    const double rsq3 = delx3 * delx3 + dely3 * dely3 + delz3 * delz3;
   +    const double r0 = sqrt(rsq0);
   +    const double r1 = sqrt(rsq1);
   +    const double r2 = sqrt(rsq2);
   +    const double r3 = sqrt(rsq3);
   +
   +    if (angle_cache_valid) {
   +      const int angle_index0 = n >> 1;
   +      const int angle_index1 = angle_index0 + 1;
   +      angle_cache_valid =
   +          fill_hot_angle_cache_entry(bond_class2_hot_angle_cache[angle_index0], anglelist[angle_index0],
   +                                     i1_0, i2_0, i1_1, i2_1, delx0, dely0, delz0, rsq0, r0, delx1,
   +                                     dely1, delz1, rsq1, r1) &&
   +          fill_hot_angle_cache_entry(bond_class2_hot_angle_cache[angle_index1], anglelist[angle_index1],
   +                                     i1_2, i2_2, i1_3, i2_3, delx2, dely2, delz2, rsq2, r2, delx3,
   +                                     dely3, delz3, rsq3, r3);
   +    }
   +
   +    const double dr0 = r0 - coeff_r0;
   +    const double dr1 = r1 - coeff_r0;
   +    const double dr2 = r2 - coeff_r0;
   +    const double dr3 = r3 - coeff_r0;
   +    const double dr20 = dr0 * dr0;
   +    const double dr21 = dr1 * dr1;
   +    const double dr22 = dr2 * dr2;
   +    const double dr23 = dr3 * dr3;
   +    const double dr30 = dr20 * dr0;
   +    const double dr31 = dr21 * dr1;
   +    const double dr32 = dr22 * dr2;
   +    const double dr33 = dr23 * dr3;
   +
   +    const double de_bond0 = coeff_2k2 * dr0 + coeff_3k3 * dr20 + coeff_4k4 * dr30;
   +    const double de_bond1 = coeff_2k2 * dr1 + coeff_3k3 * dr21 + coeff_4k4 * dr31;
   +    const double de_bond2 = coeff_2k2 * dr2 + coeff_3k3 * dr22 + coeff_4k4 * dr32;
   +    const double de_bond3 = coeff_2k2 * dr3 + coeff_3k3 * dr23 + coeff_4k4 * dr33;
   +    double fbond0, fbond1, fbond2, fbond3;
   +    if (r0 > 0.0) fbond0 = -de_bond0 / r0;
   +    else fbond0 = 0.0;
   +    if (r1 > 0.0) fbond1 = -de_bond1 / r1;
   +    else fbond1 = 0.0;
   +    if (r2 > 0.0) fbond2 = -de_bond2 / r2;
   +    else fbond2 = 0.0;
   +    if (r3 > 0.0) fbond3 = -de_bond3 / r3;
   +    else fbond3 = 0.0;
   +
   +    const double fx0 = delx0 * fbond0;
   +    const double fy0 = dely0 * fbond0;
   +    const double fz0 = delz0 * fbond0;
   +    f[i1_0].x += fx0;
   +    f[i1_0].y += fy0;
   +    f[i1_0].z += fz0;
   +    f[i2_0].x -= fx0;
   +    f[i2_0].y -= fy0;
   +    f[i2_0].z -= fz0;
   +
   +    const double fx1 = delx1 * fbond1;
   +    const double fy1 = dely1 * fbond1;
   +    const double fz1 = delz1 * fbond1;
   +    f[i1_1].x += fx1;
   +    f[i1_1].y += fy1;
   +    f[i1_1].z += fz1;
   +    f[i2_1].x -= fx1;
   +    f[i2_1].y -= fy1;
   +    f[i2_1].z -= fz1;
   +
   +    const double fx2 = delx2 * fbond2;
   +    const double fy2 = dely2 * fbond2;
   +    const double fz2 = delz2 * fbond2;
   +    f[i1_2].x += fx2;
   +    f[i1_2].y += fy2;
   +    f[i1_2].z += fz2;
   +    f[i2_2].x -= fx2;
   +    f[i2_2].y -= fy2;
   +    f[i2_2].z -= fz2;
   +
   +    const double fx3 = delx3 * fbond3;
   +    const double fy3 = dely3 * fbond3;
   +    const double fz3 = delz3 * fbond3;
   +    f[i1_3].x += fx3;
   +    f[i1_3].y += fy3;
   +    f[i1_3].z += fz3;
   +    f[i2_3].x -= fx3;
   +    f[i2_3].y -= fy3;
   +    f[i2_3].z -= fz3;
   +  }
   +
   +  // Two-way unrolling still handles the final aligned bond pair exactly.
      for (; n + 1 < nbondlist; n += 2) {
        const int i1_0 = bondlist[n].a;
        const int i2_0 = bondlist[n].b;
   @@ -166,40 +317,10 @@ void BondClass2::eval_one_type_hot(int nbondlist)
    
        if (angle_cache_valid) {
          const int angle_index = n >> 1;
   -      const int angle_i1 = anglelist[angle_index].a;
   -      const int angle_i2 = anglelist[angle_index].b;
   -      const int angle_i3 = anglelist[angle_index].c;
   -
   -      if (i1_0 == angle_i2 && i1_1 == angle_i2) {
   -        auto &cache = bond_class2_hot_angle_cache[angle_index];
   -        if (i2_0 == angle_i1 && i2_1 == angle_i3) {
   -          cache.delx1 = -delx0;
   -          cache.dely1 = -dely0;
   -          cache.delz1 = -delz0;
   -          cache.rsq1 = rsq0;
   -          cache.r1 = r0;
   -          cache.delx2 = -delx1;
   -          cache.dely2 = -dely1;
   -          cache.delz2 = -delz1;
   -          cache.rsq2 = rsq1;
   -          cache.r2 = r1;
   -        } else if (i2_0 == angle_i3 && i2_1 == angle_i1) {
   -          cache.delx1 = -delx1;
   -          cache.dely1 = -dely1;
   -          cache.delz1 = -delz1;
   -          cache.rsq1 = rsq1;
   -          cache.r1 = r1;
   -          cache.delx2 = -delx0;
   -          cache.dely2 = -dely0;
   -          cache.delz2 = -delz0;
   -          cache.rsq2 = rsq0;
   -          cache.r2 = r0;
   -        } else {
   -          angle_cache_valid = false;
   -        }
   -      } else {
   -        angle_cache_valid = false;
   -      }
   +      angle_cache_valid =
   +          fill_hot_angle_cache_entry(bond_class2_hot_angle_cache[angle_index], anglelist[angle_index],
   +                                     i1_0, i2_0, i1_1, i2_1, delx0, dely0, delz0, rsq0, r0, delx1,
   +                                     dely1, delz1, rsq1, r1);
        }
    
        const double dr0 = r0 - coeff_r0;
   diff --git a/src/CLASS2/bond_class2.h b/src/CLASS2/bond_class2.h
   index 6bcddc6e9e..3540a00c45 100644
   --- a/src/CLASS2/bond_class2.h
   +++ b/src/CLASS2/bond_class2.h
   @@ -27,8 +27,8 @@ BondStyle(class2,BondClass2);
    namespace LAMMPS_NS {
    
    struct BondClass2AngleHotEntry {
   -  double delx1, dely1, delz1, rsq1, r1;
   -  double delx2, dely2, delz2, rsq2, r2;
   +  double delx1, dely1, delz1, rsq1;
   +  double delx2, dely2, delz2, rsq2, r12;
    };
    
    extern std::vector<BondClass2AngleHotEntry> bond_class2_hot_angle_cache;
   diff --git a/src/MOLECULE/angle_harmonic.cpp b/src/MOLECULE/angle_harmonic.cpp
   index 0e42dc22eb..2c82c27a3a 100644
   --- a/src/MOLECULE/angle_harmonic.cpp
   +++ b/src/MOLECULE/angle_harmonic.cpp
   @@ -137,11 +137,10 @@ void AngleHarmonic::eval_one_type_hot(int nanglelist)
          const double delz2_0 = cache0.delz2;
          const double rsq1_0 = cache0.rsq1;
          const double rsq2_0 = cache0.rsq2;
   -      const double r1_0 = cache0.r1;
   -      const double r2_0 = cache0.r2;
   +      const double r12_0 = cache0.r12;
    
          double c0 = delx1_0 * delx2_0 + dely1_0 * dely2_0 + delz1_0 * delz2_0;
   -      c0 /= r1_0 * r2_0;
   +      c0 /= r12_0;
          if (c0 > 1.0) c0 = 1.0;
          if (c0 < -1.0) c0 = -1.0;
          double s0 = sqrt(1.0 - c0 * c0);
   @@ -151,7 +150,7 @@ void AngleHarmonic::eval_one_type_hot(int nanglelist)
          const double tk0 = coeff_k * dtheta0;
          const double a0 = -2.0 * tk0 * s0;
          const double a11_0 = a0 * c0 / rsq1_0;
   -      const double a12_0 = -a0 / (r1_0 * r2_0);
   +      const double a12_0 = -a0 / r12_0;
          const double a22_0 = a0 * c0 / rsq2_0;
          const double f1x0 = a11_0 * delx1_0 + a12_0 * delx2_0;
          const double f1y0 = a11_0 * dely1_0 + a12_0 * dely2_0;
   @@ -177,11 +176,10 @@ void AngleHarmonic::eval_one_type_hot(int nanglelist)
          const double delz2_1 = cache1.delz2;
          const double rsq1_1 = cache1.rsq1;
          const double rsq2_1 = cache1.rsq2;
   -      const double r1_1 = cache1.r1;
   -      const double r2_1 = cache1.r2;
   +      const double r12_1 = cache1.r12;
    
          double c1 = delx1_1 * delx2_1 + dely1_1 * dely2_1 + delz1_1 * delz2_1;
   -      c1 /= r1_1 * r2_1;
   +      c1 /= r12_1;
          if (c1 > 1.0) c1 = 1.0;
          if (c1 < -1.0) c1 = -1.0;
          double s1 = sqrt(1.0 - c1 * c1);
   @@ -191,7 +189,7 @@ void AngleHarmonic::eval_one_type_hot(int nanglelist)
          const double tk1 = coeff_k * dtheta1;
          const double a1 = -2.0 * tk1 * s1;
          const double a11_1 = a1 * c1 / rsq1_1;
   -      const double a12_1 = -a1 / (r1_1 * r2_1);
   +      const double a12_1 = -a1 / r12_1;
          const double a22_1 = a1 * c1 / rsq2_1;
          const double f1x1 = a11_1 * delx1_1 + a12_1 * delx2_1;
          const double f1y1 = a11_1 * dely1_1 + a12_1 * dely2_1;
   @@ -224,11 +222,10 @@ void AngleHarmonic::eval_one_type_hot(int nanglelist)
          const double delz2 = cache.delz2;
          const double rsq1 = cache.rsq1;
          const double rsq2 = cache.rsq2;
   -      const double r1 = cache.r1;
   -      const double r2 = cache.r2;
   +      const double r12 = cache.r12;
    
          double c = delx1 * delx2 + dely1 * dely2 + delz1 * delz2;
   -      c /= r1 * r2;
   +      c /= r12;
          if (c > 1.0) c = 1.0;
          if (c < -1.0) c = -1.0;
    
   @@ -240,7 +237,7 @@ void AngleHarmonic::eval_one_type_hot(int nanglelist)
          const double tk = coeff_k * dtheta;
          const double a = -2.0 * tk * s;
          const double a11 = a * c / rsq1;
   -      const double a12 = -a / (r1 * r2);
   +      const double a12 = -a / r12;
          const double a22 = a * c / rsq2;
    
          const double f1x = a11 * delx1 + a12 * delx2;


.. _iter-0025-page-head-bf07f28c4118: https://github.com/skilled-scipkg/lammps/commit/bf07f28c4118
.. _iter-0025-guardrails-incumbent-30584dbe2ab0c1af8b2a7a566a91b31999295371: https://github.com/skilled-scipkg/lammps/commit/30584dbe2ab0c1af8b2a7a566a91b31999295371
.. _iter-0025-guardrails-candidate-bf07f28c411869f7ae2ab685092e6fe0372178f5: https://github.com/skilled-scipkg/lammps/commit/bf07f28c411869f7ae2ab685092e6fe0372178f5