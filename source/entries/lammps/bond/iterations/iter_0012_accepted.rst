Iteration 0012 — 0826ffd39444 (accepted)
========================================


GitHub commit: `0826ffd39444 <iter-0012-page-head-0826ffd39444_>`_
Published branch: `fermilink-optimize/lammps-bond <https://github.com/skilled-scipkg/lammps/tree/fermilink-optimize%2Flammps-bond>`_

Change summary
--------------


Cache one-type bond_class2 bond geometry and reuse paired bond data in angle_harmonic's dominant no-evflag, Newton-on hot path, with direct fallback when local ordering does not match.

Acceptance rationale
--------------------


Correctness passed and weighted_median_bond_seconds improved from 0.24606 to 0.23436 (+4.75% vs incumbent), clearing the 2% promotion threshold.

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
| incumbent commit | `672cbe901cae <iter-0012-guardrails-incumbent-672cbe901cae5b8e7041fb6f4c3ab883fd8e6576_>`_ |
+------------------+--------------------------------------------------------------------------------------------+
| candidate commit | `0826ffd39444 <iter-0012-guardrails-candidate-0826ffd3944439fee85546aa3c9671ac4478beb4_>`_ |
+------------------+--------------------------------------------------------------------------------------------+
| incumbent metric | 0.24606                                                                                    |
+------------------+--------------------------------------------------------------------------------------------+
| candidate metric | 0.23436                                                                                    |
+------------------+--------------------------------------------------------------------------------------------+
| baseline metric  | 0.32627                                                                                    |
+------------------+--------------------------------------------------------------------------------------------+
| Δ vs incumbent   | +4.755% (lower-is-better sign)                                                             |
+------------------+--------------------------------------------------------------------------------------------+
| changed files    | src/CLASS2/bond_class2.cpp, src/MOLECULE/angle_harmonic.cpp                                |
+------------------+--------------------------------------------------------------------------------------------+


Diffstat
--------


.. code-block:: text

    src/CLASS2/bond_class2.cpp      |  54 +++++++++++++++++++++
    src/MOLECULE/angle_harmonic.cpp | 102 +++++++++++++++++++++++++++++++++++++++-
    2 files changed, 155 insertions(+), 1 deletion(-)

Diff
----


:download:`download full diff <_diffs/iter_0012_0826ffd39444.diff>`

.. code-block:: diff

   diff --git a/src/CLASS2/bond_class2.cpp b/src/CLASS2/bond_class2.cpp
   index 0984b0a7e1..358375a6dc 100644
   --- a/src/CLASS2/bond_class2.cpp
   +++ b/src/CLASS2/bond_class2.cpp
   @@ -20,6 +20,7 @@
    
    #include "atom.h"
    #include "neighbor.h"
   +#include "update.h"
    #include "comm.h"
    #include "force.h"
    #include "memory.h"
   @@ -27,6 +28,7 @@
    
    #include <cmath>
    #include <cstring>
   +#include <vector>
    
    using namespace LAMMPS_NS;
    
   @@ -38,6 +40,19 @@ using int3_t = struct {
      int a, b, t;
    };
    
   +namespace LAMMPS_NS {
   +std::vector<int> bond_class2_hot_i1;
   +std::vector<int> bond_class2_hot_i2;
   +std::vector<double> bond_class2_hot_delx;
   +std::vector<double> bond_class2_hot_dely;
   +std::vector<double> bond_class2_hot_delz;
   +std::vector<double> bond_class2_hot_rsq;
   +std::vector<double> bond_class2_hot_r;
   +bigint bond_class2_hot_timestep = -1;
   +int bond_class2_hot_nbondlist = 0;
   +bool bond_class2_hot_valid = false;
   +}    // namespace LAMMPS_NS
   +
    /* ---------------------------------------------------------------------- */
    
    BondClass2::BondClass2(LAMMPS *lmp) : Bond(lmp)
   @@ -68,6 +83,8 @@ void BondClass2::compute(int eflag, int vflag)
      int newton_bond = force->newton_bond;
      const int one_type = atom->nbondtypes == 1;
    
   +  bond_class2_hot_valid = false;
   +
      ev_init(eflag,vflag);
    
      if (nbondlist <= 0) return;
   @@ -116,6 +133,16 @@ void BondClass2::eval_one_type_hot(int nbondlist)
      const double coeff_3k3 = 3.0 * k3[1];
      const double coeff_4k4 = 4.0 * k4[1];
    
   +  bond_class2_hot_i1.resize(nbondlist);
   +  bond_class2_hot_i2.resize(nbondlist);
   +  bond_class2_hot_delx.resize(nbondlist);
   +  bond_class2_hot_dely.resize(nbondlist);
   +  bond_class2_hot_delz.resize(nbondlist);
   +  bond_class2_hot_rsq.resize(nbondlist);
   +  bond_class2_hot_r.resize(nbondlist);
   +  bond_class2_hot_nbondlist = nbondlist;
   +  bond_class2_hot_timestep = update->ntimestep;
   +
      int n = 0;
    
      // Two-way unrolling gives the compiler more independent exact work around the sqrt/div pair.
   @@ -136,6 +163,22 @@ void BondClass2::eval_one_type_hot(int nbondlist)
        const double rsq1 = delx1 * delx1 + dely1 * dely1 + delz1 * delz1;
        const double r0 = sqrt(rsq0);
        const double r1 = sqrt(rsq1);
   +
   +    bond_class2_hot_i1[n] = i1_0;
   +    bond_class2_hot_i2[n] = i2_0;
   +    bond_class2_hot_delx[n] = delx0;
   +    bond_class2_hot_dely[n] = dely0;
   +    bond_class2_hot_delz[n] = delz0;
   +    bond_class2_hot_rsq[n] = rsq0;
   +    bond_class2_hot_r[n] = r0;
   +    bond_class2_hot_i1[n + 1] = i1_1;
   +    bond_class2_hot_i2[n + 1] = i2_1;
   +    bond_class2_hot_delx[n + 1] = delx1;
   +    bond_class2_hot_dely[n + 1] = dely1;
   +    bond_class2_hot_delz[n + 1] = delz1;
   +    bond_class2_hot_rsq[n + 1] = rsq1;
   +    bond_class2_hot_r[n + 1] = r1;
   +
        const double dr0 = r0 - coeff_r0;
        const double dr1 = r1 - coeff_r0;
        const double dr20 = dr0 * dr0;
   @@ -182,6 +225,15 @@ void BondClass2::eval_one_type_hot(int nbondlist)
    
        const double rsq = delx * delx + dely * dely + delz * delz;
        const double r = sqrt(rsq);
   +
   +    bond_class2_hot_i1[n] = i1;
   +    bond_class2_hot_i2[n] = i2;
   +    bond_class2_hot_delx[n] = delx;
   +    bond_class2_hot_dely[n] = dely;
   +    bond_class2_hot_delz[n] = delz;
   +    bond_class2_hot_rsq[n] = rsq;
   +    bond_class2_hot_r[n] = r;
   +
        const double dr = r - coeff_r0;
        const double dr2 = dr * dr;
        const double dr3 = dr2 * dr;
   @@ -201,6 +253,8 @@ void BondClass2::eval_one_type_hot(int nbondlist)
        f[i2].y -= fy;
        f[i2].z -= fz;
      }
   +
   +  bond_class2_hot_valid = true;
    }
    
    template <int EVFLAG, int EFLAG, int NEWTON_BOND>
   diff --git a/src/MOLECULE/angle_harmonic.cpp b/src/MOLECULE/angle_harmonic.cpp
   index 44458e2813..af1e8d9451 100644
   --- a/src/MOLECULE/angle_harmonic.cpp
   +++ b/src/MOLECULE/angle_harmonic.cpp
   @@ -21,9 +21,11 @@
    #include "math_const.h"
    #include "memory.h"
    #include "neighbor.h"
   +#include "update.h"
    
    #include <cmath>
    #include <cstring>
   +#include <vector>
    
    using namespace LAMMPS_NS;
    using MathConst::DEG2RAD;
   @@ -39,6 +41,19 @@ using int4_t = struct {
      int a, b, c, t;
    };
    
   +namespace LAMMPS_NS {
   +extern std::vector<int> bond_class2_hot_i1;
   +extern std::vector<int> bond_class2_hot_i2;
   +extern std::vector<double> bond_class2_hot_delx;
   +extern std::vector<double> bond_class2_hot_dely;
   +extern std::vector<double> bond_class2_hot_delz;
   +extern std::vector<double> bond_class2_hot_rsq;
   +extern std::vector<double> bond_class2_hot_r;
   +extern bigint bond_class2_hot_timestep;
   +extern int bond_class2_hot_nbondlist;
   +extern bool bond_class2_hot_valid;
   +}    // namespace LAMMPS_NS
   +
    /* ---------------------------------------------------------------------- */
    
    AngleHarmonic::AngleHarmonic(LAMMPS *_lmp) : Angle(_lmp)
   @@ -113,8 +128,93 @@ void AngleHarmonic::eval_one_type_hot(int nanglelist)
      const double coeff_k = k[1];
      const double coeff_theta0 = theta0[1];
    
   -  int n = 0;
   +  if (bond_class2_hot_valid && bond_class2_hot_timestep == update->ntimestep &&
   +      bond_class2_hot_nbondlist == 2 * nanglelist) {
   +    for (int n = 0; n < nanglelist; n++) {
   +      const int i1 = anglelist[n].a;
   +      const int i2 = anglelist[n].b;
   +      const int i3 = anglelist[n].c;
   +      const int b0 = 2 * n;
   +      const int b1 = b0 + 1;
   +
   +      double delx1 = 0.0, dely1 = 0.0, delz1 = 0.0;
   +      double delx2 = 0.0, dely2 = 0.0, delz2 = 0.0;
   +      double rsq1 = 0.0, rsq2 = 0.0, r1 = 0.0, r2 = 0.0;
   +      if (bond_class2_hot_i1[b0] == i2 && bond_class2_hot_i2[b0] == i1 && bond_class2_hot_i1[b1] == i2 &&
   +          bond_class2_hot_i2[b1] == i3) {
   +        delx1 = -bond_class2_hot_delx[b0];
   +        dely1 = -bond_class2_hot_dely[b0];
   +        delz1 = -bond_class2_hot_delz[b0];
   +        rsq1 = bond_class2_hot_rsq[b0];
   +        r1 = bond_class2_hot_r[b0];
   +        delx2 = -bond_class2_hot_delx[b1];
   +        dely2 = -bond_class2_hot_dely[b1];
   +        delz2 = -bond_class2_hot_delz[b1];
   +        rsq2 = bond_class2_hot_rsq[b1];
   +        r2 = bond_class2_hot_r[b1];
   +      } else if (bond_class2_hot_i1[b0] == i2 && bond_class2_hot_i2[b0] == i3 &&
   +                 bond_class2_hot_i1[b1] == i2 && bond_class2_hot_i2[b1] == i1) {
   +        delx1 = -bond_class2_hot_delx[b1];
   +        dely1 = -bond_class2_hot_dely[b1];
   +        delz1 = -bond_class2_hot_delz[b1];
   +        rsq1 = bond_class2_hot_rsq[b1];
   +        r1 = bond_class2_hot_r[b1];
   +        delx2 = -bond_class2_hot_delx[b0];
   +        dely2 = -bond_class2_hot_dely[b0];
   +        delz2 = -bond_class2_hot_delz[b0];
   +        rsq2 = bond_class2_hot_rsq[b0];
   +        r2 = bond_class2_hot_r[b0];
   +      } else {
   +        delx1 = x[i1].x - x[i2].x;
   +        dely1 = x[i1].y - x[i2].y;
   +        delz1 = x[i1].z - x[i2].z;
   +        delx2 = x[i3].x - x[i2].x;
   +        dely2 = x[i3].y - x[i2].y;
   +        delz2 = x[i3].z - x[i2].z;
   +        rsq1 = delx1 * delx1 + dely1 * dely1 + delz1 * delz1;
   +        rsq2 = delx2 * delx2 + dely2 * dely2 + delz2 * delz2;
   +        r1 = sqrt(rsq1);
   +        r2 = sqrt(rsq2);
   +      }
    
   +      const double r12 = r1 * r2;
   +      double c = delx1 * delx2 + dely1 * dely2 + delz1 * delz2;
   +      c /= r12;
   +      if (c > 1.0) c = 1.0;
   +      if (c < -1.0) c = -1.0;
   +
   +      double s = sqrt(1.0 - c * c);
   +      if (s < SMALL) s = SMALL;
   +      s = 1.0 / s;
   +
   +      const double dtheta = acos(c) - coeff_theta0;
   +      const double tk = coeff_k * dtheta;
   +      const double a = -2.0 * tk * s;
   +      const double a11 = a * c / rsq1;
   +      const double a12 = -a / r12;
   +      const double a22 = a * c / rsq2;
   +
   +      const double f1x = a11 * delx1 + a12 * delx2;
   +      const double f1y = a11 * dely1 + a12 * dely2;
   +      const double f1z = a11 * delz1 + a12 * delz2;
   +      const double f3x = a22 * delx2 + a12 * delx1;
   +      const double f3y = a22 * dely2 + a12 * dely1;
   +      const double f3z = a22 * delz2 + a12 * delz1;
   +
   +      f[i1].x += f1x;
   +      f[i1].y += f1y;
   +      f[i1].z += f1z;
   +      f[i2].x -= f1x + f3x;
   +      f[i2].y -= f1y + f3y;
   +      f[i2].z -= f1z + f3z;
   +      f[i3].x += f3x;
   +      f[i3].y += f3y;
   +      f[i3].z += f3z;
   +    }
   +    return;
   +  }
   +
   +  int n = 0;
      // Two-way unrolling exposes more independent exact work around the sqrt/acos latency.
      for (; n + 1 < nanglelist; n += 2) {
        const int i1_0 = anglelist[n].a;


.. _iter-0012-page-head-0826ffd39444: https://github.com/skilled-scipkg/lammps/commit/0826ffd39444
.. _iter-0012-guardrails-incumbent-672cbe901cae5b8e7041fb6f4c3ab883fd8e6576: https://github.com/skilled-scipkg/lammps/commit/672cbe901cae5b8e7041fb6f4c3ab883fd8e6576
.. _iter-0012-guardrails-candidate-0826ffd3944439fee85546aa3c9671ac4478beb4: https://github.com/skilled-scipkg/lammps/commit/0826ffd3944439fee85546aa3c9671ac4478beb4