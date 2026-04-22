Iteration 0010 — 64da7449ef1c (accepted)
========================================


GitHub commit: `64da7449ef1c <iter-0010-page-head-64da7449ef1c_>`_
Published branch: `fermilink-optimize/pyscf-davidson <https://github.com/skilled-scipkg/pyscf/tree/fermilink-optimize%2Fpyscf-davidson>`_

Change summary
--------------


Add a 0.05 Hartree real_eig correction preconditioner spectral shift to reduce late-cycle B3LYP TDDFT response work

Acceptance rationale
--------------------


Primary metric improved 4.96% vs incumbent with correctness passing and no hard-reject or forbidden-path condition.

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
| incumbent commit | `c1e86bf88a7c <iter-0010-guardrails-incumbent-c1e86bf88a7c6a1869e07db4610aac679be74c73_>`_ |
+------------------+--------------------------------------------------------------------------------------------+
| candidate commit | `64da7449ef1c <iter-0010-guardrails-candidate-64da7449ef1ca3c6e649efc7b8dba9823eb3eff2_>`_ |
+------------------+--------------------------------------------------------------------------------------------+
| incumbent metric | 31.5845                                                                                    |
+------------------+--------------------------------------------------------------------------------------------+
| candidate metric | 30.0181                                                                                    |
+------------------+--------------------------------------------------------------------------------------------+
| baseline metric  | 128.753                                                                                    |
+------------------+--------------------------------------------------------------------------------------------+
| Δ vs incumbent   | +4.960% (lower-is-better sign)                                                             |
+------------------+--------------------------------------------------------------------------------------------+
| changed files    | pyscf/tdscf/_lr_eig.py                                                                     |
+------------------+--------------------------------------------------------------------------------------------+


Diffstat
--------


.. code-block:: text

    pyscf/tdscf/_lr_eig.py | 3 ++-
    1 file changed, 2 insertions(+), 1 deletion(-)

Diff
----


:download:`download full diff <_diffs/iter_0010_64da7449ef1c.diff>`

.. code-block:: diff

   diff --git a/pyscf/tdscf/_lr_eig.py b/pyscf/tdscf/_lr_eig.py
   index 41206dead..0937fed8f 100644
   --- a/pyscf/tdscf/_lr_eig.py
   +++ b/pyscf/tdscf/_lr_eig.py
   @@ -694,7 +694,8 @@ def real_eig(aop, x0, precond, tol_residual=1e-5, nroots=1, x0sym=None, pick=Non
                break
    
            r_index = r_norms > tol_residual
   -        XY_new = precond(np.vstack([R_x[:,r_index], -R_y[:,r_index]]).T, w[r_index])
   +        XY_new = precond(np.vstack([R_x[:,r_index], -R_y[:,r_index]]).T,
   +                         w[r_index] + 5e-2)
            X_new = XY_new[:,:A_size].T
            Y_new = XY_new[:,A_size:].T
            if x0sym is None:


.. _iter-0010-page-head-64da7449ef1c: https://github.com/skilled-scipkg/pyscf/commit/64da7449ef1c
.. _iter-0010-guardrails-incumbent-c1e86bf88a7c6a1869e07db4610aac679be74c73: https://github.com/skilled-scipkg/pyscf/commit/c1e86bf88a7c6a1869e07db4610aac679be74c73
.. _iter-0010-guardrails-candidate-64da7449ef1ca3c6e649efc7b8dba9823eb3eff2: https://github.com/skilled-scipkg/pyscf/commit/64da7449ef1ca3c6e649efc7b8dba9823eb3eff2