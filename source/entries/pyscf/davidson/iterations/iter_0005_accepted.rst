Iteration 0005 — ccc6bedc3559 (accepted)
========================================


GitHub commit: `ccc6bedc3559 <iter-0005-page-head-ccc6bedc3559_>`_
Published branch: `fermilink-optimize/pyscf-davidson <https://github.com/skilled-scipkg/pyscf/tree/fermilink-optimize%2Fpyscf-davidson>`_

Change summary
--------------


Correct real TDDFT Davidson preconditioning to pass the full lower LR residual block by using `-R_y` in `_lr_eig.real_eig`

Acceptance rationale
--------------------


Primary metric improved 23.1% vs incumbent with correctness passing, no hard reject, and no forbidden-path edits.

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
| incumbent commit | `e159c73ee967 <iter-0005-guardrails-incumbent-e159c73ee967af1e8b1845e69ccc5649c7f5fe33_>`_ |
+------------------+--------------------------------------------------------------------------------------------+
| candidate commit | `ccc6bedc3559 <iter-0005-guardrails-candidate-ccc6bedc355910e090c9f6b2ff95f4ec41a3611c_>`_ |
+------------------+--------------------------------------------------------------------------------------------+
| incumbent metric | 58.7037                                                                                    |
+------------------+--------------------------------------------------------------------------------------------+
| candidate metric | 45.1467                                                                                    |
+------------------+--------------------------------------------------------------------------------------------+
| baseline metric  | 128.753                                                                                    |
+------------------+--------------------------------------------------------------------------------------------+
| Δ vs incumbent   | +23.094% (lower-is-better sign)                                                            |
+------------------+--------------------------------------------------------------------------------------------+
| changed files    | pyscf/tdscf/_lr_eig.py                                                                     |
+------------------+--------------------------------------------------------------------------------------------+


Diffstat
--------


.. code-block:: text

    pyscf/tdscf/_lr_eig.py | 2 +-
    1 file changed, 1 insertion(+), 1 deletion(-)

Diff
----


:download:`download full diff <_diffs/iter_0005_ccc6bedc3559.diff>`

.. code-block:: diff

   diff --git a/pyscf/tdscf/_lr_eig.py b/pyscf/tdscf/_lr_eig.py
   index 87fdadfea..e9a3b6656 100644
   --- a/pyscf/tdscf/_lr_eig.py
   +++ b/pyscf/tdscf/_lr_eig.py
   @@ -693,7 +693,7 @@ def real_eig(aop, x0, precond, tol_residual=1e-5, nroots=1, x0sym=None, pick=Non
                break
    
            r_index = r_norms > tol_residual
   -        XY_new = precond(np.vstack([R_x[:,r_index], R_y[:,r_index]]).T, w[r_index])
   +        XY_new = precond(np.vstack([R_x[:,r_index], -R_y[:,r_index]]).T, w[r_index])
            X_new = XY_new[:,:A_size].T
            Y_new = XY_new[:,A_size:].T
            if x0sym is None:


.. _iter-0005-page-head-ccc6bedc3559: https://github.com/skilled-scipkg/pyscf/commit/ccc6bedc3559
.. _iter-0005-guardrails-incumbent-e159c73ee967af1e8b1845e69ccc5649c7f5fe33: https://github.com/skilled-scipkg/pyscf/commit/e159c73ee967af1e8b1845e69ccc5649c7f5fe33
.. _iter-0005-guardrails-candidate-ccc6bedc355910e090c9f6b2ff95f4ec41a3611c: https://github.com/skilled-scipkg/pyscf/commit/ccc6bedc355910e090c9f6b2ff95f4ec41a3611c