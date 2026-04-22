Iteration 0002 — e159c73ee967 (accepted)
========================================


GitHub commit: `e159c73ee967 <iter-0002-page-head-e159c73ee967_>`_
Published branch: `fermilink-optimize/pyscf-davidson <https://github.com/skilled-scipkg/pyscf/tree/fermilink-optimize%2Fpyscf-davidson>`_

Change summary
--------------


Limit real TDDFT Davidson expansion in `_lr_eig.real_eig` to requested roots to avoid non-target residual/preconditioner work

Acceptance rationale
--------------------


Primary metric improved 41.0% vs incumbent with correctness passing and no hard-reject condition.

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
| incumbent commit | `b93246c4bf06 <iter-0002-guardrails-incumbent-b93246c4bf061e853a348b1c61d56a491ac0eb6b_>`_ |
+------------------+--------------------------------------------------------------------------------------------+
| candidate commit | `e159c73ee967 <iter-0002-guardrails-candidate-e159c73ee967af1e8b1845e69ccc5649c7f5fe33_>`_ |
+------------------+--------------------------------------------------------------------------------------------+
| incumbent metric | 99.5699                                                                                    |
+------------------+--------------------------------------------------------------------------------------------+
| candidate metric | 58.7037                                                                                    |
+------------------+--------------------------------------------------------------------------------------------+
| baseline metric  | 128.753                                                                                    |
+------------------+--------------------------------------------------------------------------------------------+
| Δ vs incumbent   | +41.043% (lower-is-better sign)                                                            |
+------------------+--------------------------------------------------------------------------------------------+
| changed files    | pyscf/tdscf/_lr_eig.py                                                                     |
+------------------+--------------------------------------------------------------------------------------------+


Diffstat
--------


.. code-block:: text

    pyscf/tdscf/_lr_eig.py | 5 +++--
    1 file changed, 3 insertions(+), 2 deletions(-)

Diff
----


:download:`download full diff <_diffs/iter_0002_e159c73ee967.diff>`

.. code-block:: diff

   diff --git a/pyscf/tdscf/_lr_eig.py b/pyscf/tdscf/_lr_eig.py
   index 2acd5e230..87fdadfea 100644
   --- a/pyscf/tdscf/_lr_eig.py
   +++ b/pyscf/tdscf/_lr_eig.py
   @@ -534,8 +534,9 @@ def real_eig(aop, x0, precond, tol_residual=1e-5, nroots=1, x0sym=None, pick=Non
        if MAX_SPACE_INC is None:
            space_inc = nroots
        else:
   -        # Adding too many trial bases in each iteration may cause larger errors
   -        space_inc = max(nroots, min(MAX_SPACE_INC, A_size//2))
   +        # Real TDDFT response batches are expensive. Keep the expansion focused
   +        # on requested roots instead of adding non-target search roots.
   +        space_inc = nroots
    
        max_space = int(max_memory*1e6/8/(4*A_size) / 2 - space_inc)
        if max_space < nroots * 4 < A_size:


.. _iter-0002-page-head-e159c73ee967: https://github.com/skilled-scipkg/pyscf/commit/e159c73ee967
.. _iter-0002-guardrails-incumbent-b93246c4bf061e853a348b1c61d56a491ac0eb6b: https://github.com/skilled-scipkg/pyscf/commit/b93246c4bf061e853a348b1c61d56a491ac0eb6b
.. _iter-0002-guardrails-candidate-e159c73ee967af1e8b1845e69ccc5649c7f5fe33: https://github.com/skilled-scipkg/pyscf/commit/e159c73ee967af1e8b1845e69ccc5649c7f5fe33