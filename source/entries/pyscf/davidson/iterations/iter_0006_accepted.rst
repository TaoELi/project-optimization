Iteration 0006 — 7612d6ba4e26 (accepted)
========================================


GitHub commit: `7612d6ba4e26 <iter-0006-page-head-7612d6ba4e26_>`_
Published branch: `fermilink-optimize/pyscf-davidson <https://github.com/skilled-scipkg/pyscf/tree/fermilink-optimize%2Fpyscf-davidson>`_

Change summary
--------------


Cap symmetric `_lr_eig.eigh` Davidson expansion to requested roots to avoid non-target Casida residual/preconditioner work

Acceptance rationale
--------------------


Primary metric improved 17.1% vs incumbent with correctness passing and no hard-reject or forbidden-path condition.

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
| incumbent commit | `ccc6bedc3559 <iter-0006-guardrails-incumbent-ccc6bedc355910e090c9f6b2ff95f4ec41a3611c_>`_ |
+------------------+--------------------------------------------------------------------------------------------+
| candidate commit | `7612d6ba4e26 <iter-0006-guardrails-candidate-7612d6ba4e26c619e3d8b59a8862dea35e62a00d_>`_ |
+------------------+--------------------------------------------------------------------------------------------+
| incumbent metric | 45.1467                                                                                    |
+------------------+--------------------------------------------------------------------------------------------+
| candidate metric | 37.4126                                                                                    |
+------------------+--------------------------------------------------------------------------------------------+
| baseline metric  | 128.753                                                                                    |
+------------------+--------------------------------------------------------------------------------------------+
| Δ vs incumbent   | +17.131% (lower-is-better sign)                                                            |
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


:download:`download full diff <_diffs/iter_0006_7612d6ba4e26.diff>`

.. code-block:: diff

   diff --git a/pyscf/tdscf/_lr_eig.py b/pyscf/tdscf/_lr_eig.py
   index e9a3b6656..328c17319 100644
   --- a/pyscf/tdscf/_lr_eig.py
   +++ b/pyscf/tdscf/_lr_eig.py
   @@ -81,8 +81,9 @@ def eigh(aop, x0, precond, tol_residual=1e-5, lindep=1e-12, nroots=1,
        if MAX_SPACE_INC is None:
            space_inc = nroots
        else:
   -        # Adding too many trial bases in each iteration may cause larger errors
   -        space_inc = max(nroots, min(MAX_SPACE_INC, x0_size//2))
   +        # Keep TD response expansion focused on requested roots instead of
   +        # spending residual/preconditioner work on non-target Ritz vectors.
   +        space_inc = nroots
    
        max_space = int(max_memory*1e6/8/x0_size / 2 - nroots - space_inc)
        if max_space < nroots * 4 < x0_size:


.. _iter-0006-page-head-7612d6ba4e26: https://github.com/skilled-scipkg/pyscf/commit/7612d6ba4e26
.. _iter-0006-guardrails-incumbent-ccc6bedc355910e090c9f6b2ff95f4ec41a3611c: https://github.com/skilled-scipkg/pyscf/commit/ccc6bedc355910e090c9f6b2ff95f4ec41a3611c
.. _iter-0006-guardrails-candidate-7612d6ba4e26c619e3d8b59a8862dea35e62a00d: https://github.com/skilled-scipkg/pyscf/commit/7612d6ba4e26c619e3d8b59a8862dea35e62a00d