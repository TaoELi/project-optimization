Iteration 0007 — c1e86bf88a7c (accepted)
========================================


GitHub commit: `c1e86bf88a7c <iter-0007-page-head-c1e86bf88a7c_>`_
Published branch: `fermilink-optimize/pyscf-davidson <https://github.com/skilled-scipkg/pyscf/tree/fermilink-optimize%2Fpyscf-davidson>`_

Change summary
--------------


Limit symmetric `_lr_eig.eigh` residual and preconditioner candidate generation to requested target roots

Acceptance rationale
--------------------


Primary metric improved 15.6% vs incumbent with correctness passing and no hard-reject or forbidden-path condition.

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
| incumbent commit | `7612d6ba4e26 <iter-0007-guardrails-incumbent-7612d6ba4e26c619e3d8b59a8862dea35e62a00d_>`_ |
+------------------+--------------------------------------------------------------------------------------------+
| candidate commit | `c1e86bf88a7c <iter-0007-guardrails-candidate-c1e86bf88a7c6a1869e07db4610aac679be74c73_>`_ |
+------------------+--------------------------------------------------------------------------------------------+
| incumbent metric | 37.4126                                                                                    |
+------------------+--------------------------------------------------------------------------------------------+
| candidate metric | 31.5845                                                                                    |
+------------------+--------------------------------------------------------------------------------------------+
| baseline metric  | 128.753                                                                                    |
+------------------+--------------------------------------------------------------------------------------------+
| Δ vs incumbent   | +15.578% (lower-is-better sign)                                                            |
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


:download:`download full diff <_diffs/iter_0007_c1e86bf88a7c.diff>`

.. code-block:: diff

   diff --git a/pyscf/tdscf/_lr_eig.py b/pyscf/tdscf/_lr_eig.py
   index 328c17319..41206dead 100644
   --- a/pyscf/tdscf/_lr_eig.py
   +++ b/pyscf/tdscf/_lr_eig.py
   @@ -174,7 +174,7 @@ def eigh(aop, x0, precond, tol_residual=1e-5, lindep=1e-12, nroots=1,
                if x0sym is not None:
                    xs_ir = xs_ir[:row0]
    
   -        t_size = max(nroots, max_space-len(xs))
   +        t_size = min(nroots, len(xs))
            xt = -w[:t_size,None] * xs[:t_size]
            xt += ax[:t_size]
            if x0sym is not None:


.. _iter-0007-page-head-c1e86bf88a7c: https://github.com/skilled-scipkg/pyscf/commit/c1e86bf88a7c
.. _iter-0007-guardrails-incumbent-7612d6ba4e26c619e3d8b59a8862dea35e62a00d: https://github.com/skilled-scipkg/pyscf/commit/7612d6ba4e26c619e3d8b59a8862dea35e62a00d
.. _iter-0007-guardrails-candidate-c1e86bf88a7c6a1869e07db4610aac679be74c73: https://github.com/skilled-scipkg/pyscf/commit/c1e86bf88a7c6a1869e07db4610aac679be74c73