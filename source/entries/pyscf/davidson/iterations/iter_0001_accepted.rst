Iteration 0001 — b93246c4bf06 (accepted)
========================================


GitHub commit: `b93246c4bf06 <iter-0001-page-head-b93246c4bf06_>`_
Published branch: `fermilink-optimize/pyscf-davidson <https://github.com/skilled-scipkg/pyscf/tree/fermilink-optimize%2Fpyscf-davidson>`_

Change summary
--------------


Use root-specific Ritz values for LR Davidson residual preconditioning, including vectorized shifts in the shared TD preconditioner

Acceptance rationale
--------------------


22.7% primary-metric improvement with all correctness guardrails passing and no hard-reject condition.

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
| incumbent commit | `44c83aaae41f <iter-0001-guardrails-incumbent-44c83aaae41f9222bb1f78703609c8fd9d7be454_>`_ |
+------------------+--------------------------------------------------------------------------------------------+
| candidate commit | `b93246c4bf06 <iter-0001-guardrails-candidate-b93246c4bf061e853a348b1c61d56a491ac0eb6b_>`_ |
+------------------+--------------------------------------------------------------------------------------------+
| incumbent metric | 128.753                                                                                    |
+------------------+--------------------------------------------------------------------------------------------+
| candidate metric | 99.5699                                                                                    |
+------------------+--------------------------------------------------------------------------------------------+
| baseline metric  | 128.753                                                                                    |
+------------------+--------------------------------------------------------------------------------------------+
| Δ vs incumbent   | +22.666% (lower-is-better sign)                                                            |
+------------------+--------------------------------------------------------------------------------------------+
| changed files    | pyscf/tdscf/_lr_eig.py, pyscf/tdscf/rhf.py                                                 |
+------------------+--------------------------------------------------------------------------------------------+


Diffstat
--------


.. code-block:: text

    pyscf/tdscf/_lr_eig.py |  4 ++--
    pyscf/tdscf/rhf.py     | 10 ++++++++--
    2 files changed, 10 insertions(+), 4 deletions(-)

Diff
----


:download:`download full diff <_diffs/iter_0001_b93246c4bf06.diff>`

.. code-block:: diff

   diff --git a/pyscf/tdscf/_lr_eig.py b/pyscf/tdscf/_lr_eig.py
   index 1dc58e649..2acd5e230 100644
   --- a/pyscf/tdscf/_lr_eig.py
   +++ b/pyscf/tdscf/_lr_eig.py
   @@ -198,7 +198,7 @@ def eigh(aop, x0, precond, tol_residual=1e-5, lindep=1e-12, nroots=1,
            # remove subspace linear dependency
            for k, xk in enumerate(xt):
                if dx_norm[k] > tol_residual:
   -                xt[k] = precond(xk, e[0])
   +                xt[k] = precond(xk, w[k])
            xt -= xs.conj().dot(xt.T).T.dot(xs)
            xt_norm = np.linalg.norm(xt, axis=1)
    
   @@ -426,7 +426,7 @@ def eig(aop, x0, precond, tol_residual=1e-5, nroots=1, x0sym=None, pick=None,
            # remove subspace linear dependency
            for k, xk in enumerate(xt):
                if dx_norm[k] > tol_residual:
   -                xt[k] = precond(xk, e[0])
   +                xt[k] = precond(xk, w[k])
    
            xt -= xs.conj().dot(xt.T).T.dot(xs)
            c = _conj_dot(xs, xt)
   diff --git a/pyscf/tdscf/rhf.py b/pyscf/tdscf/rhf.py
   index 90b6260dc..2d33dba79 100644
   --- a/pyscf/tdscf/rhf.py
   +++ b/pyscf/tdscf/rhf.py
   @@ -875,8 +875,14 @@ class TDBase(lib.StreamObject):
        def get_precond(self, hdiag):
            def precond(x, e, *args):
                if isinstance(e, numpy.ndarray):
   -                e = e[0]
   -            diagd = hdiag - (e-self.level_shift)
   +                e = numpy.asarray(e)
   +                if e.ndim == 0 or numpy.ndim(x) == 1:
   +                    e = e.ravel()[0]
   +                    diagd = hdiag - (e-self.level_shift)
   +                else:
   +                    diagd = hdiag - (e[:,None]-self.level_shift)
   +            else:
   +                diagd = hdiag - (e-self.level_shift)
                diagd[abs(diagd)<1e-8] = 1e-8
                return x/diagd
            return precond


.. _iter-0001-page-head-b93246c4bf06: https://github.com/skilled-scipkg/pyscf/commit/b93246c4bf06
.. _iter-0001-guardrails-incumbent-44c83aaae41f9222bb1f78703609c8fd9d7be454: https://github.com/skilled-scipkg/pyscf/commit/44c83aaae41f9222bb1f78703609c8fd9d7be454
.. _iter-0001-guardrails-candidate-b93246c4bf061e853a348b1c61d56a491ac0eb6b: https://github.com/skilled-scipkg/pyscf/commit/b93246c4bf061e853a348b1c61d56a491ac0eb6b