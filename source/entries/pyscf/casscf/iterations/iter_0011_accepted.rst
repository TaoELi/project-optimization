Iteration 0011 — 334a40a0e17c (accepted)
========================================


GitHub commit: `334a40a0e17c <iter-0011-page-head-334a40a0e17c_>`_
Published branch: `fermilink-optimize/pyscf-casscf <https://github.com/skilled-scipkg/pyscf/tree/fermilink-optimize%2Fpyscf-casscf>`_

Change summary
--------------


Materialize exact Newton-CASSCF outcore ppaa/papa AO2MO datasets into guarded in-memory NumPy arrays to reduce repeated AH HDF5 row access while preserving exact CASCI tensors and convergence behavior

Acceptance rationale
--------------------


Correctness passed and the primary metric improved ~2.24% versus incumbent without persistent-cache or final-answer reuse.

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
| incumbent commit | `7193e2b76d2b <iter-0011-guardrails-incumbent-7193e2b76d2badb621d1c3087410202eead80a2c_>`_ |
+------------------+--------------------------------------------------------------------------------------------+
| candidate commit | `334a40a0e17c <iter-0011-guardrails-candidate-334a40a0e17c35677d3e6710f11c5b2238312a71_>`_ |
+------------------+--------------------------------------------------------------------------------------------+
| incumbent metric | 41.7139                                                                                    |
+------------------+--------------------------------------------------------------------------------------------+
| candidate metric | 40.7776                                                                                    |
+------------------+--------------------------------------------------------------------------------------------+
| baseline metric  | 102.444                                                                                    |
+------------------+--------------------------------------------------------------------------------------------+
| Δ vs incumbent   | +2.245% (lower-is-better sign)                                                             |
+------------------+--------------------------------------------------------------------------------------------+
| changed files    | pyscf/mcscf/newton_casscf.py                                                               |
+------------------+--------------------------------------------------------------------------------------------+


Diffstat
--------


.. code-block:: text

    pyscf/mcscf/newton_casscf.py | 22 ++++++++++++++++++++++
    1 file changed, 22 insertions(+)

Diff
----


:download:`download full diff <_diffs/iter_0011_334a40a0e17c.diff>`

.. code-block:: diff

   diff --git a/pyscf/mcscf/newton_casscf.py b/pyscf/mcscf/newton_casscf.py
   index 70458d2bf..ad8a3eef7 100644
   --- a/pyscf/mcscf/newton_casscf.py
   +++ b/pyscf/mcscf/newton_casscf.py
   @@ -442,6 +442,24 @@ def _df_jk_from_low_rank_dm(with_df, factor_pairs):
        return vj, vk
    
    
   +def _materialize_eris_ppaa_papa(eris, max_memory):
   +    ppaa = eris.ppaa
   +    papa = eris.papa
   +    if isinstance(ppaa, numpy.ndarray) and isinstance(papa, numpy.ndarray):
   +        return eris
   +
   +    try:
   +        size_mb = (ppaa.size * ppaa.dtype.itemsize +
   +                   papa.size * papa.dtype.itemsize) / 1e6
   +    except AttributeError:
   +        return eris
   +
   +    if lib.current_memory()[0] + size_mb < max_memory * .9:
   +        eris.ppaa = numpy.asarray(ppaa, order='C')
   +        eris.papa = numpy.asarray(papa, order='C')
   +    return eris
   +
   +
    def update_orb_ci(casscf, mo, ci0, eris, x0_guess=None,
                      conv_tol_grad=1e-4, max_stepsize=None, verbose=None):
        log = logger.new_logger(casscf, verbose)
   @@ -846,6 +864,10 @@ class CASSCF(mc1step.CASSCF):
        def kernel(self, mo_coeff=None, ci0=None, callback=None):
            return mc1step.CASSCF.kernel(self, mo_coeff, ci0, callback, kernel)
    
   +    def ao2mo(self, mo_coeff=None):
   +        eris = mc1step.CASSCF.ao2mo(self, mo_coeff)
   +        return _materialize_eris_ppaa_papa(eris, self.max_memory)
   +
        def casci(self, mo_coeff, ci0=None, eris=None, verbose=None, envs=None):
            log = logger.new_logger(self, verbose)
            fcasci = mc1step._fake_h_for_fast_casci(self, mo_coeff, eris)


.. _iter-0011-page-head-334a40a0e17c: https://github.com/skilled-scipkg/pyscf/commit/334a40a0e17c
.. _iter-0011-guardrails-incumbent-7193e2b76d2badb621d1c3087410202eead80a2c: https://github.com/skilled-scipkg/pyscf/commit/7193e2b76d2badb621d1c3087410202eead80a2c
.. _iter-0011-guardrails-candidate-334a40a0e17c35677d3e6710f11c5b2238312a71: https://github.com/skilled-scipkg/pyscf/commit/334a40a0e17c35677d3e6710f11c5b2238312a71