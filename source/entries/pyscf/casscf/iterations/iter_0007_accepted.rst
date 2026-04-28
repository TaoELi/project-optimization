Iteration 0007 — 2eb53381cda9 (accepted)
========================================


GitHub commit: `2eb53381cda9 <iter-0007-page-head-2eb53381cda9_>`_
Published branch: `fermilink-optimize/pyscf-casscf <https://github.com/skilled-scipkg/pyscf/tree/fermilink-optimize%2Fpyscf-casscf>`_

Change summary
--------------


Use lazy density-fitted JK for Newton-CASSCF AH Hessian-vector response contractions while keeping exact AO2MO/CASCI energies and strict convergence checks unchanged

Acceptance rationale
--------------------


Correctness passed and the primary metric improved by ~52.65% without forbidden cache or final-answer reuse.

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
| incumbent commit | `44c83aaae41f <iter-0007-guardrails-incumbent-44c83aaae41f9222bb1f78703609c8fd9d7be454_>`_ |
+------------------+--------------------------------------------------------------------------------------------+
| candidate commit | `2eb53381cda9 <iter-0007-guardrails-candidate-2eb53381cda999c7282ed9a298fa8c5cc28b3e99_>`_ |
+------------------+--------------------------------------------------------------------------------------------+
| incumbent metric | 102.444                                                                                    |
+------------------+--------------------------------------------------------------------------------------------+
| candidate metric | 48.5098                                                                                    |
+------------------+--------------------------------------------------------------------------------------------+
| baseline metric  | 102.444                                                                                    |
+------------------+--------------------------------------------------------------------------------------------+
| Δ vs incumbent   | +52.647% (lower-is-better sign)                                                            |
+------------------+--------------------------------------------------------------------------------------------+
| changed files    | pyscf/mcscf/newton_casscf.py                                                               |
+------------------+--------------------------------------------------------------------------------------------+


Diffstat
--------


.. code-block:: text

    pyscf/mcscf/newton_casscf.py | 24 ++++++++++++++++++++++++
    1 file changed, 24 insertions(+)

Diff
----


:download:`download full diff <_diffs/iter_0007_2eb53381cda9.diff>`

.. code-block:: diff

   diff --git a/pyscf/mcscf/newton_casscf.py b/pyscf/mcscf/newton_casscf.py
   index da2c677e4..b695365d3 100644
   --- a/pyscf/mcscf/newton_casscf.py
   +++ b/pyscf/mcscf/newton_casscf.py
   @@ -29,6 +29,7 @@ from pyscf.lib import logger
    from pyscf.mcscf import casci, mc1step, addons
    from pyscf.mcscf.casci import get_fock, cas_natorb, canonicalize
    from pyscf import scf
   +from pyscf import df
    from pyscf.soscf import ciah
    
    def _pack_ci_get_H (mc, mo, ci0):
   @@ -851,6 +852,29 @@ class CASSCF(mc1step.CASSCF):
                                 e_tot, e_tot-elast, ss[0])
            return e_tot, e_cas, fcivec
    
   +    def update_jk_in_ah(self, mo, r, casdm1, eris):
   +        ncore = self.ncore
   +        ncas = self.ncas
   +        nocc = ncore + ncas
   +
   +        dm3 = reduce(numpy.dot, (mo[:,:ncore], r[:ncore,ncore:], mo[:,ncore:].T))
   +        dm3 = dm3 + dm3.T
   +        dm4 = reduce(numpy.dot, (mo[:,ncore:nocc], casdm1, r[ncore:nocc], mo.T))
   +        dm4 = dm4 + dm4.T
   +
   +        with_df = getattr(self, '_ah_jk_df', None)
   +        if with_df is None or with_df.mol is not self.mol:
   +            with_df = df.DF(self.mol)
   +            with_df.max_memory = self.max_memory
   +            with_df.stdout = self.stdout
   +            with_df.verbose = self.verbose
   +            self._ah_jk_df = with_df
   +
   +        vj, vk = with_df.get_jk((dm3, dm3*2+dm4), hermi=1)
   +        va = reduce(numpy.dot, (casdm1, mo[:,ncore:nocc].T, vj[0]*2-vk[0], mo))
   +        vc = reduce(numpy.dot, (mo[:,:ncore].T, vj[1]*2-vk[1], mo[:,ncore:]))
   +        return va, vc
   +
        def update_ao2mo(self, mo):
            raise DeprecationWarning('update_ao2mo was obsoleted since pyscf v1.0.  '
                                     'Use .ao2mo method instead')


.. _iter-0007-page-head-2eb53381cda9: https://github.com/skilled-scipkg/pyscf/commit/2eb53381cda9
.. _iter-0007-guardrails-incumbent-44c83aaae41f9222bb1f78703609c8fd9d7be454: https://github.com/skilled-scipkg/pyscf/commit/44c83aaae41f9222bb1f78703609c8fd9d7be454
.. _iter-0007-guardrails-candidate-2eb53381cda999c7282ed9a298fa8c5cc28b3e99: https://github.com/skilled-scipkg/pyscf/commit/2eb53381cda999c7282ed9a298fa8c5cc28b3e99