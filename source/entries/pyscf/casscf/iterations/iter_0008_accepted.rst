Iteration 0008 — 7193e2b76d2b (accepted)
========================================


GitHub commit: `7193e2b76d2b <iter-0008-page-head-7193e2b76d2b_>`_
Published branch: `fermilink-optimize/pyscf-casscf <https://github.com/skilled-scipkg/pyscf/tree/fermilink-optimize%2Fpyscf-casscf>`_

Change summary
--------------


Factorize large Newton-CASSCF AH density-fitted JK response builds over low-rank dm3/dm4 factors, with dense DF-JK retained for small AO/rank cases and exact AO2MO/CASCI unchanged

Acceptance rationale
--------------------


Correctness passed, no hard reject or forbidden paths, and the primary metric improved ~14.01% versus the incumbent without persistent-cache or final-answer reuse.

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
| incumbent commit | `2eb53381cda9 <iter-0008-guardrails-incumbent-2eb53381cda999c7282ed9a298fa8c5cc28b3e99_>`_ |
+------------------+--------------------------------------------------------------------------------------------+
| candidate commit | `7193e2b76d2b <iter-0008-guardrails-candidate-7193e2b76d2badb621d1c3087410202eead80a2c_>`_ |
+------------------+--------------------------------------------------------------------------------------------+
| incumbent metric | 48.5098                                                                                    |
+------------------+--------------------------------------------------------------------------------------------+
| candidate metric | 41.7139                                                                                    |
+------------------+--------------------------------------------------------------------------------------------+
| baseline metric  | 102.444                                                                                    |
+------------------+--------------------------------------------------------------------------------------------+
| Δ vs incumbent   | +14.009% (lower-is-better sign)                                                            |
+------------------+--------------------------------------------------------------------------------------------+
| changed files    | pyscf/mcscf/newton_casscf.py                                                               |
+------------------+--------------------------------------------------------------------------------------------+


Diffstat
--------


.. code-block:: text

    pyscf/mcscf/newton_casscf.py | 58 +++++++++++++++++++++++++++++++++++++++-----
    1 file changed, 52 insertions(+), 6 deletions(-)

Diff
----


:download:`download full diff <_diffs/iter_0008_7193e2b76d2b.diff>`

.. code-block:: diff

   diff --git a/pyscf/mcscf/newton_casscf.py b/pyscf/mcscf/newton_casscf.py
   index b695365d3..70458d2bf 100644
   --- a/pyscf/mcscf/newton_casscf.py
   +++ b/pyscf/mcscf/newton_casscf.py
   @@ -406,6 +406,42 @@ def extract_rotation(casscf, dr, u, ci0):
        if nroots == 1: ci1 = ci1[0]
        return u, ci1
    
   +
   +def _df_jk_from_low_rank_dm(with_df, factor_pairs):
   +    '''DF J/K for symmetric densities D = A B^T + B A^T.'''
   +    nao = with_df.mol.nao_nr()
   +    nset = len(factor_pairs)
   +    nao_pair = nao * (nao+1) // 2
   +    vj = numpy.zeros((nset, nao_pair))
   +    vk = numpy.zeros((nset, nao, nao))
   +
   +    pairs = []
   +    for a, b in factor_pairs:
   +        if a.shape[1] == 0:
   +            pairs.append(None)
   +        else:
   +            pairs.append((numpy.asarray(a, order='F'),
   +                          numpy.asarray(b, order='F')))
   +
   +    for eri1 in with_df.loop():
   +        eri = lib.unpack_tril(eri1)
   +        eri2 = eri.reshape(-1, nao)
   +        rho = numpy.zeros((nset, eri1.shape[0]))
   +        for i, pair in enumerate(pairs):
   +            if pair is None:
   +                continue
   +            a, b = pair
   +            la = lib.dot(eri2, a).reshape(-1, nao, a.shape[1])
   +            lb = lib.dot(eri2, b).reshape(-1, nao, b.shape[1])
   +            rho[i] = numpy.einsum('ix,pix->p', a, lb) * 2
   +            vk1 = lib.einsum('pix,pjx->ij', la, lb)
   +            vk[i] += vk1 + vk1.T
   +        vj += lib.dot(rho, eri1)
   +
   +    vj = lib.unpack_tril(vj, 1).reshape(nset, nao, nao)
   +    return vj, vk
   +
   +
    def update_orb_ci(casscf, mo, ci0, eris, x0_guess=None,
                      conv_tol_grad=1e-4, max_stepsize=None, verbose=None):
        log = logger.new_logger(casscf, verbose)
   @@ -857,11 +893,6 @@ class CASSCF(mc1step.CASSCF):
            ncas = self.ncas
            nocc = ncore + ncas
    
   -        dm3 = reduce(numpy.dot, (mo[:,:ncore], r[:ncore,ncore:], mo[:,ncore:].T))
   -        dm3 = dm3 + dm3.T
   -        dm4 = reduce(numpy.dot, (mo[:,ncore:nocc], casdm1, r[ncore:nocc], mo.T))
   -        dm4 = dm4 + dm4.T
   -
            with_df = getattr(self, '_ah_jk_df', None)
            if with_df is None or with_df.mol is not self.mol:
                with_df = df.DF(self.mol)
   @@ -870,7 +901,22 @@ class CASSCF(mc1step.CASSCF):
                with_df.verbose = self.verbose
                self._ah_jk_df = with_df
    
   -        vj, vk = with_df.get_jk((dm3, dm3*2+dm4), hermi=1)
   +        a3 = mo[:,:ncore]
   +        b3 = numpy.dot(mo[:,ncore:], r[:ncore,ncore:].T)
   +        a4 = numpy.dot(mo[:,ncore:nocc], casdm1)
   +        b4 = numpy.dot(mo, r[ncore:nocc].T)
   +
   +        if mo.shape[0] > (ncore + ncas) * 4:
   +            vj_lr, vk_lr = _df_jk_from_low_rank_dm(
   +                with_df, ((a3, b3), (a4, b4)))
   +            vj = numpy.asarray((vj_lr[0], vj_lr[0]*2 + vj_lr[1]))
   +            vk = numpy.asarray((vk_lr[0], vk_lr[0]*2 + vk_lr[1]))
   +        else:
   +            dm3 = numpy.dot(a3, b3.T)
   +            dm3 = dm3 + dm3.T
   +            dm4 = numpy.dot(a4, b4.T)
   +            dm4 = dm4 + dm4.T
   +            vj, vk = with_df.get_jk((dm3, dm3*2+dm4), hermi=1)
            va = reduce(numpy.dot, (casdm1, mo[:,ncore:nocc].T, vj[0]*2-vk[0], mo))
            vc = reduce(numpy.dot, (mo[:,:ncore].T, vj[1]*2-vk[1], mo[:,ncore:]))
            return va, vc


.. _iter-0008-page-head-7193e2b76d2b: https://github.com/skilled-scipkg/pyscf/commit/7193e2b76d2b
.. _iter-0008-guardrails-incumbent-2eb53381cda999c7282ed9a298fa8c5cc28b3e99: https://github.com/skilled-scipkg/pyscf/commit/2eb53381cda999c7282ed9a298fa8c5cc28b3e99
.. _iter-0008-guardrails-candidate-7193e2b76d2badb621d1c3087410202eead80a2c: https://github.com/skilled-scipkg/pyscf/commit/7193e2b76d2badb621d1c3087410202eead80a2c