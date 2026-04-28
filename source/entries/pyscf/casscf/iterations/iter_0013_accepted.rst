Iteration 0013 — 911063b081d2 (accepted)
========================================


GitHub commit: `911063b081d2 <iter-0013-page-head-911063b081d2_>`_
Published branch: `fermilink-optimize/pyscf-casscf <https://github.com/skilled-scipkg/pyscf/tree/fermilink-optimize%2Fpyscf-casscf>`_

Change summary
--------------


Cache per-AO2MO low-rank DF A-side transforms and use cached materialized ppaa/papa slices to vectorize Newton-CASSCF AH H_co contractions while preserving exact AO2MO/CASCI energies and tolerances

Acceptance rationale
--------------------


Correctness passed and the 39.745700547s primary metric is ~2.53% faster than the incumbent without persistent-cache or final-answer reuse.

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
| incumbent commit | `334a40a0e17c <iter-0013-guardrails-incumbent-334a40a0e17c35677d3e6710f11c5b2238312a71_>`_ |
+------------------+--------------------------------------------------------------------------------------------+
| candidate commit | `911063b081d2 <iter-0013-guardrails-candidate-911063b081d23927561808e6a192015c338acc19_>`_ |
+------------------+--------------------------------------------------------------------------------------------+
| incumbent metric | 40.7776                                                                                    |
+------------------+--------------------------------------------------------------------------------------------+
| candidate metric | 39.7457                                                                                    |
+------------------+--------------------------------------------------------------------------------------------+
| baseline metric  | 102.444                                                                                    |
+------------------+--------------------------------------------------------------------------------------------+
| Δ vs incumbent   | +2.531% (lower-is-better sign)                                                             |
+------------------+--------------------------------------------------------------------------------------------+
| changed files    | pyscf/mcscf/newton_casscf.py                                                               |
+------------------+--------------------------------------------------------------------------------------------+


Diffstat
--------


.. code-block:: text

    pyscf/mcscf/newton_casscf.py | 124 +++++++++++++++++++++++++++++++++++++------
    1 file changed, 107 insertions(+), 17 deletions(-)

Diff
----


:download:`download full diff <_diffs/iter_0013_911063b081d2.diff>`

.. code-block:: diff

   diff --git a/pyscf/mcscf/newton_casscf.py b/pyscf/mcscf/newton_casscf.py
   index ad8a3eef7..3958495b4 100644
   --- a/pyscf/mcscf/newton_casscf.py
   +++ b/pyscf/mcscf/newton_casscf.py
   @@ -324,17 +324,32 @@ def gen_g_hop(casscf, mo, ci0, eris, verbose=None):
            tdm1 = tdm1 + tdm1.T
            tdm2 = tdm2 + tdm2.transpose(1,0,3,2)
            tdm2 =(tdm2 + tdm2.transpose(2,3,0,1)) * .5
   -        vhf_a = numpy.empty((nmo,ncore))
   -        paaa = numpy.empty((nmo,ncas,ncas,ncas))
   -        jk = 0
   -        for i in range(nmo):
   -            jbuf = eris.ppaa[i]
   -            kbuf = eris.papa[i]
   -            paaa[i] = jbuf[ncore:nocc]
   -            vhf_a[i] = numpy.einsum('quv,uv->q', jbuf[:ncore], tdm1)
   -            vhf_a[i]-= numpy.einsum('uqv,uv->q', kbuf[:,:ncore], tdm1) * .5
   -            jk += numpy.einsum('quv,q->uv', jbuf, ddm_c[i])
   -            jk -= numpy.einsum('uqv,q->uv', kbuf, ddm_c[i]) * .5
   +        hco_slices = _materialized_hco_slices(eris, ncore, nocc,
   +                                              casscf.max_memory)
   +        if hco_slices is not None:
   +            ppaa = eris.ppaa
   +            papa = eris.papa
   +            paaa, ppaa_core, papa_core = hco_slices
   +            vhf_a = numpy.einsum('iquv,uv->iq', ppaa_core, tdm1)
   +            vhf_a-= numpy.einsum('iuqv,uv->iq', papa_core, tdm1) * .5
   +            jk = numpy.zeros((ncas,ncas), dtype=h1e_mo.dtype)
   +            if ncore > 0:
   +                jk += numpy.einsum('iquv,iq->uv', ppaa_core, rc) * 2
   +                jk += numpy.einsum('iquv,qi->uv', ppaa[:ncore], rc) * 2
   +                jk -= numpy.einsum('iuqv,iq->uv', papa_core, rc)
   +                jk -= numpy.einsum('iuqv,qi->uv', papa[:ncore], rc)
   +        else:
   +            vhf_a = numpy.empty((nmo,ncore))
   +            paaa = numpy.empty((nmo,ncas,ncas,ncas))
   +            jk = 0
   +            for i in range(nmo):
   +                jbuf = eris.ppaa[i]
   +                kbuf = eris.papa[i]
   +                paaa[i] = jbuf[ncore:nocc]
   +                vhf_a[i] = numpy.einsum('quv,uv->q', jbuf[:ncore], tdm1)
   +                vhf_a[i]-= numpy.einsum('uqv,uv->q', kbuf[:,:ncore], tdm1) * .5
   +                jk += numpy.einsum('quv,q->uv', jbuf, ddm_c[i])
   +                jk -= numpy.einsum('uqv,q->uv', kbuf, ddm_c[i]) * .5
            g_dm2 = numpy.einsum('puwx,wxuv->pv', paaa, tdm2)
            aaaa = numpy.dot(ra.T, paaa.reshape(nmo,-1)).reshape([ncas]*4)
            aaaa = aaaa + aaaa.transpose(1,0,2,3)
   @@ -407,7 +422,39 @@ def extract_rotation(casscf, dr, u, ci0):
        return u, ci1
    
    
   -def _df_jk_from_low_rank_dm(with_df, factor_pairs):
   +def _df_lr_a_transforms(with_df, factors):
   +    nao = with_df.mol.nao_nr()
   +    arrays = []
   +    total_rank = 0
   +    for a in factors:
   +        if a.shape[1] == 0:
   +            arrays.append(None)
   +        else:
   +            a = numpy.asarray(a, order='F')
   +            arrays.append(a)
   +            total_rank += a.shape[1]
   +
   +    if total_rank == 0:
   +        return tuple(None for a in arrays)
   +
   +    naux = with_df.get_naoaux()
   +    dtype = next(a.dtype for a in arrays if a is not None)
   +    size_mb = naux * nao * total_rank * dtype.itemsize / 1e6
   +    avail_mb = getattr(with_df, 'max_memory', 0) - lib.current_memory()[0]
   +    if size_mb > max(200, avail_mb * .25):
   +        return tuple(None for a in arrays)
   +
   +    out = [[] if a is not None else None for a in arrays]
   +    for eri1 in with_df.loop():
   +        eri2 = lib.unpack_tril(eri1).reshape(-1, nao)
   +        naux_blk = eri1.shape[0]
   +        for i, a in enumerate(arrays):
   +            if a is not None:
   +                out[i].append(lib.dot(eri2, a).reshape(naux_blk, nao, a.shape[1]))
   +    return tuple(out)
   +
   +
   +def _df_jk_from_low_rank_dm(with_df, factor_pairs, a_transforms=None):
        '''DF J/K for symmetric densities D = A B^T + B A^T.'''
        nao = with_df.mol.nao_nr()
        nset = len(factor_pairs)
   @@ -423,7 +470,10 @@ def _df_jk_from_low_rank_dm(with_df, factor_pairs):
                pairs.append((numpy.asarray(a, order='F'),
                              numpy.asarray(b, order='F')))
    
   -    for eri1 in with_df.loop():
   +    if a_transforms is None:
   +        a_transforms = (None,) * nset
   +
   +    for iblk, eri1 in enumerate(with_df.loop()):
            eri = lib.unpack_tril(eri1)
            eri2 = eri.reshape(-1, nao)
            rho = numpy.zeros((nset, eri1.shape[0]))
   @@ -431,7 +481,13 @@ def _df_jk_from_low_rank_dm(with_df, factor_pairs):
                if pair is None:
                    continue
                a, b = pair
   -            la = lib.dot(eri2, a).reshape(-1, nao, a.shape[1])
   +            la_blocks = a_transforms[i]
   +            if la_blocks is None or iblk >= len(la_blocks):
   +                la = lib.dot(eri2, a).reshape(-1, nao, a.shape[1])
   +            else:
   +                la = la_blocks[iblk]
   +                if la.shape[0] != eri1.shape[0]:
   +                    la = lib.dot(eri2, a).reshape(-1, nao, a.shape[1])
                lb = lib.dot(eri2, b).reshape(-1, nao, b.shape[1])
                rho[i] = numpy.einsum('ix,pix->p', a, lb) * 2
                vk1 = lib.einsum('pix,pjx->ij', la, lb)
   @@ -460,6 +516,32 @@ def _materialize_eris_ppaa_papa(eris, max_memory):
        return eris
    
    
   +def _materialized_hco_slices(eris, ncore, nocc, max_memory):
   +    ppaa = eris.ppaa
   +    papa = eris.papa
   +    if not (isinstance(ppaa, numpy.ndarray) and isinstance(papa, numpy.ndarray)):
   +        return None
   +
   +    cache = getattr(eris, '_ah_hco_slices', None)
   +    if (cache is not None and cache[0] is ppaa and cache[1] is papa and
   +        cache[2] == ncore and cache[3] == nocc):
   +        return cache[4:]
   +
   +    paaa = ppaa[:,ncore:nocc]
   +    ppaa_core = ppaa[:,:ncore]
   +    papa_core = papa[:,:,:ncore]
   +    size_mb = (paaa.size * paaa.dtype.itemsize +
   +               ppaa_core.size * ppaa_core.dtype.itemsize +
   +               papa_core.size * papa_core.dtype.itemsize) / 1e6
   +    if lib.current_memory()[0] + size_mb < max_memory * .9:
   +        paaa = numpy.asarray(paaa, order='C')
   +        ppaa_core = numpy.asarray(ppaa_core, order='C')
   +        papa_core = numpy.asarray(papa_core, order='C')
   +
   +    eris._ah_hco_slices = (ppaa, papa, ncore, nocc, paaa, ppaa_core, papa_core)
   +    return paaa, ppaa_core, papa_core
   +
   +
    def update_orb_ci(casscf, mo, ci0, eris, x0_guess=None,
                      conv_tol_grad=1e-4, max_stepsize=None, verbose=None):
        log = logger.new_logger(casscf, verbose)
   @@ -923,17 +1005,25 @@ class CASSCF(mc1step.CASSCF):
                with_df.verbose = self.verbose
                self._ah_jk_df = with_df
    
   -        a3 = mo[:,:ncore]
            b3 = numpy.dot(mo[:,ncore:], r[:ncore,ncore:].T)
   -        a4 = numpy.dot(mo[:,ncore:nocc], casdm1)
            b4 = numpy.dot(mo, r[ncore:nocc].T)
    
            if mo.shape[0] > (ncore + ncas) * 4:
   +            a_cache = getattr(eris, '_ah_jk_lr_a_cache', None)
   +            if a_cache is None:
   +                a3 = mo[:,:ncore]
   +                a4 = numpy.dot(mo[:,ncore:nocc], casdm1)
   +                a_cache = (a3, a4,
   +                           _df_lr_a_transforms(with_df, (a3, a4)))
   +                eris._ah_jk_lr_a_cache = a_cache
   +            a3, a4, a_transforms = a_cache
                vj_lr, vk_lr = _df_jk_from_low_rank_dm(
   -                with_df, ((a3, b3), (a4, b4)))
   +                with_df, ((a3, b3), (a4, b4)), a_transforms)
                vj = numpy.asarray((vj_lr[0], vj_lr[0]*2 + vj_lr[1]))
                vk = numpy.asarray((vk_lr[0], vk_lr[0]*2 + vk_lr[1]))
            else:
   +            a3 = mo[:,:ncore]
   +            a4 = numpy.dot(mo[:,ncore:nocc], casdm1)
                dm3 = numpy.dot(a3, b3.T)
                dm3 = dm3 + dm3.T
                dm4 = numpy.dot(a4, b4.T)


.. _iter-0013-page-head-911063b081d2: https://github.com/skilled-scipkg/pyscf/commit/911063b081d2
.. _iter-0013-guardrails-incumbent-334a40a0e17c35677d3e6710f11c5b2238312a71: https://github.com/skilled-scipkg/pyscf/commit/334a40a0e17c35677d3e6710f11c5b2238312a71
.. _iter-0013-guardrails-candidate-911063b081d23927561808e6a192015c338acc19: https://github.com/skilled-scipkg/pyscf/commit/911063b081d23927561808e6a192015c338acc19