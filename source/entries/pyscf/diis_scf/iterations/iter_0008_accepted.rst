Iteration 0008 — a5198fab0782 (accepted)
========================================


GitHub commit: `a5198fab0782 <iter-0008-page-head-a5198fab0782_>`_
Published branch: `fermilink-optimize/pyscf-diis_scf <https://github.com/skilled-scipkg/pyscf/tree/fermilink-optimize%2Fpyscf-diis_scf>`_

Change summary
--------------


Combine AO-basis CDIIS error vectors for restricted near-zero initial HOMO-LUMO gaps with deferred per-cycle writes for the default temporary chkfile when no callback is active.

Acceptance rationale
--------------------


Correctness passed and the 2.18% incumbent improvement clears the 2% promotion threshold without persistent caching or answer replay.

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
| incumbent commit | `3fc6d11adee9 <iter-0008-guardrails-incumbent-3fc6d11adee97f0e08d83121c94db014cddc0be1_>`_ |
+------------------+--------------------------------------------------------------------------------------------+
| candidate commit | `a5198fab0782 <iter-0008-guardrails-candidate-a5198fab07828e455dcde5705ce3dd3324adc61a_>`_ |
+------------------+--------------------------------------------------------------------------------------------+
| incumbent metric | 1.13042                                                                                    |
+------------------+--------------------------------------------------------------------------------------------+
| candidate metric | 1.1058                                                                                     |
+------------------+--------------------------------------------------------------------------------------------+
| baseline metric  | 1.22637                                                                                    |
+------------------+--------------------------------------------------------------------------------------------+
| Δ vs incumbent   | +2.177% (lower-is-better sign)                                                             |
+------------------+--------------------------------------------------------------------------------------------+
| changed files    | pyscf/scf/hf.py                                                                            |
+------------------+--------------------------------------------------------------------------------------------+


Diffstat
--------


.. code-block:: text

    pyscf/scf/hf.py | 21 ++++++++++++++++++---
    1 file changed, 18 insertions(+), 3 deletions(-)

Diff
----


:download:`download full diff <_diffs/iter_0008_a5198fab0782.diff>`

.. code-block:: diff

   diff --git a/pyscf/scf/hf.py b/pyscf/scf/hf.py
   index d5d18d94d..c62c220c8 100644
   --- a/pyscf/scf/hf.py
   +++ b/pyscf/scf/hf.py
   @@ -153,7 +153,15 @@ Keyword argument "init_dm" is replaced by "dm0"''')
            # Since the ingredients for the Fock matrix has already been built, we can
            # just go ahead and use it to determine the orthonormal basis vectors.
            fock = mf.get_fock(h1e, s1e, vhf, dm)
   -        _, mf_diis.Corth = mf.eig(fock, s1e)
   +        mo_energy_diis, corth = mf.eig(fock, s1e)
   +        mo_energy_diis = numpy.asarray(mo_energy_diis)
   +        nocc = mol.nelectron // 2
   +        if (mo_energy_diis.ndim == 1 and 0 < nocc < mo_energy_diis.size and
   +                abs(numpy.sort(mo_energy_diis.real)[nocc] -
   +                    numpy.sort(mo_energy_diis.real)[nocc-1]) < 1e-8):
   +            mf_diis.Corth = None
   +        else:
   +            mf_diis.Corth = corth
        else:
            mf_diis = None
    
   @@ -161,6 +169,10 @@ Keyword argument "init_dm" is replaced by "dm0"''')
            # Explicit overwrite the mol object in chkfile
            # Note in pbc.scf, mf.mol == mf.cell, cell is saved under key "mol"
            chkfile.save_mol(mol, mf.chkfile)
   +    tmp_chkfile = getattr(getattr(mf, '_chkfile', None), 'name', None)
   +    standard_dump_chk = getattr(mf.dump_chk, '__func__', None) is SCF.dump_chk
   +    defer_chkfile = (dump_chk and mf.chkfile and not callable(callback) and
   +                     mf.chkfile == tmp_chkfile and standard_dump_chk)
    
        # A preprocessing hook before the SCF iteration
        mf.pre_kernel(locals())
   @@ -196,7 +208,7 @@ Keyword argument "init_dm" is replaced by "dm0"''')
            elif abs(e_tot-last_hf_e) < conv_tol and norm_gorb < conv_tol_grad:
                scf_conv = True
    
   -        if dump_chk and mf.chkfile:
   +        if dump_chk and mf.chkfile and not defer_chkfile:
                mf.dump_chk(locals())
    
            if callable(callback):
   @@ -245,9 +257,12 @@ Keyword argument "init_dm" is replaced by "dm0"''')
                scf_conv = False
            logger.info(mf, 'Extra cycle  E= %.15g  delta_E= %4.3g  |g|= %4.3g  |ddm|= %4.3g',
                        e_tot, e_tot-last_hf_e, norm_gorb, norm_ddm)
   -        if dump_chk and mf.chkfile:
   +        if dump_chk and mf.chkfile and not defer_chkfile:
                mf.dump_chk(locals())
    
   +    if defer_chkfile:
   +        mf.dump_chk(locals())
   +
        logger.timer(mf, 'scf_cycle', *cput0)
        # A post-processing hook before return
        mf.post_kernel(locals())


.. _iter-0008-page-head-a5198fab0782: https://github.com/skilled-scipkg/pyscf/commit/a5198fab0782
.. _iter-0008-guardrails-incumbent-3fc6d11adee97f0e08d83121c94db014cddc0be1: https://github.com/skilled-scipkg/pyscf/commit/3fc6d11adee97f0e08d83121c94db014cddc0be1
.. _iter-0008-guardrails-candidate-a5198fab07828e455dcde5705ce3dd3324adc61a: https://github.com/skilled-scipkg/pyscf/commit/a5198fab07828e455dcde5705ce3dd3324adc61a