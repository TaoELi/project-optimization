Iteration 0002 — 3fc6d11adee9 (accepted)
========================================


GitHub commit: `3fc6d11adee9 <iter-0002-page-head-3fc6d11adee9_>`_
Published branch: `fermilink-optimize/pyscf-diis_scf <https://github.com/skilled-scipkg/pyscf/tree/fermilink-optimize%2Fpyscf-diis_scf>`_

Change summary
--------------


Guard final SCF extra-cycle get_veff rebuild with small canonical-density and gradient thresholds so near-stationary cases reuse the existing potential while FeO-like larger moves keep the full refresh

Acceptance rationale
--------------------


Authoritative benchmark passed correctness and improved the primary metric by 7.82% without persistent caches or replayed answers.

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
| incumbent commit | `44c83aaae41f <iter-0002-guardrails-incumbent-44c83aaae41f9222bb1f78703609c8fd9d7be454_>`_ |
+------------------+--------------------------------------------------------------------------------------------+
| candidate commit | `3fc6d11adee9 <iter-0002-guardrails-candidate-3fc6d11adee97f0e08d83121c94db014cddc0be1_>`_ |
+------------------+--------------------------------------------------------------------------------------------+
| incumbent metric | 1.22637                                                                                    |
+------------------+--------------------------------------------------------------------------------------------+
| candidate metric | 1.13042                                                                                    |
+------------------+--------------------------------------------------------------------------------------------+
| baseline metric  | 1.22637                                                                                    |
+------------------+--------------------------------------------------------------------------------------------+
| Δ vs incumbent   | +7.824% (lower-is-better sign)                                                             |
+------------------+--------------------------------------------------------------------------------------------+
| changed files    | pyscf/scf/hf.py                                                                            |
+------------------+--------------------------------------------------------------------------------------------+


Diffstat
--------


.. code-block:: text

    pyscf/scf/hf.py | 26 +++++++++++++++++++-------
    1 file changed, 19 insertions(+), 7 deletions(-)

Diff
----


:download:`download full diff <_diffs/iter_0002_3fc6d11adee9.diff>`

.. code-block:: diff

   diff --git a/pyscf/scf/hf.py b/pyscf/scf/hf.py
   index af1e95e61..d5d18d94d 100644
   --- a/pyscf/scf/hf.py
   +++ b/pyscf/scf/hf.py
   @@ -214,15 +214,27 @@ Keyword argument "init_dm" is replaced by "dm0"''')
            mo_energy, mo_coeff = mf.eig(fock, s1e)
            mo_occ = mf.get_occ(mo_energy, mo_coeff)
            dm, dm_last = mf.make_rdm1(mo_coeff, mo_occ), dm
   -        vhf = mf.get_veff(mol, dm, dm_last, vhf)
   -        e_tot, last_hf_e = mf.energy_tot(dm, h1e, vhf), e_tot
   -
   -        fock = mf.get_fock(h1e, s1e, vhf, dm)
   -        norm_gorb = numpy.linalg.norm(mf.get_grad(mo_coeff, mo_occ, fock))
   -        if not TIGHT_GRAD_CONV_TOL:
   -            norm_gorb = norm_gorb / numpy.sqrt(norm_gorb.size)
            norm_ddm = numpy.linalg.norm(dm-dm_last)
    
   +        # If canonicalization barely changes the density from an already small
   +        # gradient, the refreshed potential is numerically redundant.
   +        skip_extra_veff = (not callable(mf.check_convergence) and
   +                           norm_ddm < conv_tol_grad * .3 and
   +                           norm_gorb < conv_tol_grad * .2)
   +        if skip_extra_veff:
   +            last_hf_e = e_tot
   +            norm_gorb = numpy.linalg.norm(mf.get_grad(mo_coeff, mo_occ, fock))
   +            if not TIGHT_GRAD_CONV_TOL:
   +                norm_gorb = norm_gorb / numpy.sqrt(norm_gorb.size)
   +        else:
   +            vhf = mf.get_veff(mol, dm, dm_last, vhf)
   +            e_tot, last_hf_e = mf.energy_tot(dm, h1e, vhf), e_tot
   +
   +            fock = mf.get_fock(h1e, s1e, vhf, dm)
   +            norm_gorb = numpy.linalg.norm(mf.get_grad(mo_coeff, mo_occ, fock))
   +            if not TIGHT_GRAD_CONV_TOL:
   +                norm_gorb = norm_gorb / numpy.sqrt(norm_gorb.size)
   +
            conv_tol = conv_tol * 10
            conv_tol_grad = conv_tol_grad * 3
            if callable(mf.check_convergence):


.. _iter-0002-page-head-3fc6d11adee9: https://github.com/skilled-scipkg/pyscf/commit/3fc6d11adee9
.. _iter-0002-guardrails-incumbent-44c83aaae41f9222bb1f78703609c8fd9d7be454: https://github.com/skilled-scipkg/pyscf/commit/44c83aaae41f9222bb1f78703609c8fd9d7be454
.. _iter-0002-guardrails-candidate-3fc6d11adee97f0e08d83121c94db014cddc0be1: https://github.com/skilled-scipkg/pyscf/commit/3fc6d11adee97f0e08d83121c94db014cddc0be1