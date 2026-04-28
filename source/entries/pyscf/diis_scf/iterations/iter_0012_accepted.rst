Iteration 0012 — 89eac5279c1a (accepted)
========================================


GitHub commit: `89eac5279c1a <iter-0012-page-head-89eac5279c1a_>`_
Published branch: `fermilink-optimize/pyscf-diis_scf <https://github.com/skilled-scipkg/pyscf/tree/fermilink-optimize%2Fpyscf-diis_scf>`_

Change summary
--------------


Cache reusable DFT numerical-integration AO block-loop outputs for unchanged RKS/UKS molecular grids, preserving existing NumInt integration semantics while avoiding repeated AO evaluation across SCF cycles.

Acceptance rationale
--------------------


Correctness passed and the per-object AO block cache improved the primary metric by 22.20% over incumbent without persistent caching or answer replay.

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
| incumbent commit | `a5198fab0782 <iter-0012-guardrails-incumbent-a5198fab07828e455dcde5705ce3dd3324adc61a_>`_ |
+------------------+--------------------------------------------------------------------------------------------+
| candidate commit | `89eac5279c1a <iter-0012-guardrails-candidate-89eac5279c1aaa4cbf77e1e7fe59e268c3f1f23e_>`_ |
+------------------+--------------------------------------------------------------------------------------------+
| incumbent metric | 1.1058                                                                                     |
+------------------+--------------------------------------------------------------------------------------------+
| candidate metric | 0.860287                                                                                   |
+------------------+--------------------------------------------------------------------------------------------+
| baseline metric  | 1.22637                                                                                    |
+------------------+--------------------------------------------------------------------------------------------+
| Δ vs incumbent   | +22.203% (lower-is-better sign)                                                            |
+------------------+--------------------------------------------------------------------------------------------+
| changed files    | pyscf/dft/rks.py, pyscf/dft/uks.py                                                         |
+------------------+--------------------------------------------------------------------------------------------+


Diffstat
--------


.. code-block:: text

    pyscf/dft/rks.py | 63 +++++++++++++++++++++++++++++++++++++++++++++++++++++++-
    pyscf/dft/uks.py |  4 +++-
    2 files changed, 65 insertions(+), 2 deletions(-)

Diff
----


:download:`download full diff <_diffs/iter_0012_89eac5279c1a.diff>`

.. code-block:: diff

   diff --git a/pyscf/dft/rks.py b/pyscf/dft/rks.py
   index fbb2c9357..349f813ec 100644
   --- a/pyscf/dft/rks.py
   +++ b/pyscf/dft/rks.py
   @@ -34,6 +34,63 @@ from pyscf.dft import gen_grid
    from pyscf.dft import numint
    from pyscf import __config__
    
   +_AO_CACHE_MAX_MEMORY = getattr(__config__, 'dft_rks_ao_cache_max_memory', 512)
   +
   +def _cached_numint_call(ks, ni, mol, grids, max_memory, fn):
   +    coords = getattr(grids, 'coords', None)
   +    weights = getattr(grids, 'weights', None)
   +    if coords is None or weights is None or max_memory <= 0:
   +        return fn()
   +
   +    ao_cache = getattr(ks, '_numint_ao_cache', None)
   +    if ao_cache is None:
   +        ao_cache = ks._numint_ao_cache = {}
   +
   +    block_loop = ni.block_loop
   +
   +    def cached_block_loop(mol1, grids1, nao=None, deriv=0, max_memory=2000,
   +                          non0tab=None, blksize=None, buf=None):
   +        if mol1 is not mol or grids1 is not grids or non0tab is not None or buf is not None:
   +            yield from block_loop(mol1, grids1, nao, deriv, max_memory,
   +                                  non0tab, blksize, buf)
   +            return
   +
   +        if nao is None:
   +            nao = mol1.nao
   +        ngrids = grids1.coords.shape[0]
   +        comp = (deriv + 1) * (deriv + 2) * (deriv + 3) // 6
   +        cache_mb = comp * ngrids * nao * numpy.dtype('float64').itemsize / 1e6
   +        max_cache_mb = min(_AO_CACHE_MAX_MEMORY, max_memory * .25)
   +        if cache_mb > max_cache_mb:
   +            yield from block_loop(mol1, grids1, nao, deriv, max_memory,
   +                                  non0tab, blksize, buf)
   +            return
   +
   +        non0tab_key = getattr(grids1, 'non0tab', None)
   +        key = (id(mol1), id(grids1), grids1.coords.ctypes.data,
   +               grids1.weights.ctypes.data, None if non0tab_key is None else id(non0tab_key),
   +               grids1.coords.shape, grids1.weights.shape, nao, deriv, blksize)
   +        blocks = ao_cache.get(key)
   +        if blocks is None:
   +            blocks = []
   +            for ao, mask, weight, coords in block_loop(mol1, grids1, nao, deriv,
   +                                                       max_memory, non0tab,
   +                                                       blksize, buf):
   +                if mask is not None:
   +                    mask = numpy.array(mask, copy=True)
   +                blocks.append((numpy.array(ao, copy=True), mask, weight, coords))
   +            ao_cache.clear()
   +            ao_cache[key] = blocks
   +
   +        for ao, mask, weight, coords in blocks:
   +            yield ao, mask, weight, coords
   +
   +    ni.block_loop = cached_block_loop
   +    try:
   +        return fn()
   +    finally:
   +        ni.block_loop = block_loop
   +
    def get_veff(ks, mol=None, dm=None, dm_last=0, vhf_last=0, hermi=1):
        '''Coulomb + XC functional
    
   @@ -77,7 +134,9 @@ def get_veff(ks, mol=None, dm=None, dm_last=0, vhf_last=0, hermi=1):
            n, exc, vxc = 0, 0, 0
        else:
            max_memory = ks.max_memory - lib.current_memory()[0]
   -        n, exc, vxc = ni.nr_rks(mol, ks.grids, ks.xc, dm, max_memory=max_memory)
   +        n, exc, vxc = _cached_numint_call(
   +            ks, ni, mol, ks.grids, max_memory,
   +            lambda: ni.nr_rks(mol, ks.grids, ks.xc, dm, max_memory=max_memory))
            logger.debug(ks, 'nelec by numeric integration = %s', n)
            if ks.do_nlc():
                if ni.libxc.is_nlc(ks.xc):
   @@ -343,6 +402,7 @@ class KohnShamDFT:
    ##################################################
    # don't modify the following attributes, they are not input options
            self._numint = numint.NumInt()
   +        self._numint_ao_cache = {}
    
        @property
        def omega(self):
   @@ -482,6 +542,7 @@ class KohnShamDFT:
            hf.SCF.reset(self, mol)
            self.grids.reset(mol)
            self.nlcgrids.reset(mol)
   +        self._numint_ao_cache = {}
            return self
    
        def check_sanity(self):
   diff --git a/pyscf/dft/uks.py b/pyscf/dft/uks.py
   index 3c47e4447..ad4b4ed4d 100644
   --- a/pyscf/dft/uks.py
   +++ b/pyscf/dft/uks.py
   @@ -49,7 +49,9 @@ def get_veff(ks, mol=None, dm=None, dm_last=0, vhf_last=0, hermi=1):
            n, exc, vxc = (0,0), 0, 0
        else:
            max_memory = ks.max_memory - lib.current_memory()[0]
   -        n, exc, vxc = ni.nr_uks(mol, ks.grids, ks.xc, dm, max_memory=max_memory)
   +        n, exc, vxc = rks._cached_numint_call(
   +            ks, ni, mol, ks.grids, max_memory,
   +            lambda: ni.nr_uks(mol, ks.grids, ks.xc, dm, max_memory=max_memory))
            logger.debug(ks, 'nelec by numeric integration = %s', n)
            if ks.do_nlc():
                if ni.libxc.is_nlc(ks.xc):


.. _iter-0012-page-head-89eac5279c1a: https://github.com/skilled-scipkg/pyscf/commit/89eac5279c1a
.. _iter-0012-guardrails-incumbent-a5198fab07828e455dcde5705ce3dd3324adc61a: https://github.com/skilled-scipkg/pyscf/commit/a5198fab07828e455dcde5705ce3dd3324adc61a
.. _iter-0012-guardrails-candidate-89eac5279c1aaa4cbf77e1e7fe59e268c3f1f23e: https://github.com/skilled-scipkg/pyscf/commit/89eac5279c1aaa4cbf77e1e7fe59e268c3f1f23e