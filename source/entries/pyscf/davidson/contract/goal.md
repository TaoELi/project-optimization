# Optimization Goal

## Package
pyscf

## Language
python

## Target
Optimize the Davidson-style subspace eigensolver used by PySCF TDDFT/TDA, with primary focus on `pyscf/tdscf/_lr_eig.py` and the TD response call sites in `pyscf/tdscf/rhf.py`, `pyscf/tdscf/rks.py`, `pyscf/tdscf/uhf.py`, and `pyscf/tdscf/uks.py`.

Target optimization opportunities include:
- more efficient preconditioner strategy for reduced davidson cycles
- lower-cost projected-subspace construction and update in `eigh`, `eig`, and `real_eig`

Do not treat this as a local Python micro-optimization task. The goal is materially faster TDDFT/TDA eigensolver behavior through better Davidson/subspace algorithm choices.

## Editable Scope
- pyscf/tdscf/_lr_eig.py
- pyscf/tdscf/rhf.py
- pyscf/tdscf/rks.py
- pyscf/tdscf/uhf.py
- pyscf/tdscf/uks.py
- pyscf/lib/linalg_helper.py

## Performance Metric
Minimize end-to-end TDDFT/TDA kernel time.

Primary objective should be weighted median total wall-clock time across all benchmark cases. Secondary objective should be lower Davidson iteration count or fewer matrix-vector applications when the benchmark runner can expose those metrics.

## Correctness Constraints
- Excitation energies absolute delta <= 5e-6 Hartree vs incumbent baseline for every reported root
- Oscillator strengths absolute delta <= 1e-4 for singlet closed-shell cases where the benchmark exposes them
- Exact match of the values of transition dipole moments is not required as gauge change may flip the sign of transition dipoles
- All requested roots must converge, and root ordering should remain consistent with the incumbent baseline
- Do not loosen SCF `conv_tol`, TD solver `conv_tol`, `lindep`, `max_cycle`, `positive_eig_threshold`, `deg_eia_thresh`, `nstates`, or symmetry filtering
- Do not replace TDDFT with TDA/Casida, reduce the number of roots, change functionals/basis sets, or alter DFT grid settings to gain speed
- No case-specific shortcuts keyed on molecule identity, spin state, functional family, or whether the case is train vs test

## Representative Workloads
- train-rks-bp86-casida-benzene: benzene geometry from `examples/2-benchmark/bz.py` but with smaller basis / 6-31g / RKS / `xc='b88,p86'` / `CasidaTDDFT` / singlet / `nstates=12`
- train-rks-b3lyp-tddft-benzene: benzene geometry from `examples/2-benchmark/bz.py` but with smaller basis / 6-31g / RKS / `xc='b3lyp5'` / `TDDFT` / singlet / `nstates=10`
- train-uks-bp86-casida-allyl: allyl radical geometry from `examples/mp/12-dfump2-natorbs.py` but with smaller basis / def2-svp / spin=1 / UKS / `xc='b88,p86'` / `CasidaTDDFT` / `nstates=8`
- test-rks-bp86-casida-benzene-631gss: benzene geometry from `examples/2-benchmark/bz.py` / 6-31g** / RKS / `xc='b88,p86'` / `CasidaTDDFT` / singlet / `nstates=12`
- test-rks-b3lyp-tddft-benzene-631gss: benzene geometry from `examples/2-benchmark/bz.py` / 6-31g** / RKS / `xc='b3lyp5'` / `TDDFT` / singlet / `nstates=10`
- test-uks-bp86-casida-allyl-def2tzvp: allyl radical geometry from `examples/mp/12-dfump2-natorbs.py` / def2-TZVP / spin=1 / UKS / `xc='b88,p86'` / `CasidaTDDFT` / `nstates=8`

## Build
```bash
export SOURCE_REPO_ROOT="$(cd "$(git rev-parse --git-common-dir)/.." && pwd)"
export VENV="/anvil/scratch/x-tli22/fermilink_optimize/project_pyscf/venvs/fermilink-optimize/pyscf-davidson"
source "$VENV/bin/activate"
module remove cmake
cd pyscf/lib
mkdir -p build
cd build
cmake ..
cmake --build . -j4
cd ../../../
python -m pip install -e .
```

## Notes
- Base the benchmark setups on the larger single-machine geometries already shipped in the local PySCF tree:
  - benzene from `examples/2-benchmark/bz.py`
  - allyl radical from `examples/mp/12-dfump2-natorbs.py`
- Prefer a smaller number of materially larger cases over many toy test cases, so the benchmark is dominated by Davidson/subspace work rather than Python overhead or SCF startup noise.
- For DFT cases, mirror the upstream test setup with `dft.radi.ATOM_SPECIFIC_TREUTLER_GRIDS = False` and `mf.grids.prune = None` so the benchmark is dominated by TDDFT/TDA solver behavior instead of grid-noise differences.
- Keep benchmark behavior deterministic across repeated runs.
- If the benchmark runner can expose them, record per-case Davidson iteration count, matrix-vector application count, and total TD kernel wall time.
- Keep all workloads runnable on a single workstation-class machine with BLAS thread counts pinned to 1; prefer increasing molecular size or `nstates` only until TD kernel time clearly dominates SCF time.
- In the generated benchmark YAML, include a top-level split block:
  ```yaml
  split:
    train_case_ids:
      - train-rks-bp86-casida-benzene
      - train-rks-b3lyp-tddft-benzene
      - train-uks-bp86-casida-allyl
  ```
