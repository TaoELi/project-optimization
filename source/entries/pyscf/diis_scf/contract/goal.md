# Optimization Goal

## Package
pyscf

## Language
python

## Target
Optimize the molecular Kohn-Sham DFT SCF loop in PySCF, with primary focus on
DIIS-family acceleration and early-cycle stabilization in the conventional
`mf.kernel()` path.

Focus on the shared molecular RKS / UKS hot path:
- `pyscf/scf/hf.py` -- generic SCF kernel, DIIS setup, Fock update call sites, and convergence checks
- `pyscf/scf/diis.py` -- CDIIS / ADIIS / EDIIS logic, error-vector construction, and switching behavior
- `pyscf/lib/diis.py` -- DIIS subspace storage, conditioning, extrapolation, and rollback behavior
- `pyscf/dft/rks.py`, `pyscf/dft/uks.py` -- DFT effective-potential path and class wiring
- `pyscf/scf/uhf.py` only when a change is required by the same shared DIIS / SCF machinery

In scope:
- more reliable DIIS subspace rejection or rollback when the extrapolation problem is ill-conditioned
- better CDIIS / ADIIS / EDIIS handoff in early or oscillatory cycles
- adaptive use of already-supported damping / level-shift controls without changing their public semantics
- reductions in SCF cycle count or DIIS-phase cost on hard but incumbent-converged molecular DFT cases

Out of scope:
- switching workloads to Newton / SOSCF
- loosening `conv_tol` / `conv_tol_grad` or increasing `max_cycle`
- changing molecule, basis, XC functional, grid level, charge, spin, occupations, symmetry, or initial guess
- speedups from unrelated code outside the SCF / DIIS path

## Editable Scope
- pyscf/lib/diis.py
- pyscf/scf/diis.py
- pyscf/scf/hf.py
- pyscf/scf/uhf.py
- pyscf/dft/rks.py
- pyscf/dft/uks.py

## Performance Metric
Minimize `weighted_median_scf_kernel_seconds`, defined as the weighted median of
per-case end-to-end `mf.kernel()` wall-clock time across all representative
workloads under pinned single-thread execution.

Secondary objectives:
- lower `scf_cycles`
- lower `diis_update_seconds`, `get_fock_seconds`, and `eig_seconds` when the runner can expose them
- no loss of convergence on any representative case

## Correctness Constraints
- All representative workloads must converge in the incumbent baseline and in accepted candidates under the stated `conv_tol`, `conv_tol_grad`, and `max_cycle` settings.
- Total SCF energy absolute delta <= `5e-8` Hartree for RKS cases and <= `1e-7` Hartree for UKS cases versus the incumbent baseline.
- Final orbital-gradient norm must remain <= the workload `conv_tol_grad` threshold when the runner exposes it.
- Molecular-orbital energies RMS delta for occupied orbitals and the lowest 10 virtual orbitals <= `2e-5` Hartree when the runner exposes comparable orbitals.
- Density-matrix RMS delta <= `2e-6` versus the incumbent baseline when the runner exposes density matrices.
- Preserve user-facing semantics for `mf.diis`, `mf.DIIS`, `diis_space`, `diis_start_cycle`, `diis_space_rollback`, `diis_damp`, `damp`, `level_shift`, `conv_tol`, `conv_tol_grad`, and `max_cycle`.
- Do not disable DIIS, silently enable Newton / SOSCF / smearing / fractional occupations, reduce DFT grid quality, or increase thread count to gain speed.
- Easy regression cases must remain converged with no more than 1 additional SCF cycle versus the incumbent baseline.
- No case-specific shortcuts keyed on molecule identity, charge, spin, basis, XC functional, or whether a case is `train-` or `test-`.

## Representative Workloads
## Representative Workloads
- train-rks-h4-square: H4 square, side length 1.70 Angstrom, charge=0, spin=0, RKS, `xc='pbe0'`, `basis='cc-pVDZ'`, `init_guess='minao'`, `grids.level=3`, `conv_tol=1e-9`, `max_cycle=80`, `diis_space=8`, no smearing.
- train-rks-stretched-n2: N2 at bond length 2.40 Angstrom, charge=0, spin=0, RKS, `xc='b3lyp'`, `basis='def2-SVP'`, `init_guess='minao'`, `grids.level=3`, `conv_tol=1e-9`, `max_cycle=80`, `diis_space=8`, no smearing.
- test-uks-no2-radical: bent NO2 radical, N-O=1.20 Angstrom, O-N-O angle=134 degrees, charge=0, spin=1, UKS, `xc='b3lyp'`, `basis='6-31g*'`, `init_guess='minao'`, `grids.level=3`, `conv_tol=1e-9`, `max_cycle=80`, `diis_space=8`.
- test-uks-feo: linear FeO, Fe-O=1.62 Angstrom, charge=0, spin=4, UKS, `xc='pbe0'`, `basis='def2-SVP'`, `init_guess='atom'`, `grids.level=3`, `conv_tol=1e-8`, `max_cycle=100`, `diis_space=10`.
- test-rks-benzene-anion-diffuse: benzene radical anion, standard planar D6h geometry, charge=-1, spin=1, UKS, `xc='b3lyp'`, `basis='6-31+g*'`, `init_guess='minao'`, `grids.level=3`, `conv_tol=1e-9`, `max_cycle=100`, `diis_space=10`.
- test-rks-h6-ring: H6 regular hexagon, side length 1.45 Angstrom, charge=0, spin=0, RKS, `xc='pbe0'`, `basis='cc-pVDZ'`, `init_guess='minao'`, `grids.level=3`, `conv_tol=1e-9`, `max_cycle=80`, `diis_space=8`.
- test-rks-stretched-co: CO at bond length 2.30 Angstrom, charge=0, spin=0, RKS, `xc='b3lyp'`, `basis='def2-SVP'`, `init_guess='minao'`, `grids.level=3`, `conv_tol=1e-9`, `max_cycle=80`, `diis_space=8`.
- test-uks-o2-dimer: two O2 fragments separated by 3.20 Angstrom, each O-O=1.21 Angstrom, total charge=0, spin=4, UKS, `xc='pbe'`, `basis='def2-SVP'`, `init_guess='atom'`, `grids.level=3`, `conv_tol=1e-8`, `max_cycle=100`, `diis_space=10`.

## Build
```bash
export SOURCE_REPO_ROOT="$(cd "$(git rev-parse --git-common-dir)/.." && pwd)"
export VENV="/anvil/scratch/x-tli22/fermilink_optimize/project_pyscf/venvs/fermilink-optimize/pyscf-diis_scf"
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
- Keep benchmark behavior deterministic with `OMP_NUM_THREADS=1`, `MKL_NUM_THREADS=1`, `OPENBLAS_NUM_THREADS=1`, and `NUMEXPR_NUM_THREADS=1`.
- Run the `## Build` commands from the PySCF repo root inside the campaign's active Python environment.
- Preserve the `train-` / `test-` case ids directly in generated benchmark cases and let FermiLink infer the split from the prefixes.
- Treat the cases above as hard but incumbent-converged workloads. Do not intentionally include baseline-nonconverged cases in the initial benchmark suite; baseline correctness must pass before the optimize campaign can start.
- If the runner can expose them, record per-case `scf_kernel_seconds`, `scf_cycles`, `diis_update_seconds`, `get_fock_seconds`, `eig_seconds`, `converged`, `e_tot`, `norm_gorb`, and density-matrix change norm.
- Assume the campaign is launched with `bin/fermilink-optimize-python` or an already activated environment. Do not hard-code a site-specific absolute venv path in generated runtime commands.
- If a case proves too fragile to pass preflight reproducibly, replace it with a nearby molecule/system of similar AO size and SCF difficulty rather than weakening the correctness gates.
