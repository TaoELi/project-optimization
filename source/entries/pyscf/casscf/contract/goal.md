# Optimization Goal

## Package
pyscf

## Language
python

## Target
Optimize PySCF CASSCF convergence behavior and end-to-end runtime, with primary focus on the matrix-free second-order / CIAH solver in `pyscf/mcscf/newton_casscf.py` and the one-step / two-step macro-micro iteration control flow in `pyscf/mcscf/mc1step.py`.

Prioritize algorithmic and numerical improvements that reduce wasted iterations without changing the scientific target:
- more stable augmented-Hessian / CIAH preconditioning using inexpensive diagonal, orbital-gap, active-space, or response information without forming the full orbital Hessian
- more reliable handling of near-redundant orbital rotations, tiny denominators, and indefinite local curvature
- better trust-region, step-acceptance, and keyframe-recovery behavior so difficult cases reach the same stationary point with fewer macroiterations and AH microiterations
- lower wasted CASCI, Hessian-vector, and restart work through smarter reuse of existing solver state when that does not change the mathematical problem
- safer state-specific and state-averaged root handling through overlap-aware diagnostics and root tracking

Out of scope:
- loosening CASSCF, AH, or FCI tolerances or increasing iteration limits to win the benchmark
- replacing CASSCF with CASCI, DMRG, selected CI, reduced active spaces, or molecule-specific shortcuts
- broad unrelated changes to SCF, integral generation, or non-MCSCF modules unless required by the same solver path

## Editable Scope
- pyscf/mcscf/newton_casscf.py
- pyscf/mcscf/newton_casscf_symm.py
- pyscf/mcscf/mc1step.py
- pyscf/mcscf/mc1step_symm.py
- pyscf/mcscf/addons.py
- pyscf/mcscf/casci.py
- pyscf/soscf/ciah.py
- pyscf/lib/linalg_helper.py

## Performance Metric
Minimize `weighted_median_casscf_kernel_seconds`, defined as the weighted median of per-case end-to-end solver wall time measured from the selected CASSCF call (`mc.mc1step()`, `mc.mc2step()`, or `mc.newton().kernel()`) until convergence, after molecule construction, RHF or ROHF setup, and active-orbital selection are complete.

Secondary objectives:
- lower `macro_cycles`, `micro_cycles`, `ah_seconds`, `h_op_calls`, `casci_seconds`, `ao2mo_seconds`, `keyframe_restarts`, and `rejected_steps` when the runner can expose them
- preserve convergence and root identity on every representative case

## Correctness Constraints
- All benchmark cases must converge in the incumbent baseline and in accepted candidates under the stated `conv_tol`, `conv_tol_grad`, `max_cycle_macro`, `max_cycle_micro`, and AH settings.
- Single-state final total CASSCF energy absolute delta must be <= `5e-8` Hartree versus the incumbent baseline.
- State-averaged final total energy absolute delta must be <= `1e-7` Hartree, and per-state energy absolute delta must be <= `5e-6` Hartree, when the runner exposes `e_states`.
- Final orbital-gradient norm and CI-gradient norm must remain no worse than the workload convergence thresholds when the runner exposes them.
- Preserve the same molecule, basis, charge, spin, symmetry flag, active electron count, active orbital count, state weights, target root, solver family, and initial orbital sorting or projection path as the baseline workload.
- Root matching for state-specific excited-root or state-averaged cases must use overlaps when comparable vectors are available; sorted-energy-only matching is not sufficient on crowded cases.
- Do not loosen `conv_tol`, `conv_tol_grad`, `fcisolver.conv_tol`, `max_cycle_macro`, `max_cycle_micro`, `max_stepsize`, `ah_conv_tol`, `ah_lindep`, `ah_level_shift`, `ah_start_tol`, `ah_start_cycle`, `ah_max_cycle`, `nroots`, state weights, `fix_spin_`, `wfnsym`, or `internal_rotation`.
- Do not replace Newton-CASSCF with an easier solver, reduce the number of states, alter active-space selection, disable scheduler callbacks, or change public API behavior to gain speed.
- No case-specific shortcuts keyed on molecule identity, bond length, basis, active space, root number, or whether a case is `train-` or `test-`.

## Representative Workloads
- train-newton-benzene-6-31g-cas66-pi-sort: Benzene geometry and active-orbital sort from `examples/mcscf/17-approx_orbital_hessian.py`; RHF; `basis='6-31g'`; `symmetry=True`; CAS(6o,6e); fixed active-orbital sort `[17,20,21,22,23,30]`; run `mc.newton().kernel(mo)`; optional comparison against `mcscf.approx_hessian(mcscf.CASSCF(...)).kernel(mo)`; purpose: larger AO/virtual orbital space with a small exact-FCI CAS, useful for testing matrix-free AH/preconditioner behavior without making the CI solve the bottleneck.
```python
atom = """
C -0.65830719  0.61123287 -0.00800148
C  0.73685281  0.61123287 -0.00800148
C  1.43439081  1.81898387 -0.00800148
C  0.73673681  3.02749287 -0.00920048
C -0.65808819  3.02741487 -0.00967948
C -1.35568919  1.81920887 -0.00868348
H -1.20806619 -0.34108413 -0.00755148
H  1.28636081 -0.34128013 -0.00668648
H  2.53407081  1.81906387 -0.00736748
H  1.28693681  3.97963587 -0.00925948
H -1.20821019  3.97969587 -0.01063248
H -2.45529319  1.81939187 -0.00886348
"""
```
- test-newton-benzene-ccpvtz-cas66-pi-sort: Benzene geometry and active-orbital sort from `examples/mcscf/17-approx_orbital_hessian.py`; RHF; `basis='ccpvtz'`; `symmetry=True`; CAS(6o,6e); fixed active-orbital sort `[17,20,21,22,23,30]`; run `mc.newton().kernel(mo)`; optional comparison against `mcscf.approx_hessian(mcscf.CASSCF(...)).kernel(mo)`; purpose: larger AO/virtual orbital space with a small exact-FCI CAS, useful for testing matrix-free AH/preconditioner behavior without making the CI solve the bottleneck.
```python
atom = """
C -0.65830719  0.61123287 -0.00800148
C  0.73685281  0.61123287 -0.00800148
C  1.43439081  1.81898387 -0.00800148
C  0.73673681  3.02749287 -0.00920048
C -0.65808819  3.02741487 -0.00967948
C -1.35568919  1.81920887 -0.00868348
H -1.20806619 -0.34108413 -0.00755148
H  1.28636081 -0.34128013 -0.00668648
H  2.53407081  1.81906387 -0.00736748
H  1.28693681  3.97963587 -0.00925948
H -1.20821019  3.97969587 -0.01063248
H -2.45529319  1.81939187 -0.00886348
"""
```

## Build
```bash
export SOURCE_REPO_ROOT="$(cd "$(git rev-parse --git-common-dir)/.." && pwd)"
export VENV="/anvil/scratch/x-tli22/fermilink_optimize/project_pyscf_casscf/venvs/fermilink-optimize/pyscf-casscf"
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
- Keep benchmark behavior deterministic with `OMP_NUM_THREADS=16`, `OPENBLAS_NUM_THREADS=1`, `MKL_NUM_THREADS=1`, and `NUMEXPR_NUM_THREADS=1`.
- Run the `## Build` commands from the PySCF repo root inside the campaign's active environment; the generated benchmark should use explicit build pre-commands derived from this section.
- Preserve the `train-` and `test-` ids directly and let FermiLink infer the split from the prefixes; do not add a manual `split` block when every case already uses those prefixes.
- Time only the selected CASSCF solver call after molecule setup, SCF setup, and active-orbital selection are complete.
- If the runner can expose them, record per-case `casscf_kernel_seconds`, `macro_cycles`, `micro_cycles`, `ah_seconds`, `h_op_calls`, `casci_seconds`, `ao2mo_seconds`, `keyframe_restarts`, `rejected_steps`, `converged`, `e_tot`, `e_states`, `norm_gorb`, `norm_gci`, and root-overlap diagnostics.
- Keep the initial optimize benchmark to fixed-geometry cases that finish reproducibly in repeated local runs; do not include whole geometry scans, multistage homotopy schedules, or very heavy Cr2-style stress cases in the initial autogen suite.
- If a listed case proves baseline-nonreproducible on the target machine, replace it before the campaign starts with a nearby fixed-geometry case in the same CASSCF regime rather than weakening tolerances or broadening the editable scope.
