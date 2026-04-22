# FermiLink Optimize Worker Memory

- package_id: pyscf
- benchmark_id: pyscf-tddft-davidson-lr-eig-20260421
- benchmark_path: .fermilink-optimize/benchmark.worker.yaml
- program_path: .fermilink-optimize/program.md
- controller_memory_path: .fermilink-optimize/memory.md
- results_path: .fermilink-optimize/results.tsv
- worker_iteration: 17
- reset_at_utc: 2026-04-22T00:20:02.223813Z

This file is reset at the start of each outer optimize iteration and archived under `.fermilink-optimize/runs/iter_XXXX/worker_memory.md`.

## Short-Term Memory (Operational)
### Current objective
- Prepare exactly one candidate that is ready for authoritative benchmark evaluation.
### Plan
- [x] Read benchmark/program/controller memory/results/skills.
- [x] Inspect the current implementation and pick one optimization hypothesis.
- [x] Apply focused edits only within benchmark-approved paths.
- [x] Run quick local checks or launch/poll long-running worker jobs if needed.
- [x] Update this memory with progress and finish only when the candidate is benchmark-ready.
### Progress log
- Worker iteration initialized.
- Read benchmark contract, optimize program, controller memory, worker memory, results ledger, and PySCF skill files at start of loop.
- Confirmed worktree HEAD is incumbent `64da7449e`; no tracked source changes are present before this experiment.
- Recent rejected candidates show overhead-only TD response and preconditioner batching changes did not reduce total Davidson matvecs enough; this loop will target solver work generation instead.
- Local B3LYP `real_eig` probe on incumbent used response-vector batches `[10, 10, 10, 10, 8, 2]` (50 vectors); nearby fixed shifts, small nonzero initial `Y`, small extra search roots, hybrid Casida, Hermitian response splitting, and extra late correction directions did not reduce this count without changing energies or increasing work.
- Applied final source candidate in `pyscf/tdscf/rhf.py`: restricted full TDDFT now passes occupied-only MO data into response-cache setup and avoids repeated conjugation, diagonal `einsum` scaling, and `hstack` allocation in the hot matvec closure.
- Quick checks passed: `python -m py_compile pyscf/tdscf/rhf.py pyscf/tdscf/_lr_eig.py`; direct B3LYP `real_eig` smoke converged all 10 roots with batches `[10, 10, 10, 10, 8, 2]` (50 vectors) and excitation energies matching the incumbent probe.

## Tactical Notes
### Job tracking
- none yet
### Candidate summary
- Candidate hypothesis: reduce B3LYP restricted full TDDFT kernel overhead without perturbing Davidson convergence by using occupied-only response-cache inputs and cheaper local array operations in `rhf._gen_tdhf_operation`.
### Debug notes
- The final candidate intentionally leaves `_lr_eig.real_eig` convergence behavior unchanged after local probes showed count-oriented tuning was unfavorable.
