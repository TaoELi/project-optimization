# FermiLink Optimize Worker Memory

- package_id: pyscf
- benchmark_id: pyscf-diis-scf-kernel-v1
- benchmark_path: .fermilink-optimize/benchmark.worker.yaml
- program_path: .fermilink-optimize/program.md
- controller_memory_path: .fermilink-optimize/memory.md
- results_path: .fermilink-optimize/results.tsv
- worker_iteration: 27
- reset_at_utc: 2026-04-27T07:35:22.383291Z

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
- Read required benchmark contract, optimize program, controller memory, worker memory, results ledger, and PySCF skills index/API notes.
- Current incumbent is 8d8a60d41b1b at weighted median 0.812246244983; a promotable candidate needs about 0.796001 seconds or lower.
- Recent rejected experiments show UKS-only spin-summed Coulomb and checkpoint-I/O recombinations are sub-threshold; selected a shared RKS/UKS DFT NumInt cached-path experiment instead.
- Edited `pyscf/dft/rks.py` only: added a helper mirroring stock NumInt block-size selection and allowed the accepted AO cache to retain the original AO tensor when a stock block loop will produce exactly one block; multi-block/custom-loop paths still copy AO arrays as before.
- Quick checks passed: `python -m py_compile pyscf/dft/rks.py pyscf/dft/uks.py`; direct H4/N2 RKS NumInt comparisons against uncached direct calls gave zero nelec/exc/vmat deltas; visible H4 and N2 RKS kernels converged with unchanged 9 and 6 cycles; a small UKS O2 smoke converged in 7 cycles.
- Local A/B timing by wrapping `block_loop` to force the old copy path was noisy but directionally neutral-to-faster for N2 and unchanged for H4; no authoritative benchmark command was run.

## Tactical Notes
### Job tracking
- none yet
### Candidate summary
- Experiment candidate: avoid the initial AO tensor copy when the accepted cached NumInt path is populated from a single stock `NumInt.block_loop` block, preserving the existing copy behavior for multi-block/custom-loop cases.
### Debug notes
- The current accepted AO cache wraps `NumInt.block_loop` in `pyscf/dft/rks.py` and already caches AO blocks; remaining `nr_rks`/`nr_uks` paths still rebuild `pair_mask`, `ao_loc`, `nbins`, and rho-evaluator screening metadata each call.
- Direct probes on the visible H4/N2 RKS grids showed each GGA AO cache stores one block around 22-23 MB; a dense RKS contraction path was numerically identical but too small by itself, so the source edit targets first-call AO cache population overhead instead.
