# FermiLink Optimize Worker Memory

- package_id: lammps
- benchmark_id: lammps-tip4p-long-nve-pair-kspace
- benchmark_path: .fermilink-optimize/benchmark.worker.yaml
- program_path: .fermilink-optimize/program.md
- controller_memory_path: .fermilink-optimize/memory.md
- results_path: .fermilink-optimize/results.tsv
- worker_iteration: 39
- reset_at_utc: 2026-04-13T06:37:32.234661Z

This file is reset at the start of each outer optimize iteration and archived under `.fermilink-optimize/runs/iter_XXXX/worker_memory.md`.

## Short-Term Memory (Operational)
### Current objective
- Prepare exactly one candidate that is ready for authoritative benchmark evaluation.
### Plan
- [x] Read benchmark/program/controller memory/results/skills.
- [x] Inspect the current implementation and pick one optimization hypothesis.
- [x] Apply a single focused experiment: split TIP4P cache hits into `xM`-only vs hydrogen-index reads, and memoize the pair-side oxygen reach scalar while preserving exact force/tally geometry.
- [x] Rebuild and run a short local `timer full` smoke check on the 16-rank TIP4P long input.
- [x] Update this memory with the local outcome and finish only when the candidate is benchmark-ready.
### Progress log
- Worker iteration initialized.
- Read the benchmark contract, optimize program, controller memory, worker memory, results ledger, and local LAMMPS skills at the start of the turn.
- Reviewed the current incumbent source plus recent controller notes; the strongest post-incumbent direction remains semantics-preserving hot-path access cleanup, while eager all-oxygen preparation, broader PPPM cache expansion, and oxygen-geometry rewrites are explicitly ruled out.
- Ran a short local 16-rank `timer full` baseline on `in.tip4p_nve_long` with `run 400`; it completed in 6.158 s total with Pair at 65.18% and Kspace at 21.56%, confirming that the next candidate should still lean pair-heavy while trimming paired PPPM cache-hit overhead where it is nearly free.
- Implemented the single candidate in the approved TIP4P pair/PPPM files: pair cache hits now separate `xM` coordinate reads from hydrogen-index reads, pair prechecks reuse a memoized per-site `cut_coul + |OM|` scalar without changing the exact surviving M-site separation math, and PPPM `particle_map()` now uses an `xM`-only cached accessor.
- Rebuilt `build/lmp` successfully after the source changes.
- Re-ran the same short local 16-rank `timer full` check on `in.tip4p_nve_long` with `run 400`; the thermo checkpoints stayed identical at steps 0/100/200/300/400, total loop time improved slightly from 6.15814 s to 6.15144 s, Pair improved from 4.0139 s to 4.0027 s, and Kspace stayed effectively flat at 1.3274 s -> 1.3280 s.

## Tactical Notes
### Job tracking
- none yet
### Candidate summary
- Ready for controller benchmark: add `xM`-only TIP4P cache accessors in pair/PPPM, defer pair hydrogen-index reads until exact oxygen survivors need force distribution, and cache the pair-side oxygen reach scalar with the existing per-step `xM` update so prechecks stop recomputing that square root.
### Debug notes
- Local profiling/tooling note: `perf`, `valgrind`, and `gprofng` are unavailable in this environment; local direction-finding used a shortened MPI run with `timer full`.
