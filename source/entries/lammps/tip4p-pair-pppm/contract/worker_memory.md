# FermiLink Optimize Worker Memory

- package_id: lammps
- benchmark_id: lammps-tip4p-water-nve-v1
- benchmark_path: .fermilink-optimize/benchmark.worker.yaml
- program_path: .fermilink-optimize/program.md
- controller_memory_path: .fermilink-optimize/memory.md
- results_path: .fermilink-optimize/results.tsv
- worker_iteration: 18
- reset_at_utc: 2026-04-09T14:47:15.677888Z

This file is reset at the start of each outer optimize iteration and archived under `.fermilink-optimize/runs/iter_XXXX/worker_memory.md`.

## Short-Term Memory (Operational)
### Current objective
- Prepare exactly one candidate that is ready for authoritative benchmark evaluation.
### Plan
- Read benchmark/program/controller memory/results/skills.
- Inspect the current implementation and pick one optimization hypothesis.
- Apply focused edits only within benchmark-approved paths.
- Run quick local checks or launch/poll long-running worker jobs if needed.
- Update this memory with progress and finish only when the candidate is benchmark-ready.
### Progress log
- Worker iteration initialized.

## Tactical Notes
### Job tracking
- none yet
### Candidate summary
- pending
### Debug notes
- none yet
