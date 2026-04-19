# FermiLink Optimize Program

- package_id: lammps
- benchmark_id: lammps-tip4p-long-nve-pair-kspace

## Purpose
- Search the code implementation space for better performance while preserving benchmark correctness.
- Use one experiment at a time.
- Keep accepted changes simple, reproducible, and benchmark-backed.

## Workflow
1. Read the benchmark contract, skills, memory, and recent results.
2. Propose exactly one candidate change.
3. Edit only benchmark-approved source files.
4. Optionally run quick local checks, but do not run the authoritative benchmark command.
5. Stop after the candidate change is ready.

## Heuristics
- Prefer simpler accepted changes when gains are marginal.
- Avoid broad refactors that blur the causal source of improvement.
- When stuck, change one dominant hypothesis instead of many knobs at once.
- Never weaken tolerances or special-case benchmark inputs.
