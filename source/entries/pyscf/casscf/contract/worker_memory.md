# FermiLink Optimize Worker Memory

- package_id: pyscf
- benchmark_id: pyscf-casscf-ciah-mc1step-autogen-v1
- benchmark_path: .fermilink-optimize/benchmark.worker.yaml
- program_path: .fermilink-optimize/program.md
- controller_memory_path: .fermilink-optimize/memory.md
- results_path: .fermilink-optimize/results.tsv
- worker_iteration: 21
- reset_at_utc: 2026-04-27T07:56:32.527366Z

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
- 2026-04-27T07:57:53Z: Read the benchmark contract, optimize program, controller memory, worker memory, results ledger, and local PySCF skills index/advanced-topic notes. Current incumbent is iteration 13 (`911063b081d2`) at `weighted_median_casscf_kernel_seconds=39.745700547`.
- 2026-04-27T07:57:53Z: Current `pyscf/mcscf/newton_casscf.py` already includes the accepted same-run materialized `ppaa`/`papa` path, vectorized H_co contractions, and cached low-rank DF A-side transforms. The symmetry-enabled benchmark delegates through `newton_casscf_symm.CASSCF` to `newton_casscf.kernel`.
- 2026-04-27T08:11:25Z: Custom non-authoritative cc-pVTZ run of the incumbent reproduced 3 macro cycles with per-macro AH counts 17/17/6 H-op and 3/5/2 keyframes, with about 23.9s in AO2MO and 15.4s in AH. This confirms both exact AO2MO and repeated AH DF block unpacking remain relevant.
- 2026-04-27T08:11:25Z: Candidate experiment applied: cache same-run unpacked DF AO blocks on the AH DF object when memory allows, and cache exact per-ERIS hcore MO transforms for reuse in `gen_g_hop` and fast CASCI setup. Edited only `pyscf/mcscf/newton_casscf.py`, `pyscf/mcscf/mc1step.py`, and this worker memory.
- 2026-04-27T08:42:35Z: Final candidate additionally uses the exact ERIS `_h1e_mo` plus materialized `ppaa`/`papa` tensors to build the canonicalization Fock directly in the current MO basis when `cas_natorb` is disabled, avoiding the AO Fock back-transform and immediate MO projection while preserving the normal fallback.
- 2026-04-27T08:42:35Z: Syntax check passed with `python -m py_compile pyscf/mcscf/newton_casscf.py pyscf/mcscf/mc1step.py pyscf/mcscf/casci.py`.
- 2026-04-27T08:42:35Z: Quick non-authoritative 6-31G smoke converged with `e_tot=-230.69789482604548`, `e_cas=-6.470336536383968`, `norm_gorb=4.700029891618346e-06`, and 0.703s kernel time.
- 2026-04-27T08:42:35Z: Final non-authoritative cc-pVTZ smoke converged in 3 macro cycles, 10 keyframes, and 40 Hx steps with `e_tot=-230.84849342127535`, `e_cas=-6.472898855638846`, final logged `|grad|=4.3e-05`, and 38.928s kernel time. This preserves the incumbent solver counts and exact energy fields in the custom check while landing near the 2% threshold versus the 39.745700547s incumbent.

## Tactical Notes
### Job tracking
- none yet
### Candidate summary
- Same-run invariant-transform caching and exact MO-basis canonicalization for Newton-CASSCF: guarded unpacked DF AO block cache for low-rank AH response contractions, exact hcore MO-transform reuse on ERIS objects, fast CASCI setup from `_h1e_mo`, and non-natorb canonicalization from the exact materialized `ppaa`/`papa` MO-basis Fock.
### Debug notes
- `.git` is hidden in this worktree as `.git.fermilink-hidden`; source inspection is proceeding without relying on git metadata.
