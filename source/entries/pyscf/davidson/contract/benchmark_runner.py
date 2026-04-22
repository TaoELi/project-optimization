#!/usr/bin/env python3
from __future__ import annotations

import argparse
import contextlib
import gc
import json
import os
import re
import resource
import sys
import time
import traceback
from pathlib import Path
from typing import Any, Callable, Iterator

np: Any = None
yaml: Any = None

THREAD_ENV_VARS = (
    "OMP_NUM_THREADS",
    "MKL_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "BLIS_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
)


def _load_runtime_modules() -> None:
    global np, yaml
    if np is None:
        import numpy as _np

        np = _np
    if yaml is None:
        import yaml as _yaml

        yaml = _yaml

BENZENE_ATOM = """
c   1.217739890298750 -0.703062453466927  0.000000000000000
h   2.172991468538160 -1.254577209307266  0.000000000000000
c   1.217739890298750  0.703062453466927  0.000000000000000
h   2.172991468538160  1.254577209307266  0.000000000000000
c   0.000000000000000  1.406124906933854  0.000000000000000
h   0.000000000000000  2.509154418614532  0.000000000000000
c  -1.217739890298750  0.703062453466927  0.000000000000000
h  -2.172991468538160  1.254577209307266  0.000000000000000
c  -1.217739890298750 -0.703062453466927  0.000000000000000
h  -2.172991468538160 -1.254577209307266  0.000000000000000
c   0.000000000000000 -1.406124906933854  0.000000000000000
h   0.000000000000000 -2.509154418614532  0.000000000000000
"""

ALLYL_ATOM = """
C    -1.1528    -0.1151    -0.4645
C     0.2300    -0.1171    -0.3508
C     0.9378     0.2246     0.7924
H     0.4206     0.5272     1.7055
H     2.0270     0.2021     0.8159
H    -1.6484    -0.3950    -1.3937
H    -1.7866     0.1687     0.3784
H     0.8086    -0.4120    -1.2337
"""


def _configure_thread_env() -> None:
    for name in THREAD_ENV_VARS:
        os.environ.setdefault(name, "1")


def _input_root() -> Path:
    raw = os.environ.get("FERMILINK_GOAL_INPUT_ROOT")
    if raw:
        return Path(raw).expanduser().resolve()
    return Path.cwd().resolve()


def _resolve_input_path(path_value: Any, input_root: Path) -> Path | None:
    if not path_value:
        return None
    path = Path(str(path_value)).expanduser()
    if path.is_absolute():
        return path.resolve()
    return (input_root / path).resolve()


def _extract_atom_block(source_path: Path) -> str | None:
    try:
        text = source_path.read_text(encoding="utf-8")
    except OSError:
        return None
    match = re.search(r"(?:mol\.)?atom\s*=\s*'''(.*?)'''", text, re.DOTALL)
    if match is None:
        return None
    atom = match.group(1).strip()
    return atom or None


def _case_atom(case: dict[str, Any], input_root: Path) -> tuple[str, str]:
    if case.get("atom"):
        return str(case["atom"]), "case.atom"

    source_path = _resolve_input_path(case.get("geometry_source"), input_root)
    if source_path is not None:
        atom = _extract_atom_block(source_path)
        if atom is not None:
            return atom, str(source_path)

    geometry_name = str(case.get("geometry_name") or "").lower()
    if geometry_name == "benzene":
        return BENZENE_ATOM.strip(), "embedded:benzene"
    if geometry_name == "allyl":
        return ALLYL_ATOM.strip(), "embedded:allyl"
    raise ValueError(f"unknown geometry for case {case.get('id')!r}")


def _to_float_list(value: Any) -> list[float]:
    if value is None:
        return []
    arr = np.asarray(value)
    if np.iscomplexobj(arr):
        arr = arr.real
    return [float(x) for x in arr.reshape(-1)]


def _to_bool_list(value: Any, expected: int) -> list[bool]:
    arr = np.asarray(value)
    if arr.shape == ():
        return [bool(arr.item())] * expected
    return [bool(x) for x in arr.reshape(-1)]


def _transition_dipole_norm(td_obj: Any) -> list[float]:
    dip = np.asarray(td_obj.transition_dipole(), dtype=float)
    if dip.size == 0:
        return []
    return [float(x) for x in np.linalg.norm(dip.reshape(dip.shape[0], -1), axis=1)]


def _weighted_median(values_and_weights: list[tuple[float, float]]) -> float:
    pairs = [(float(v), float(w)) for v, w in values_and_weights if float(w) > 0]
    if not pairs:
        return 0.0
    pairs.sort(key=lambda item: item[0])
    total_weight = sum(weight for _, weight in pairs)
    cutoff = total_weight * 0.5
    cumulative = 0.0
    for value, weight in pairs:
        cumulative += weight
        if cumulative >= cutoff:
            return value
    return pairs[-1][0]


def _rss_mb() -> float:
    usage = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    if sys.platform == "darwin":
        return float(usage) / (1024.0 * 1024.0)
    return float(usage) / 1024.0


class LRSolverStats:
    def __init__(self) -> None:
        self.calls: list[dict[str, Any]] = []

    @property
    def davidson_iterations(self) -> int:
        return int(sum(int(call.get("aop_calls", 0)) for call in self.calls))

    @property
    def matvec_applications(self) -> int:
        return int(sum(int(call.get("matvec_applications", 0)) for call in self.calls))


def _trial_count(zs: Any) -> int:
    arr = np.asarray(zs)
    if arr.ndim == 0:
        return 1
    return int(arr.shape[0]) if arr.ndim > 1 else 1


def _wrap_solver(name: str, solver: Callable[..., Any], stats: LRSolverStats) -> Callable[..., Any]:
    def wrapped(aop: Callable[[Any], Any], x0: Any, precond: Callable[..., Any], *args: Any, **kwargs: Any) -> Any:
        call_stats: dict[str, Any] = {
            "solver": name,
            "aop_calls": 0,
            "matvec_applications": 0,
            "elapsed_seconds": 0.0,
        }

        def counted_aop(zs: Any) -> Any:
            count = _trial_count(zs)
            call_stats["aop_calls"] += 1
            call_stats["matvec_applications"] += count
            return aop(zs)

        started = time.perf_counter()
        try:
            return solver(counted_aop, x0, precond, *args, **kwargs)
        finally:
            call_stats["elapsed_seconds"] = time.perf_counter() - started
            stats.calls.append(call_stats)

    return wrapped


@contextlib.contextmanager
def _patched_lr_solvers(stats: LRSolverStats) -> Iterator[None]:
    from pyscf.tdscf import _lr_eig, rhf, rks, uhf, uks

    targets = [
        (_lr_eig, "eigh", "eigh"),
        (_lr_eig, "eig", "eig"),
        (_lr_eig, "real_eig", "real_eig"),
        (rhf, "lr_eigh", "eigh"),
        (rhf, "lr_eig", "eig"),
        (rhf, "real_eig", "real_eig"),
        (uhf, "lr_eigh", "eigh"),
        (uhf, "lr_eig", "eig"),
        (uhf, "real_eig", "real_eig"),
        (rks, "lr_eigh", "eigh"),
        (uks, "lr_eigh", "eigh"),
    ]
    originals: list[tuple[Any, str, Any]] = []
    try:
        for module, attr, solver_name in targets:
            if hasattr(module, attr):
                original = getattr(module, attr)
                originals.append((module, attr, original))
                setattr(module, attr, _wrap_solver(solver_name, original, stats))
        yield
    finally:
        for module, attr, original in reversed(originals):
            setattr(module, attr, original)


def _build_scf(case: dict[str, Any], atom: str) -> Any:
    from pyscf import dft, gto, lib, scf

    lib.num_threads(1)
    mol = gto.M(
        atom=atom,
        basis=str(case["basis"]),
        charge=int(case.get("charge", 0)),
        spin=int(case.get("spin", 0)),
        symmetry=bool(case.get("symmetry", False)),
        verbose=0,
    )
    scf_method = str(case["scf_method"])
    if scf_method in {"RKS", "UKS"}:
        mf = getattr(dft, scf_method)(mol)
        mf.xc = str(case["xc"])
        mf.grids.prune = None
    else:
        mf = getattr(scf, scf_method)(mol)
    mf.max_memory = float(case.get("max_memory", 4000))
    mf.conv_tol = float(case.get("scf_conv_tol", 1.0e-10))
    mf.verbose = 0
    mf.chkfile = None
    return mf


def _build_td(case: dict[str, Any], mf: Any) -> Any:
    td_method = str(case["td_method"])
    frozen = case.get("frozen")
    td = getattr(mf, td_method)(frozen)
    td.nstates = int(case["nstates"])
    td.conv_tol = float(case.get("td_conv_tol", case.get("conv_tol", td.conv_tol)))
    td.lindep = float(case.get("lindep", td.lindep))
    td.max_cycle = int(case.get("max_cycle", td.max_cycle))
    td.max_memory = float(case.get("max_memory", td.max_memory))
    td.positive_eig_threshold = float(
        case.get("positive_eig_threshold", td.positive_eig_threshold)
    )
    td.deg_eia_thresh = float(case.get("deg_eia_thresh", td.deg_eia_thresh))
    if case.get("wfnsym") is not None:
        td.wfnsym = case["wfnsym"]
    if case.get("singlet") is not None and hasattr(td, "singlet"):
        td.singlet = bool(case["singlet"])
    td.verbose = 0
    td.chkfile = None
    return td


def _run_case(case: dict[str, Any], input_root: Path) -> dict[str, Any]:
    from pyscf import dft

    dft.radi.ATOM_SPECIFIC_TREUTLER_GRIDS = False

    case_id = str(case.get("id") or "case")
    weight = float(case.get("weight", 1.0))
    started = time.perf_counter()
    stats = LRSolverStats()
    try:
        atom, geometry_resolved_from = _case_atom(case, input_root)
        mf = _build_scf(case, atom)

        scf_started = time.perf_counter()
        scf_energy = mf.kernel()
        scf_wall = time.perf_counter() - scf_started
        if not bool(getattr(mf, "converged", False)):
            raise RuntimeError("ground-state SCF did not converge")

        td = _build_td(case, mf)
        td_started = time.perf_counter()
        with _patched_lr_solvers(stats):
            e, _xy = td.kernel()
        td_wall = time.perf_counter() - td_started

        expected_roots = int(case["nstates"])
        energies = _to_float_list(e)
        conv = _to_bool_list(getattr(td, "converged", []), len(energies))
        roots_ok = len(energies) >= expected_roots and len(conv) >= expected_roots
        roots_ok = roots_ok and all(conv[:expected_roots])

        oscillator_strength: list[float] = []
        transition_dipole_norm: list[float] = []
        if bool(case.get("oscillator_strength", False)):
            oscillator_strength = _to_float_list(td.oscillator_strength())
            transition_dipole_norm = _transition_dipole_norm(td)
        else:
            oscillator_strength = [0.0]
            transition_dipole_norm = [0.0]

        total_seconds = time.perf_counter() - started
        result = {
            "id": case_id,
            "converged": bool(roots_ok),
            "wall_seconds": float(td_wall),
            "total_seconds": float(total_seconds),
            "td_kernel_wall_time_s": float(td_wall),
            "scf_wall_time_s": float(scf_wall),
            "scf_energy": float(scf_energy),
            "e": energies,
            "converged_roots": conv,
            "root_count": int(len(energies)),
            "oscillator_strength": oscillator_strength,
            "transition_dipole_norm": transition_dipole_norm,
            "davidson_iterations": stats.davidson_iterations,
            "matvec_applications": stats.matvec_applications,
            "lr_solver_calls": stats.calls,
            "weight": weight,
            "geometry_resolved_from": geometry_resolved_from,
            "error": "" if roots_ok else "not all requested TD roots converged",
        }
        return result
    except BaseException as exc:
        total_seconds = time.perf_counter() - started
        error = "".join(traceback.format_exception_only(type(exc), exc)).strip()
        return {
            "id": case_id,
            "converged": False,
            "wall_seconds": 0.0,
            "total_seconds": float(total_seconds),
            "td_kernel_wall_time_s": 0.0,
            "e": [],
            "converged_roots": [],
            "root_count": 0,
            "oscillator_strength": [],
            "transition_dipole_norm": [],
            "davidson_iterations": stats.davidson_iterations,
            "matvec_applications": stats.matvec_applications,
            "lr_solver_calls": stats.calls,
            "weight": weight,
            "error": error,
        }
    finally:
        gc.collect()


def _json_default(value: Any) -> Any:
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, complex):
        return float(value.real)
    return str(value)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--benchmark", required=True)
    parser.add_argument("--emit-json", action="store_true")
    args = parser.parse_args()

    _configure_thread_env()
    _load_runtime_modules()
    benchmark_path = Path(args.benchmark).expanduser().resolve()
    payload = yaml.safe_load(benchmark_path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise SystemExit("benchmark file must be a YAML object")
    cases_raw = payload.get("cases")
    if not isinstance(cases_raw, list) or not cases_raw:
        raise SystemExit("benchmark file requires a non-empty cases list")

    input_root = _input_root()
    case_results = [_run_case(case, input_root) for case in cases_raw if isinstance(case, dict)]
    primary_metric = str(
        payload.get("controller", {})
        .get("objective", {})
        .get("primary_metric", "weighted_median_td_kernel_wall_seconds")
    )
    weighted_wall = _weighted_median(
        [
            (float(result.get("td_kernel_wall_time_s") or 0.0), float(result.get("weight") or 1.0))
            for result in case_results
        ]
    )
    failures = sum(1 for result in case_results if not bool(result.get("converged")))
    total_iterations = sum(int(result.get("davidson_iterations") or 0) for result in case_results)
    total_matvecs = sum(int(result.get("matvec_applications") or 0) for result in case_results)
    output = {
        "benchmark_id": str(payload.get("benchmark_id") or "pyscf-tddft-davidson"),
        "correctness_ok": failures == 0,
        "summary_metrics": {
            primary_metric: float(weighted_wall),
            "peak_rss_mb": _rss_mb(),
            "total_failures": int(failures),
            "total_davidson_iterations": int(total_iterations),
            "total_matvec_applications": int(total_matvecs),
        },
        "cases": case_results,
    }
    print(json.dumps(output, sort_keys=True, default=_json_default))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())