#!/usr/bin/env python3
from __future__ import annotations

import argparse
import contextlib
import json
import os
import resource
import sys
import time
import traceback
from pathlib import Path
from typing import Any, Iterator

for _name in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_name, "1")
os.environ.setdefault("PYTHONHASHSEED", "0")

import yaml


CASE_META_KEYS = {
    "id",
    "weight",
    "description",
    "timeout_seconds",
    "command",
    "parameters",
}


def _load_benchmark(path: Path) -> dict[str, Any]:
    payload = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError("benchmark file must contain a YAML object")
    return payload


def _case_parameters(case: dict[str, Any]) -> dict[str, Any]:
    params: dict[str, Any] = {}
    nested = case.get("parameters")
    if isinstance(nested, dict):
        params.update(nested)
    for key, value in case.items():
        if key not in CASE_META_KEYS:
            params[key] = value
    return params


def _input_root() -> Path:
    raw = os.environ.get("FERMILINK_GOAL_INPUT_ROOT")
    if raw:
        return Path(raw).expanduser().resolve()
    return Path.cwd().resolve()


def _resolve_input_path(value: Any, input_root: Path) -> Path:
    path = Path(str(value)).expanduser()
    if path.is_absolute():
        return path.resolve()
    return (input_root / path).resolve()


def _case_atom(params: dict[str, Any], input_root: Path) -> Any:
    for key in ("atom_file", "geometry_file", "xyz_file"):
        if key in params and params[key]:
            return _resolve_input_path(params[key], input_root).read_text(encoding="utf-8")
    if "atom" not in params:
        raise ValueError("case requires atom or atom_file")
    return params["atom"]


def _to_jsonable_array(value: Any) -> Any:
    import numpy as np

    arr = np.asarray(value)
    if np.iscomplexobj(arr):
        if arr.size == 0 or float(np.max(np.abs(arr.imag))) < 1.0e-12:
            arr = arr.real
        else:
            arr = np.stack((arr.real, arr.imag), axis=-1)
    return arr.tolist()


def _real_float(value: Any) -> float:
    import numpy as np

    arr = np.asarray(value)
    if np.iscomplexobj(arr):
        value = arr.real
    return float(value)


def _energy_window_one_spin(mo_energy: Any, mo_occ: Any) -> list[float]:
    import numpy as np

    energies = np.asarray(mo_energy, dtype=float)
    occ = np.asarray(mo_occ)
    occupied = np.sort(energies[occ > 0])
    virtual = np.sort(energies[occ <= 0])[:10]
    return np.concatenate((occupied, virtual)).astype(float).tolist()


def _mo_energy_window(mo_energy: Any, mo_occ: Any) -> Any:
    import numpy as np

    energies = np.asarray(mo_energy, dtype=object)
    occ = np.asarray(mo_occ, dtype=object)
    if energies.ndim >= 2 and len(energies) == 2:
        return [
            _energy_window_one_spin(mo_energy[0], mo_occ[0]),
            _energy_window_one_spin(mo_energy[1], mo_occ[1]),
        ]
    return _energy_window_one_spin(mo_energy, mo_occ)


def _capture_norm_gorb(state: dict[str, float], grad: Any) -> None:
    import numpy as np
    from pyscf.scf import hf as scf_hf

    arr = np.asarray(grad)
    norm = float(np.linalg.norm(arr))
    if arr.size and not scf_hf.TIGHT_GRAD_CONV_TOL:
        norm /= float(np.sqrt(arr.size))
    state["norm_gorb"] = norm


def _wrap_timed_method(obj: Any, name: str, metric: str, timers: dict[str, float]) -> Any:
    original = getattr(obj, name)

    def wrapped(*args: Any, **kwargs: Any) -> Any:
        started = time.perf_counter()
        try:
            return original(*args, **kwargs)
        finally:
            timers[metric] += max(0.0, time.perf_counter() - started)

    setattr(obj, name, wrapped)
    return original


def _wrap_get_grad(obj: Any, state: dict[str, float]) -> Any:
    original = getattr(obj, "get_grad")

    def wrapped(*args: Any, **kwargs: Any) -> Any:
        grad = original(*args, **kwargs)
        _capture_norm_gorb(state, grad)
        return grad

    setattr(obj, "get_grad", wrapped)
    return original


@contextlib.contextmanager
def _time_diis_updates(timers: dict[str, float]) -> Iterator[None]:
    from pyscf.scf import diis as scf_diis

    patched: list[tuple[type[Any], Any]] = []
    for class_name in ("CDIIS", "ADIIS", "EDIIS"):
        cls = getattr(scf_diis, class_name, None)
        original = getattr(cls, "update", None)
        if cls is None or not callable(original):
            continue

        def make_wrapper(func: Any) -> Any:
            def wrapped(self: Any, *args: Any, **kwargs: Any) -> Any:
                started = time.perf_counter()
                try:
                    return func(self, *args, **kwargs)
                finally:
                    timers["diis_update_seconds"] += max(0.0, time.perf_counter() - started)

            return wrapped

        setattr(cls, "update", make_wrapper(original))
        patched.append((cls, original))
    try:
        yield
    finally:
        for cls, original in reversed(patched):
            setattr(cls, "update", original)


def _configure_mean_field(params: dict[str, Any], input_root: Path) -> Any:
    from pyscf import dft, gto, lib

    lib.num_threads(1)
    atom = _case_atom(params, input_root)
    mol = gto.M(
        atom=atom,
        unit=str(params.get("unit", "Angstrom")),
        basis=params["basis"],
        charge=int(params.get("charge", 0)),
        spin=int(params.get("spin", 0)),
        verbose=0,
    )
    method = str(params.get("method", "RKS")).upper()
    if method == "RKS":
        mf = dft.RKS(mol)
    elif method == "UKS":
        mf = dft.UKS(mol)
    else:
        raise ValueError(f"unsupported method {method!r}")

    mf.verbose = 0
    mf.xc = str(params["xc"])
    mf.init_guess = str(params.get("init_guess", "minao"))
    mf.grids.level = int(params.get("grids_level", params.get("grids.level", 3)))
    mf.conv_tol = float(params.get("conv_tol", mf.conv_tol))
    if params.get("conv_tol_grad") is not None:
        mf.conv_tol_grad = float(params["conv_tol_grad"])
    mf.max_cycle = int(params.get("max_cycle", mf.max_cycle))
    mf.diis_space = int(params.get("diis_space", mf.diis_space))

    for attr in ("diis_start_cycle", "diis_space_rollback"):
        if params.get(attr) is not None:
            setattr(mf, attr, int(params[attr]))
    for attr in ("diis_damp", "damp", "level_shift"):
        if params.get(attr) is not None:
            setattr(mf, attr, params[attr])
    if params.get("diis") is not None:
        mf.diis = bool(params["diis"])
    if params.get("direct_scf") is not None:
        mf.direct_scf = bool(params["direct_scf"])
    return mf


def _run_case(case: dict[str, Any], input_root: Path) -> dict[str, Any]:
    case_id = str(case.get("id") or "case")
    weight = float(case.get("weight", 1.0) or 1.0)
    result: dict[str, Any] = {
        "id": case_id,
        "weight": weight,
        "converged": False,
        "wall_seconds": 0.0,
        "total_seconds": 0.0,
        "scf_kernel_seconds": 0.0,
        "diis_update_seconds": 0.0,
        "get_fock_seconds": 0.0,
        "eig_seconds": 0.0,
        "error": "",
    }
    total_started = time.perf_counter()
    timers = {
        "diis_update_seconds": 0.0,
        "get_fock_seconds": 0.0,
        "eig_seconds": 0.0,
    }
    state: dict[str, float] = {}
    restored: list[tuple[Any, str, Any]] = []
    try:
        params = _case_parameters(case)
        mf = _configure_mean_field(params, input_root)
        for name, metric in (("get_fock", "get_fock_seconds"), ("eig", "eig_seconds")):
            restored.append((mf, name, _wrap_timed_method(mf, name, metric, timers)))
        restored.append((mf, "get_grad", _wrap_get_grad(mf, state)))

        kernel_started = time.perf_counter()
        try:
            with _time_diis_updates(timers):
                e_tot = mf.kernel()
        finally:
            result["scf_kernel_seconds"] = max(0.0, time.perf_counter() - kernel_started)
            for obj, name, original in reversed(restored):
                setattr(obj, name, original)

        density_matrix = mf.make_rdm1()
        result.update(
            {
                "converged": bool(mf.converged),
                "e_tot": _real_float(e_tot),
                "scf_cycles": int(getattr(mf, "cycles", 0) or 0),
                "norm_gorb": float(state.get("norm_gorb", 0.0)),
                "mo_energy_window": _mo_energy_window(mf.mo_energy, mf.mo_occ),
                "density_matrix": _to_jsonable_array(density_matrix),
                "wall_seconds": float(result["scf_kernel_seconds"]),
            }
        )
        if not result["converged"]:
            result["error"] = "SCF did not converge"
    except Exception as exc:
        result["converged"] = False
        result["error"] = "".join(traceback.format_exception_only(type(exc), exc)).strip()
        for obj, name, original in reversed(restored):
            try:
                setattr(obj, name, original)
            except Exception:
                pass
    finally:
        total_seconds = max(0.0, time.perf_counter() - total_started)
        result["total_seconds"] = total_seconds
        if not result.get("wall_seconds"):
            result["wall_seconds"] = total_seconds
        for key, value in timers.items():
            result[key] = float(value)
    return result


def _weighted_median(values_and_weights: list[tuple[float, float]]) -> float:
    clean = [(float(value), max(0.0, float(weight))) for value, weight in values_and_weights]
    clean = [(value, weight) for value, weight in clean if weight > 0.0]
    if not clean:
        return 0.0
    clean.sort(key=lambda item: item[0])
    half = sum(weight for _, weight in clean) * 0.5
    cumulative = 0.0
    for value, weight in clean:
        cumulative += weight
        if cumulative >= half:
            return float(value)
    return float(clean[-1][0])


def _primary_summary_metric(primary_metric: str, cases: list[dict[str, Any]]) -> float:
    if primary_metric.startswith("weighted_median_"):
        field = primary_metric[len("weighted_median_") :]
    else:
        field = "scf_kernel_seconds"
    return _weighted_median(
        [
            (float(case.get(field, case.get("wall_seconds", 0.0)) or 0.0), float(case.get("weight", 1.0) or 1.0))
            for case in cases
        ]
    )


def _peak_rss_mb() -> float:
    rss = float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    if sys.platform == "darwin":
        return rss / (1024.0 * 1024.0)
    return rss / 1024.0


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--benchmark", required=True)
    parser.add_argument("--emit-json", action="store_true")
    args = parser.parse_args()

    benchmark_path = Path(args.benchmark).expanduser().resolve()
    payload = _load_benchmark(benchmark_path)
    cases_raw = payload.get("cases")
    cases = [item for item in cases_raw if isinstance(item, dict)] if isinstance(cases_raw, list) else []
    primary_metric = str(
        (((payload.get("controller") or {}).get("objective") or {}).get("primary_metric"))
        or "weighted_median_scf_kernel_seconds"
    )

    input_root = _input_root()
    case_results = [_run_case(case, input_root) for case in cases]
    failures = sum(1 for item in case_results if not bool(item.get("converged")) or bool(item.get("error")))
    summary_metrics = {
        primary_metric: _primary_summary_metric(primary_metric, case_results),
        "peak_rss_mb": _peak_rss_mb(),
        "total_failures": int(failures),
        "weighted_median_diis_update_seconds": _primary_summary_metric(
            "weighted_median_diis_update_seconds", case_results
        ),
        "weighted_median_get_fock_seconds": _primary_summary_metric(
            "weighted_median_get_fock_seconds", case_results
        ),
        "weighted_median_eig_seconds": _primary_summary_metric("weighted_median_eig_seconds", case_results),
    }
    output = {
        "benchmark_id": str(payload.get("benchmark_id") or "pyscf-diis-scf-kernel"),
        "correctness_ok": failures == 0,
        "summary_metrics": summary_metrics,
        "cases": case_results,
    }
    print(json.dumps(output, sort_keys=True, separators=(",", ":")))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())