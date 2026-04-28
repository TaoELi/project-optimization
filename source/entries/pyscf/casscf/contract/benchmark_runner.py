#!/usr/bin/env python3
from __future__ import annotations

import argparse
import contextlib
import json
import math
import os
import resource
import sys
import time
import traceback
from pathlib import Path
from typing import Any, Iterator

import yaml


PINNED_PYTHON = (
    "/anvil/scratch/x-tli22/fermilink_optimize/project_pyscf_casscf/"
    "venvs/fermilink-optimize/pyscf-casscf/bin/python"
)


def _jsonable(value: Any) -> Any:
    try:
        import numpy
    except Exception:
        numpy = None

    if value is None or isinstance(value, (str, bool, int)):
        return value
    if isinstance(value, float):
        return value if math.isfinite(value) else None
    if isinstance(value, complex):
        if abs(value.imag) <= 1.0e-12:
            return float(value.real)
        return {"real": float(value.real), "imag": float(value.imag)}
    if numpy is not None:
        if isinstance(value, numpy.generic):
            return _jsonable(value.item())
        if isinstance(value, numpy.ndarray):
            return _jsonable(value.tolist())
    if isinstance(value, dict):
        return {str(k): _jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(v) for v in value]
    return str(value)


def _as_float(value: Any) -> float | None:
    if value is None:
        return None
    try:
        out = float(value)
    except (TypeError, ValueError):
        return None
    return out if math.isfinite(out) else None


def _as_int(value: Any) -> int | None:
    if value is None:
        return None
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def _weighted_median(values_and_weights: list[tuple[float, float]]) -> float:
    pairs = sorted(
        (float(value), float(weight))
        for value, weight in values_and_weights
        if math.isfinite(float(value)) and math.isfinite(float(weight)) and float(weight) > 0.0
    )
    if not pairs:
        return 0.0
    halfway = sum(weight for _, weight in pairs) * 0.5
    running = 0.0
    for value, weight in pairs:
        running += weight
        if running >= halfway:
            return value
    return pairs[-1][0]


def _peak_rss_mb() -> float:
    rss = float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    if sys.platform == "darwin":
        return rss / (1024.0 * 1024.0)
    return rss / 1024.0


def _find_project_root(start: Path) -> Path:
    current = start.resolve()
    if current.is_file():
        current = current.parent
    for path in (current, *current.parents):
        if (path / "pyscf").is_dir() and (path / "setup.py").is_file():
            return path
    raise RuntimeError(f"could not locate PySCF project root from {start}")


def _configure_threads() -> None:
    for key, value in (
        ("OMP_NUM_THREADS", "16"),
        ("OPENBLAS_NUM_THREADS", "1"),
        ("MKL_NUM_THREADS", "1"),
        ("NUMEXPR_NUM_THREADS", "1"),
        ("PYTHONHASHSEED", "0"),
    ):
        os.environ.setdefault(key, value)
    try:
        from pyscf import lib
    except Exception:
        return
    try:
        threads = int(os.environ.get("OMP_NUM_THREADS") or "16")
    except ValueError:
        threads = 16
    if threads > 0:
        lib.num_threads(threads)


def _resolve_case_paths(value: Any, input_root: Path | None, key: str = "") -> Any:
    if input_root is None:
        return value
    if isinstance(value, dict):
        return {k: _resolve_case_paths(v, input_root, str(k)) for k, v in value.items()}
    if isinstance(value, list):
        return [_resolve_case_paths(v, input_root, key) for v in value]
    if isinstance(value, str) and (key.endswith("_path") or key.endswith("_file")):
        path = Path(value)
        if not path.is_absolute():
            return str((input_root / path).resolve())
    return value


def _make_mol(case: dict[str, Any]) -> Any:
    from pyscf import gto

    atom = case.get("atom")
    if not atom:
        atom_file = case.get("atom_file") or case.get("geometry_file")
        if atom_file:
            atom = Path(str(atom_file)).read_text(encoding="utf-8")
    if not atom:
        raise ValueError("case must define atom or atom_file")
    mol = gto.Mole()
    mol.atom = atom
    mol.basis = case.get("basis")
    mol.charge = int(case.get("charge", 0))
    mol.spin = int(case.get("spin", 0))
    mol.symmetry = bool(case.get("symmetry", False))
    mol.unit = str(case.get("unit", "Angstrom"))
    mol.verbose = int(case.get("mol_verbose", 0))
    if case.get("max_memory") is not None:
        mol.max_memory = int(case["max_memory"])
    mol.build()
    return mol


def _make_mf(mol: Any, case: dict[str, Any]) -> Any:
    from pyscf import scf

    mf_type = str(case.get("mf_type") or "RHF").upper()
    if mf_type == "RHF":
        mf = scf.RHF(mol)
    elif mf_type == "ROHF":
        mf = scf.ROHF(mol)
    elif mf_type == "UHF":
        mf = scf.UHF(mol)
    else:
        raise ValueError(f"unsupported mf_type {mf_type!r}")
    mf.verbose = int(case.get("mf_verbose", 0))
    if case.get("mf_conv_tol") is not None:
        mf.conv_tol = float(case["mf_conv_tol"])
    if case.get("mf_max_cycle") is not None:
        mf.max_cycle = int(case["mf_max_cycle"])
    mf.kernel()
    if not bool(getattr(mf, "converged", False)):
        raise RuntimeError("mean-field reference did not converge")
    return mf


def _make_state_average_mix(base_mc: Any, mol: Any, case: dict[str, Any]) -> Any:
    from pyscf import fci
    from pyscf.mcscf import addons

    weights = case.get("state_weights") or [0.25, 0.25, 0.25, 0.25]
    mix = case.get("state_average_mix") or {}
    mode = str(mix.get("mode") or "pointgroup_sa4")
    if mode != "pointgroup_sa4":
        raise ValueError(f"unsupported state_average_mix mode {mode!r}")

    irreps = list(mix.get("irreps") or ["A1", "B1"])
    nroots = list(mix.get("nroots_per_solver") or [2, 2])
    singlet = bool(mix.get("singlet", False))
    solvers = []
    for wfnsym, nroot in zip(irreps, nroots):
        solver = fci.solver(mol, symm=bool(case.get("symmetry", False)), singlet=singlet)
        solver.nroots = int(nroot)
        solver.wfnsym = wfnsym
        solvers.append(solver)
    return addons.state_average_mix(base_mc, solvers, weights)


def _select_solver(base_mc: Any, mol: Any, case: dict[str, Any]) -> Any:
    from pyscf import mcscf

    family = str(case.get("solver_family") or "newton")
    weights = case.get("state_weights")
    if family == "newton_state_average_mix":
        return _make_state_average_mix(base_mc, mol, case).newton()
    if family == "state_average_mix":
        return _make_state_average_mix(base_mc, mol, case)
    if weights:
        averaged = base_mc.state_average(weights, case.get("wfnsym"))
        return averaged.newton() if family == "newton_state_average" else averaged
    if case.get("target_root") is not None:
        base_mc = base_mc.state_specific_(state=int(case["target_root"]))
    if family == "newton":
        return base_mc.newton()
    if family in ("mc1step", "mc2step"):
        return base_mc
    if family == "approx_hessian":
        return mcscf.approx_hessian(base_mc)
    raise ValueError(f"unsupported solver_family {family!r}")


def _apply_solver_options(solver: Any, case: dict[str, Any]) -> None:
    option_keys = (
        "conv_tol",
        "conv_tol_grad",
        "max_cycle_macro",
        "max_cycle_micro",
        "max_stepsize",
        "ah_level_shift",
        "ah_conv_tol",
        "ah_lindep",
        "ah_start_tol",
        "ah_start_cycle",
        "ah_max_cycle",
        "internal_rotation",
    )
    for key in option_keys:
        if key in case and case[key] is not None and hasattr(solver, key):
            setattr(solver, key, case[key])

    if case.get("wfnsym") is not None and hasattr(solver, "fcisolver"):
        solver.fcisolver.wfnsym = case["wfnsym"]

    if case.get("fcisolver_conv_tol") is not None and hasattr(solver, "fcisolver"):
        conv_tol = float(case["fcisolver_conv_tol"])
        fcisolver = solver.fcisolver
        if hasattr(fcisolver, "fcisolvers"):
            for subsolver in fcisolver.fcisolvers:
                subsolver.conv_tol = conv_tol
        else:
            fcisolver.conv_tol = conv_tol

    solver.verbose = int(case.get("casscf_verbose", 0))


@contextlib.contextmanager
def _time_instance_method(obj: Any, name: str, metric: str, stats: dict[str, Any]) -> Iterator[None]:
    if not hasattr(obj, name):
        yield
        return
    original = getattr(obj, name)

    def wrapped(*args: Any, **kwargs: Any) -> Any:
        started = time.perf_counter()
        try:
            return original(*args, **kwargs)
        finally:
            stats[metric] = float(stats.get(metric, 0.0)) + max(0.0, time.perf_counter() - started)

    setattr(obj, name, wrapped)
    try:
        yield
    finally:
        setattr(obj, name, original)


@contextlib.contextmanager
def _time_newton_ah(stats: dict[str, Any]) -> Iterator[None]:
    try:
        from pyscf.mcscf import newton_casscf
    except Exception:
        yield
        return
    original = newton_casscf.update_orb_ci

    def wrapped(*args: Any, **kwargs: Any) -> Any:
        started = time.perf_counter()
        try:
            return original(*args, **kwargs)
        finally:
            stats["ah_seconds"] = float(stats.get("ah_seconds", 0.0)) + max(
                0.0, time.perf_counter() - started
            )

    newton_casscf.update_orb_ci = wrapped
    try:
        yield
    finally:
        newton_casscf.update_orb_ci = original


@contextlib.contextmanager
def _time_rotate_orb_cc(solver: Any, stats: dict[str, Any]) -> Iterator[None]:
    if not hasattr(solver, "rotate_orb_cc"):
        yield
        return
    original = solver.rotate_orb_cc

    class TimedGenerator:
        def __init__(self, generator: Any):
            self._generator = generator

        def __iter__(self) -> "TimedGenerator":
            return self

        def __next__(self) -> Any:
            started = time.perf_counter()
            try:
                return next(self._generator)
            finally:
                stats["ah_seconds"] = float(stats.get("ah_seconds", 0.0)) + max(
                    0.0, time.perf_counter() - started
                )

        def close(self) -> Any:
            close = getattr(self._generator, "close", None)
            if callable(close):
                return close()
            return None

    def wrapped(*args: Any, **kwargs: Any) -> Any:
        return TimedGenerator(original(*args, **kwargs))

    setattr(solver, "rotate_orb_cc", wrapped)
    try:
        yield
    finally:
        setattr(solver, "rotate_orb_cc", original)


def _make_callback(stats: dict[str, Any]) -> Any:
    micro_by_macro: dict[int, int] = {}

    def callback(envs: dict[str, Any]) -> None:
        imacro = _as_int(envs.get("imacro"))
        if imacro is not None:
            stats["macro_cycles"] = max(int(stats.get("macro_cycles", 0)), imacro)

        totmicro = _as_int(envs.get("totmicro"))
        if totmicro is not None:
            stats["micro_cycles"] = max(int(stats.get("micro_cycles", 0)), totmicro)
        elif imacro is not None and envs.get("stat") is not None:
            imic = _as_int(getattr(envs["stat"], "imic", None))
            if imic is not None:
                micro_by_macro[imacro] = imic
                stats["micro_cycles"] = int(sum(micro_by_macro.values()))
        else:
            imicro = _as_int(envs.get("imicro"))
            if imicro is not None:
                stats["micro_cycles"] = max(int(stats.get("micro_cycles", 0)), imicro)

        for key in ("tot_hop", "totinner", "njk"):
            value = _as_int(envs.get(key))
            if value is not None:
                stats["h_op_calls"] = max(int(stats.get("h_op_calls", 0)), value)

        tot_kf = _as_int(envs.get("tot_kf"))
        if tot_kf is not None:
            stats["keyframe_restarts"] = max(int(stats.get("keyframe_restarts", 0)), tot_kf)

        for out_key, candidates in (
            ("norm_gorb", ("norm_gorb", "norm_gorb0", "norm_gall")),
            ("norm_gci", ("norm_gci",)),
        ):
            for candidate in candidates:
                value = _as_float(envs.get(candidate))
                if value is not None:
                    stats[out_key] = value
                    break

    return callback


def _ci_vectors(ci: Any) -> list[Any]:
    try:
        import numpy
    except Exception:
        return []
    if ci is None:
        return []
    if isinstance(ci, numpy.ndarray):
        return [ci]
    if isinstance(ci, (list, tuple)):
        out = []
        for item in ci:
            if isinstance(item, numpy.ndarray):
                out.append(item)
        return out
    return []


def _root_overlap(ci: Any) -> list[list[float]]:
    try:
        import numpy
    except Exception:
        return []
    vectors = []
    for vector in _ci_vectors(ci):
        arr = numpy.asarray(vector, dtype=float).ravel()
        norm = float(numpy.linalg.norm(arr))
        if norm > 0.0:
            vectors.append(arr / norm)
    if not vectors:
        return []
    nvec = len(vectors)
    overlap = numpy.zeros((nvec, nvec))
    for i in range(nvec):
        for j in range(nvec):
            overlap[i, j] = abs(float(numpy.dot(vectors[i], vectors[j])))
    return _jsonable(overlap)


def _state_energies(solver: Any) -> list[float]:
    for obj in (solver, getattr(solver, "fcisolver", None)):
        if obj is None:
            continue
        for name in ("e_states", "e_roots", "energies"):
            value = getattr(obj, name, None)
            if value is None:
                continue
            if isinstance(value, (list, tuple)):
                return [float(x) for x in value]
            try:
                import numpy

                if isinstance(value, numpy.ndarray):
                    return [float(x) for x in value.ravel()]
            except Exception:
                pass
    return []


def _compute_orbital_grad_norm(solver: Any, fallback: float | None) -> float | None:
    try:
        import numpy

        grad = solver.get_grad(solver.mo_coeff)
        value = float(numpy.linalg.norm(grad))
        return value if math.isfinite(value) else fallback
    except Exception:
        return fallback


def _prepare_case(case: dict[str, Any]) -> tuple[Any, Any, Any]:
    from pyscf import mcscf

    mol = _make_mol(case)
    mf = _make_mf(mol, case)
    base_mc = mcscf.CASSCF(mf, int(case["ncas"]), case["nelecas"])
    base_mc.verbose = int(case.get("casscf_verbose", 0))
    active_sort = case.get("active_orbital_sort")
    if active_sort:
        mo = base_mc.sort_mo(list(active_sort), mo_coeff=mf.mo_coeff, base=int(case.get("active_sort_base", 1)))
    else:
        mo = mf.mo_coeff
    solver = _select_solver(base_mc, mol, case)
    _apply_solver_options(solver, case)
    return mol, solver, mo


def _run_solver_call(solver: Any, mo: Any, family: str, callback: Any) -> tuple[Any, Any, Any, Any, Any]:
    if family == "mc2step":
        return solver.mc2step(mo, callback=callback)
    if family == "mc1step":
        return solver.mc1step(mo, callback=callback)
    return solver.kernel(mo, callback=callback)


def _run_case(case: dict[str, Any], input_root: Path | None) -> dict[str, Any]:
    case = _resolve_case_paths(case, input_root)
    case_id = str(case.get("id") or "case")
    weight = float(case.get("weight", 1.0) or 1.0)
    total_started = time.perf_counter()
    stats: dict[str, Any] = {
        "macro_cycles": 0,
        "micro_cycles": 0,
        "ah_seconds": 0.0,
        "h_op_calls": 0,
        "casci_seconds": 0.0,
        "ao2mo_seconds": 0.0,
        "keyframe_restarts": 0,
        "rejected_steps": 0,
        "norm_gorb": None,
        "norm_gci": None,
    }
    result: dict[str, Any] = {
        "id": case_id,
        "weight": weight,
        "converged": False,
        "wall_seconds": 0.0,
        "total_seconds": 0.0,
        "casscf_kernel_seconds": 0.0,
        "error": "",
    }
    try:
        _, solver, mo = _prepare_case(case)
        family = str(case.get("solver_family") or "newton")
        callback = _make_callback(stats)

        kernel_started = time.perf_counter()
        with _time_instance_method(solver, "casci", "casci_seconds", stats):
            with _time_instance_method(solver, "ao2mo", "ao2mo_seconds", stats):
                with _time_newton_ah(stats):
                    with _time_rotate_orb_cc(solver, stats):
                        e_tot, e_cas, ci, mo_coeff, mo_energy = _run_solver_call(
                            solver, mo, family, callback
                        )
        casscf_seconds = max(0.0, time.perf_counter() - kernel_started)
        norm_gorb = _compute_orbital_grad_norm(solver, _as_float(stats.get("norm_gorb")))
        norm_gci = _as_float(stats.get("norm_gci"))
        result.update(
            {
                "converged": bool(getattr(solver, "converged", False)),
                "e_tot": _as_float(getattr(solver, "e_tot", e_tot)),
                "e_cas": _as_float(getattr(solver, "e_cas", e_cas)),
                "e_states": _state_energies(solver),
                "norm_gorb": norm_gorb,
                "norm_gci": norm_gci,
                "root_overlap": _root_overlap(getattr(solver, "ci", ci)),
                "casscf_kernel_seconds": casscf_seconds,
                "wall_seconds": casscf_seconds,
                "macro_cycles": int(stats.get("macro_cycles", 0)),
                "micro_cycles": int(stats.get("micro_cycles", 0)),
                "ah_seconds": float(stats.get("ah_seconds", 0.0)),
                "h_op_calls": int(stats.get("h_op_calls", 0)),
                "casci_seconds": float(stats.get("casci_seconds", 0.0)),
                "ao2mo_seconds": float(stats.get("ao2mo_seconds", 0.0)),
                "keyframe_restarts": int(stats.get("keyframe_restarts", 0)),
                "rejected_steps": int(stats.get("rejected_steps", 0)),
                "mo_energy": _jsonable(mo_energy),
            }
        )
        if not result["converged"]:
            result["error"] = "CASSCF did not converge"
    except Exception as exc:
        result["converged"] = False
        result["error"] = "".join(traceback.format_exception_only(type(exc), exc)).strip()
    finally:
        total_seconds = max(0.0, time.perf_counter() - total_started)
        result["total_seconds"] = total_seconds
        if not result.get("wall_seconds"):
            result["wall_seconds"] = total_seconds
    return _jsonable(result)


def _primary_summary_metric(primary_metric: str, case_results: list[dict[str, Any]]) -> float:
    if primary_metric.startswith("weighted_median_"):
        field = primary_metric[len("weighted_median_") :]
    else:
        field = "casscf_kernel_seconds"
    return _weighted_median(
        [
            (
                float(case.get(field, case.get("wall_seconds", 0.0)) or 0.0),
                float(case.get("weight", 1.0) or 1.0),
            )
            for case in case_results
        ]
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--benchmark", required=True)
    parser.add_argument("--emit-json", action="store_true")
    args = parser.parse_args()

    benchmark_path = Path(args.benchmark).expanduser().resolve()
    project_root = _find_project_root(benchmark_path)
    os.chdir(project_root)
    sys.path.insert(0, str(project_root))
    _configure_threads()

    payload = yaml.safe_load(benchmark_path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise SystemExit("benchmark file must be a YAML object")
    raw_cases = payload.get("cases")
    cases = [item for item in raw_cases if isinstance(item, dict)] if isinstance(raw_cases, list) else []
    input_root_raw = os.environ.get("FERMILINK_GOAL_INPUT_ROOT")
    input_root = Path(input_root_raw).resolve() if input_root_raw else None
    primary_metric = str(
        (((payload.get("controller") or {}).get("objective") or {}).get("primary_metric"))
        or "weighted_median_casscf_kernel_seconds"
    )

    case_results = [_run_case(case, input_root) for case in cases]
    failures = sum(1 for item in case_results if not bool(item.get("converged")) or bool(item.get("error")))
    summary_metrics = {
        primary_metric: _primary_summary_metric(primary_metric, case_results),
        "weighted_median_casscf_kernel_seconds": _primary_summary_metric(
            "weighted_median_casscf_kernel_seconds", case_results
        ),
        "weighted_median_ah_seconds": _primary_summary_metric("weighted_median_ah_seconds", case_results),
        "weighted_median_casci_seconds": _primary_summary_metric("weighted_median_casci_seconds", case_results),
        "weighted_median_ao2mo_seconds": _primary_summary_metric("weighted_median_ao2mo_seconds", case_results),
        "total_h_op_calls": int(sum(int(item.get("h_op_calls", 0) or 0) for item in case_results)),
        "total_macro_cycles": int(sum(int(item.get("macro_cycles", 0) or 0) for item in case_results)),
        "total_micro_cycles": int(sum(int(item.get("micro_cycles", 0) or 0) for item in case_results)),
        "total_failures": int(failures),
        "peak_rss_mb": _peak_rss_mb(),
    }
    output = {
        "benchmark_id": str(payload.get("benchmark_id") or "pyscf-casscf-ciah-mc1step-autogen-v1"),
        "correctness_ok": failures == 0,
        "summary_metrics": _jsonable(summary_metrics),
        "cases": case_results,
    }
    print(json.dumps(_jsonable(output), sort_keys=True, separators=(",", ":")))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())