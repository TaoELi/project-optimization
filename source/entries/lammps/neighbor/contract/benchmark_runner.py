#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

import yaml


PRIMARY_METRIC = "weighted_median_neigh_seconds"
FAILURE_PENALTY = 1.0e30
FORCE_SAMPLE_SIZE = 256
WRAPPER_DIRNAME = "runtime_wrappers"
ROOT_SENTINELS = (".fermilink-optimize", "src", "cmake")
CASE_REQUIRED_FIELDS = (
    "loop_seconds",
    "neigh_seconds",
    "pair_seconds",
    "kspace_seconds",
    "comm_seconds",
    "bond_seconds",
    "timesteps_per_second",
    "thermo_steps",
    "atom_count",
    "etotal_series",
    "pe_series",
    "ke_series",
    "temp_series",
    "press_series",
    "density_series",
    "force_sample_atom_ids",
    "forces_xyz",
    "neighbor_list_builds",
    "dangerous_builds",
)
TIMING_LABEL_MAP = {
    "Pair": "pair_seconds",
    "Bond": "bond_seconds",
    "Kspace": "kspace_seconds",
    "Neigh": "neigh_seconds",
    "Comm": "comm_seconds",
}
LOOP_RE = re.compile(
    r"Loop time of\s+([0-9eE+\-.]+)\s+on\s+(\d+)\s+procs\s+for\s+(\d+)\s+steps\s+with\s+(\d+)\s+atoms"
)
PERF_RE = re.compile(
    r"Performance:\s+([0-9eE+\-.]+)\s+([A-Za-z/]+),.*?([0-9eE+\-.]+)\s+timesteps/s(?:,\s+([0-9eE+\-.]+)\s+([A-Za-z\-\/]+))?"
)
NEIGH_BUILDS_RE = re.compile(r"Neighbor list builds =\s+([0-9eE+\-.]+)")
DANGEROUS_BUILDS_RE = re.compile(r"Dangerous builds =\s+([0-9eE+\-.]+)")
RSS_RE = re.compile(r"Maximum resident set size \(kbytes\):\s+(\d+)")


def _read_yaml(path: Path) -> dict[str, Any]:
    payload = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"benchmark file {path} must contain a YAML object")
    return payload


def _read_json(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"JSON file {path} must contain an object")
    return payload


def _discover_repo_root(benchmark_path: Path) -> Path:
    candidates: list[Path] = [Path.cwd(), benchmark_path.parent]
    candidates.extend(benchmark_path.parents)
    seen: set[Path] = set()
    for candidate in candidates:
        resolved = candidate.resolve()
        if resolved in seen:
            continue
        seen.add(resolved)
        if all((resolved / name).exists() for name in ROOT_SENTINELS):
            return resolved
    raise FileNotFoundError(f"could not discover repository root from {benchmark_path}")


def _runtime_env_defaults(payload: dict[str, Any]) -> dict[str, str]:
    runtime = payload.get("runtime")
    if not isinstance(runtime, dict):
        return {}
    env_raw = runtime.get("env")
    if not isinstance(env_raw, dict):
        return {}
    return {str(key): str(value) for key, value in env_raw.items()}


def _goal_input_candidates(repo_root: Path, base_env: dict[str, str]) -> list[Path]:
    candidates: list[Path] = []

    def add_candidate(raw: str | None) -> None:
        if not raw:
            return
        path = Path(raw)
        if not path.is_absolute():
            path = repo_root / path
        candidates.append(path.resolve())

    add_candidate(os.environ.get("FERMILINK_GOAL_INPUT_ROOT"))
    add_candidate(base_env.get("FERMILINK_GOAL_INPUT_ROOT"))

    goal_inputs_path = repo_root / ".fermilink-optimize" / "autogen" / "goal_inputs.json"
    if goal_inputs_path.is_file():
        try:
            goal_inputs = _read_json(goal_inputs_path)
        except Exception:
            goal_inputs = {}
        for key in ("worker_root_rel", "all_root_rel"):
            rel_path = goal_inputs.get(key)
            if isinstance(rel_path, str):
                add_candidate(rel_path)

    deduped: list[Path] = []
    seen: set[Path] = set()
    for candidate in candidates:
        if candidate in seen:
            continue
        seen.add(candidate)
        deduped.append(candidate)
    return deduped


def _resolve_case_input_root(repo_root: Path, case: dict[str, Any], base_env: dict[str, str]) -> Path:
    required_files = [str(case.get("input_script") or ""), str(case.get("data_file") or "")]
    if not all(required_files):
        raise ValueError(f"case {case.get('id', 'case')} is missing input_script or data_file")

    for candidate in _goal_input_candidates(repo_root, base_env):
        if candidate.is_dir() and all((candidate / rel).is_file() for rel in required_files):
            return candidate

    joined = ", ".join(required_files)
    raise FileNotFoundError(f"could not resolve staged input root containing {joined}")


def _find_mpi_launcher() -> list[str]:
    env_value = os.environ.get("FERMILINK_MPI_LAUNCHER") or os.environ.get("MPI_LAUNCHER")
    if env_value:
        launcher = env_value.strip().split()
        if launcher:
            return launcher
    for program in ("mpirun", "mpiexec", "srun"):
        if shutil.which(program):
            return [program]
    raise FileNotFoundError("no MPI launcher found; set FERMILINK_MPI_LAUNCHER to mpirun/mpiexec/srun")


def _mpi_command(launcher: list[str], ranks: int) -> list[str]:
    if not launcher:
        raise ValueError("launcher must not be empty")
    head = Path(launcher[0]).name
    if head == "srun":
        return [*launcher, "-n", str(ranks)]
    return [*launcher, "-np", str(ranks)]


def _timed_command(command: list[str]) -> list[str]:
    time_bin = shutil.which("/usr/bin/time") or shutil.which("time")
    if time_bin:
        return [time_bin, "-v", *command]
    return command


def _artifacts_for_case(repo_root: Path, case_id: str) -> tuple[Path, Path]:
    wrapper_dir = repo_root / ".fermilink-optimize" / "autogen" / WRAPPER_DIRNAME
    wrapper_dir.mkdir(parents=True, exist_ok=True)
    safe_case = re.sub(r"[^A-Za-z0-9_.-]+", "-", case_id) or "case"
    unique = f"{os.getpid()}-{time.time_ns()}"
    wrapper_path = wrapper_dir / f"{safe_case}-{unique}.in"
    force_dump_path = wrapper_dir / f"{safe_case}-{unique}.force.dump"
    return wrapper_path, force_dump_path


def _write_wrapper_input(wrapper_path: Path, input_script: Path, force_dump_path: Path) -> None:
    wrapper_text = "\n".join(
        (
            "timer full",
            f"include {input_script.name}",
            f"write_dump all custom {force_dump_path} id fx fy fz modify sort id",
        )
    )
    wrapper_path.write_text(wrapper_text + "\n", encoding="utf-8")


def _try_float(text: str) -> float | None:
    try:
        return float(text)
    except (TypeError, ValueError):
        return None


def _is_number(token: str) -> bool:
    try:
        float(token)
        return True
    except ValueError:
        return False


def _normalize_header(name: str) -> str:
    return "".join(ch for ch in name.lower() if ch.isalnum())


def _column_index(header: list[str], aliases: list[str]) -> int:
    normalized = {_normalize_header(name): idx for idx, name in enumerate(header)}
    for alias in aliases:
        match = normalized.get(alias)
        if match is not None:
            return match
    raise KeyError(f"missing thermo column for aliases {aliases!r}")


def _parse_loop_metrics(stdout_text: str) -> dict[str, Any]:
    metrics: dict[str, Any] = {}
    for line in stdout_text.splitlines():
        match = LOOP_RE.search(line)
        if match:
            metrics["loop_seconds"] = float(match.group(1))
            metrics["mpi_tasks"] = int(match.group(2))
            metrics["run_steps"] = int(match.group(3))
            metrics["atoms"] = int(match.group(4))
        match = PERF_RE.search(line)
        if match:
            metrics["throughput_per_day"] = float(match.group(1))
            metrics["throughput_per_day_unit"] = match.group(2)
            metrics["timesteps_per_second"] = float(match.group(3))
            atom_step_value = _try_float(match.group(4) or "")
            atom_step_unit = match.group(5)
            if atom_step_value is not None and atom_step_unit:
                metrics["atom_steps_per_second"] = atom_step_value
                metrics["atom_steps_per_second_unit"] = atom_step_unit
    if "timesteps_per_second" not in metrics and "loop_seconds" in metrics and metrics["loop_seconds"] > 0:
        run_steps = metrics.get("run_steps")
        if isinstance(run_steps, int) and run_steps > 0:
            metrics["timesteps_per_second"] = run_steps / metrics["loop_seconds"]
    return metrics


def _parse_timing_rows(stdout_text: str) -> dict[str, float]:
    metrics: dict[str, float] = {}
    for raw_line in stdout_text.splitlines():
        if "|" not in raw_line:
            continue
        parts = [part.strip() for part in raw_line.split("|")]
        if len(parts) < 3:
            continue
        label = parts[0]
        field_name = TIMING_LABEL_MAP.get(label)
        if not field_name:
            continue
        avg_time = _try_float(parts[2].split()[0] if parts[2] else "")
        if avg_time is None:
            continue
        metrics[field_name] = avg_time
    return metrics


def _parse_build_counters(stdout_text: str) -> dict[str, int]:
    metrics: dict[str, int] = {}
    neigh_match = NEIGH_BUILDS_RE.search(stdout_text)
    if neigh_match:
        metrics["neighbor_list_builds"] = int(float(neigh_match.group(1)))
    dangerous_match = DANGEROUS_BUILDS_RE.search(stdout_text)
    if dangerous_match:
        metrics["dangerous_builds"] = int(float(dangerous_match.group(1)))
    return metrics


def _parse_thermo_table(stdout_text: str) -> dict[str, Any]:
    tables: list[tuple[list[str], list[list[float]]]] = []
    header: list[str] | None = None
    rows: list[list[float]] = []

    for raw_line in stdout_text.splitlines():
        stripped = raw_line.strip()
        if not stripped:
            continue

        tokens = stripped.split()
        if "Step" in tokens and any(token in tokens for token in ("Temp", "TotEng", "Press", "Density", "Atoms")):
            if header and rows:
                tables.append((header, rows))
            header = tokens
            rows = []
            continue

        if header is None:
            continue

        if stripped.startswith("Loop time of"):
            if rows:
                tables.append((header, rows))
            header = None
            rows = []
            continue

        if len(tokens) == len(header) and all(_is_number(token) for token in tokens):
            rows.append([float(token) for token in tokens])

    if header and rows:
        tables.append((header, rows))

    if not tables:
        raise ValueError("missing thermo table")

    header, rows = tables[-1]
    step_idx = _column_index(header, ["step"])
    atoms_idx = _column_index(header, ["atoms"])
    temp_idx = _column_index(header, ["temp"])
    etotal_idx = _column_index(header, ["toteng", "etotal"])
    pe_idx = _column_index(header, ["poteng", "pe"])
    ke_idx = _column_index(header, ["kineng", "ke"])
    press_idx = _column_index(header, ["press"])
    density_idx = _column_index(header, ["density"])

    steps = [int(round(row[step_idx])) for row in rows]
    atoms = [int(round(row[atoms_idx])) for row in rows]
    temp_series = [float(row[temp_idx]) for row in rows]
    etotal_series = [float(row[etotal_idx]) for row in rows]
    pe_series = [float(row[pe_idx]) for row in rows]
    ke_series = [float(row[ke_idx]) for row in rows]
    press_series = [float(row[press_idx]) for row in rows]
    density_series = [float(row[density_idx]) for row in rows]

    atom_count = atoms[-1] if atoms else 0
    if len(steps) >= 2 and atom_count > 0 and steps[-1] != steps[0]:
        energy_drift = (etotal_series[-1] - etotal_series[0]) / (atom_count * (steps[-1] - steps[0]))
    else:
        energy_drift = 0.0

    return {
        "thermo_steps": steps,
        "atom_count": atom_count,
        "etotal_series": etotal_series,
        "energy_drift_per_atom_per_step": float(energy_drift),
        "pe_series": pe_series,
        "ke_series": ke_series,
        "temp_series": temp_series,
        "press_series": press_series,
        "density_series": density_series,
        "etotal": etotal_series[-1],
        "pe": pe_series[-1],
        "ke": ke_series[-1],
        "temp": temp_series[-1],
        "press": press_series[-1],
        "density": density_series[-1],
    }


def _sample_indices(total_count: int, sample_size: int) -> list[int]:
    if total_count <= sample_size:
        return list(range(total_count))
    if sample_size <= 1:
        return [0]
    return sorted({(index * (total_count - 1)) // (sample_size - 1) for index in range(sample_size)})


def _parse_force_dump(force_dump_path: Path) -> tuple[list[int], list[float], int]:
    lines = force_dump_path.read_text(encoding="utf-8", errors="replace").splitlines()
    atom_lines: list[str] = []
    field_names: list[str] = []
    atom_count = 0

    index = 0
    while index < len(lines):
        stripped = lines[index].strip()
        if stripped == "ITEM: NUMBER OF ATOMS" and index + 1 < len(lines):
            atom_count = int(lines[index + 1].strip())
            index += 2
            continue
        if stripped.startswith("ITEM: ATOMS"):
            field_names = stripped.split()[2:]
            atom_lines = lines[index + 1 :]
            break
        index += 1

    if not atom_lines or not field_names:
        raise ValueError(f"missing atom section in force dump {force_dump_path}")

    try:
        id_idx = field_names.index("id")
        fx_idx = field_names.index("fx")
        fy_idx = field_names.index("fy")
        fz_idx = field_names.index("fz")
    except ValueError as exc:
        raise ValueError(f"force dump missing required columns: {field_names}") from exc

    rows: list[tuple[int, float, float, float]] = []
    for raw_line in atom_lines:
        stripped = raw_line.strip()
        if not stripped or stripped.startswith("ITEM:"):
            break
        tokens = stripped.split()
        rows.append(
            (
                int(tokens[id_idx]),
                float(tokens[fx_idx]),
                float(tokens[fy_idx]),
                float(tokens[fz_idx]),
            )
        )

    if not rows:
        raise ValueError(f"force dump {force_dump_path} contains no atom rows")

    rows.sort(key=lambda item: item[0])
    selected = [rows[idx] for idx in _sample_indices(len(rows), FORCE_SAMPLE_SIZE)]

    atom_ids: list[int] = []
    force_components: list[float] = []
    for atom_id, fx, fy, fz in selected:
        atom_ids.append(atom_id)
        force_components.extend((fx, fy, fz))

    if atom_count <= 0:
        atom_count = len(rows)
    return atom_ids, force_components, atom_count


def _parse_peak_rss_mb(stderr_text: str) -> float:
    match = RSS_RE.search(stderr_text)
    if not match:
        return 0.0
    return float(match.group(1)) / 1024.0


def _extract_error_text(stdout_text: str, stderr_text: str) -> str:
    chunks = [text.strip() for text in (stderr_text, stdout_text) if text and text.strip()]
    if not chunks:
        return "process failed without output"
    lines = "\n".join(chunks).splitlines()
    return "\n".join(lines[-40:])


def _weighted_median(values: list[tuple[float, float]]) -> float:
    ordered = sorted((value, weight) for value, weight in values if weight > 0.0)
    total_weight = sum(weight for _, weight in ordered)
    if total_weight <= 0:
        raise ValueError("weighted median requires positive total weight")
    threshold = total_weight / 2.0
    cumulative = 0.0
    for value, weight in ordered:
        cumulative += weight
        if cumulative >= threshold:
            return value
    return ordered[-1][0]


def _primary_metric_name(payload: dict[str, Any]) -> str:
    controller = payload.get("controller")
    if not isinstance(controller, dict):
        return PRIMARY_METRIC
    objective = controller.get("objective")
    if not isinstance(objective, dict):
        return PRIMARY_METRIC
    raw_name = objective.get("primary_metric")
    if isinstance(raw_name, str) and raw_name.strip():
        return raw_name.strip()
    return PRIMARY_METRIC


def _run_case(
    repo_root: Path,
    payload: dict[str, Any],
    case: dict[str, Any],
    base_env: dict[str, str],
    launcher: list[str],
) -> dict[str, Any]:
    case_id = str(case.get("id") or "case")
    result: dict[str, Any] = {
        "id": case_id,
        "converged": False,
        "wall_seconds": 0.0,
        "total_seconds": 0.0,
        "error": "",
    }
    started = time.perf_counter()
    wrapper_path: Path | None = None
    force_dump_path: Path | None = None
    try:
        input_root = _resolve_case_input_root(repo_root, case, base_env)
        input_script = input_root / str(case["input_script"])
        if not input_script.is_file():
            raise FileNotFoundError(f"missing input script {input_script}")
        data_file = input_root / str(case["data_file"])
        if not data_file.is_file():
            raise FileNotFoundError(f"missing data file {data_file}")

        lmp_path = repo_root / "build" / "lmp"
        if not lmp_path.is_file():
            raise FileNotFoundError(f"missing LAMMPS executable {lmp_path}")

        mpi_ranks = int(case.get("mpi_ranks") or 0)
        if mpi_ranks <= 0:
            raise ValueError(f"case {case_id} must specify a positive mpi_ranks")

        wrapper_path, force_dump_path = _artifacts_for_case(repo_root, case_id)
        _write_wrapper_input(wrapper_path, input_script, force_dump_path)

        command = _timed_command(
            [
                *_mpi_command(launcher, mpi_ranks),
                str(lmp_path),
                "-log",
                "none",
                "-in",
                str(wrapper_path),
            ]
        )

        child_env = os.environ.copy()
        for key, value in base_env.items():
            child_env.setdefault(key, value)
        if case.get("omp_threads") is not None:
            child_env["OMP_NUM_THREADS"] = str(case["omp_threads"])

        controller = payload.get("controller")
        timeout_seconds: float | None = None
        if isinstance(controller, dict):
            raw_timeout = controller.get("timeout_seconds")
            if isinstance(raw_timeout, (int, float)) and float(raw_timeout) > 0:
                timeout_seconds = float(raw_timeout)

        completed = subprocess.run(
            command,
            cwd=str(input_root),
            env=child_env,
            text=True,
            capture_output=True,
            timeout=timeout_seconds,
            check=False,
        )
        elapsed = max(0.0, time.perf_counter() - started)
        result["wall_seconds"] = elapsed
        result["return_code"] = int(completed.returncode)

        stdout_text = completed.stdout or ""
        stderr_text = completed.stderr or ""

        if completed.returncode != 0:
            result["total_seconds"] = elapsed
            result["error"] = _extract_error_text(stdout_text, stderr_text)
            return result

        if force_dump_path is None or not force_dump_path.is_file():
            raise FileNotFoundError(f"missing force dump {force_dump_path}")

        parsed_metrics: dict[str, Any] = {}
        parsed_metrics.update(_parse_loop_metrics(stdout_text))
        parsed_metrics.update(_parse_timing_rows(stdout_text))
        parsed_metrics.update(_parse_build_counters(stdout_text))
        parsed_metrics.update(_parse_thermo_table(stdout_text))

        force_sample_atom_ids, forces_xyz, dumped_atoms = _parse_force_dump(force_dump_path)
        parsed_metrics["force_sample_atom_ids"] = force_sample_atom_ids
        parsed_metrics["forces_xyz"] = forces_xyz
        parsed_metrics["peak_rss_mb"] = _parse_peak_rss_mb(stderr_text)

        atom_count = int(parsed_metrics.get("atom_count") or 0)
        if atom_count <= 0:
            raise ValueError("non-positive atom_count")
        if dumped_atoms != atom_count:
            raise ValueError(
                f"atom count mismatch between thermo ({atom_count}) and force dump ({dumped_atoms})"
            )

        if "loop_seconds" in parsed_metrics:
            result["total_seconds"] = float(parsed_metrics["loop_seconds"])
        else:
            result["total_seconds"] = elapsed

        for key, value in parsed_metrics.items():
            result[key] = value

        missing_fields = [field for field in CASE_REQUIRED_FIELDS if field not in result]
        if missing_fields:
            result["error"] = "missing parsed fields: " + ", ".join(missing_fields)
        elif int(result.get("dangerous_builds", 1)) != 0:
            result["error"] = f"dangerous builds reported: {result['dangerous_builds']}"
        else:
            result["converged"] = True
            result["error"] = ""
        return result
    except subprocess.TimeoutExpired as exc:
        elapsed = max(0.0, time.perf_counter() - started)
        result["wall_seconds"] = elapsed
        result["total_seconds"] = elapsed
        result["return_code"] = 124
        result["error"] = f"timeout after {exc.timeout} seconds"
        return result
    except Exception as exc:
        elapsed = max(0.0, time.perf_counter() - started)
        result["wall_seconds"] = elapsed
        result["total_seconds"] = elapsed
        result["return_code"] = 1
        result["error"] = str(exc)
        return result
    finally:
        for path in (wrapper_path, force_dump_path):
            if path is None:
                continue
            try:
                path.unlink(missing_ok=True)
            except OSError:
                pass


def _run_benchmark(benchmark_path: Path) -> dict[str, Any]:
    payload = _read_yaml(benchmark_path)
    benchmark_id = str(payload.get("benchmark_id") or benchmark_path.stem)
    primary_metric_name = _primary_metric_name(payload)
    repo_root = _discover_repo_root(benchmark_path)

    cases_raw = payload.get("cases")
    if not isinstance(cases_raw, list) or not cases_raw:
        raise ValueError("benchmark YAML must contain a non-empty cases list")
    cases = [case for case in cases_raw if isinstance(case, dict)]
    if not cases:
        raise ValueError("benchmark YAML cases must be objects")

    base_env = _runtime_env_defaults(payload)
    launcher = _find_mpi_launcher()

    case_results = [_run_case(repo_root, payload, case, base_env, launcher) for case in cases]
    good_neigh_values: list[tuple[float, float]] = []
    good_loop_values: list[tuple[float, float]] = []
    failures = 0
    peak_rss_mb = 0.0

    for case, result in zip(cases, case_results):
        peak_rss_mb = max(peak_rss_mb, float(result.get("peak_rss_mb") or 0.0))
        weight = float(case.get("weight") or 1.0)
        if result.get("converged"):
            if isinstance(result.get("neigh_seconds"), (int, float)):
                good_neigh_values.append((float(result["neigh_seconds"]), weight))
            if isinstance(result.get("loop_seconds"), (int, float)):
                good_loop_values.append((float(result["loop_seconds"]), weight))
        else:
            failures += 1

    if good_neigh_values and failures == 0:
        primary_metric = _weighted_median(good_neigh_values)
    else:
        primary_metric = FAILURE_PENALTY

    if good_loop_values and failures == 0:
        weighted_median_loop = _weighted_median(good_loop_values)
    else:
        weighted_median_loop = FAILURE_PENALTY

    correctness_ok = failures == 0 and all(int(result.get("dangerous_builds", 0)) == 0 for result in case_results)
    return {
        "benchmark_id": benchmark_id,
        "correctness_ok": correctness_ok,
        "summary_metrics": {
            primary_metric_name: float(primary_metric),
            "weighted_median_loop_seconds": float(weighted_median_loop),
            "peak_rss_mb": float(peak_rss_mb),
            "total_failures": int(failures),
        },
        "cases": case_results,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--benchmark", required=True)
    parser.add_argument("--emit-json", action="store_true")
    args = parser.parse_args()

    benchmark_path = Path(args.benchmark).resolve()
    try:
        payload = _run_benchmark(benchmark_path)
    except Exception as exc:
        payload = {
            "benchmark_id": benchmark_path.stem,
            "correctness_ok": False,
            "summary_metrics": {
                PRIMARY_METRIC: FAILURE_PENALTY,
                "peak_rss_mb": 0.0,
            },
            "cases": [
                {
                    "id": "benchmark",
                    "converged": False,
                    "wall_seconds": 0.0,
                    "total_seconds": 0.0,
                    "error": str(exc),
                }
            ],
        }

    if args.emit_json:
        json.dump(payload, sys.stdout, indent=2, sort_keys=True)
        sys.stdout.write("\n")
    else:
        sys.stdout.write(json.dumps(payload) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())