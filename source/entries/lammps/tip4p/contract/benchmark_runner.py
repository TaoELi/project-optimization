#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import os
import re
import resource
import shutil
import subprocess
import tempfile
import time
from pathlib import Path
from typing import Any

import yaml

DEFAULT_PRIMARY_METRIC = "weighted_median_pair_plus_kspace_seconds"
FAILURE_PENALTY = 1.0e30
FORCE_SAMPLE_SIZE = 256
TIME_RSS_RE = re.compile(r"FERMILINK_MAXRSS_KB=(\d+)")
LOOP_RE = re.compile(
    r"Loop time of\s+([0-9eE+.\-]+)\s+on\s+(\d+)\s+procs\s+for\s+(\d+)\s+steps\s+with\s+(\d+)\s+atoms"
)


def _find_repo_root(benchmark_path: Path) -> Path:
    start = benchmark_path.resolve().parent
    for candidate in (start, *start.parents):
        if (
            (candidate / ".fermilink-optimize").is_dir()
            and (candidate / "src").is_dir()
            and (candidate / "cmake").is_dir()
        ):
            return candidate
    raise FileNotFoundError(f"unable to locate repository root from {benchmark_path}")


def _load_yaml(path: Path) -> dict[str, Any]:
    payload = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError("benchmark file must be a YAML mapping")
    return payload


def _load_goal_inputs(benchmark_dir: Path) -> dict[str, Any]:
    goal_inputs_path = benchmark_dir / "goal_inputs.json"
    if not goal_inputs_path.is_file():
        return {}
    payload = json.loads(goal_inputs_path.read_text(encoding="utf-8"))
    return payload if isinstance(payload, dict) else {}


def _resolve_relative_to_repo(repo_root: Path, value: str) -> Path:
    path = Path(value).expanduser()
    if path.is_absolute():
        return path.resolve()
    return (repo_root / path).resolve()


def _resolve_input_root(repo_root: Path, benchmark_dir: Path) -> Path:
    env_value = os.environ.get("FERMILINK_GOAL_INPUT_ROOT", "").strip()
    if env_value:
        return _resolve_relative_to_repo(repo_root, env_value)

    goal_inputs = _load_goal_inputs(benchmark_dir)
    all_root_rel = goal_inputs.get("all_root_rel")
    if isinstance(all_root_rel, str) and all_root_rel:
        return _resolve_relative_to_repo(repo_root, all_root_rel)

    goal_dir = goal_inputs.get("goal_dir")
    if isinstance(goal_dir, str) and goal_dir:
        return Path(goal_dir).expanduser().resolve()

    raise FileNotFoundError(
        "unable to resolve staged goal input root; set FERMILINK_GOAL_INPUT_ROOT or provide goal_inputs.json"
    )


def _resolve_lammps_binary(repo_root: Path) -> Path:
    env_value = os.environ.get("FERMILINK_LAMMPS_BIN", "").strip()
    if env_value:
        return _resolve_relative_to_repo(repo_root, env_value)
    return (repo_root / "build" / "lmp").resolve()


def _as_float(value: Any, default: float = 0.0) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def _as_int(value: Any, default: int = 0) -> int:
    try:
        return int(value)
    except (TypeError, ValueError):
        return default


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


def _parse_thermo_table(log_text: str) -> dict[str, Any]:
    tables: list[tuple[list[str], list[list[float]]]] = []
    header: list[str] | None = None
    rows: list[list[float]] = []

    for raw_line in log_text.splitlines():
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
    }


def _parse_loop_and_timing(log_text: str) -> dict[str, float]:
    loop_matches = list(LOOP_RE.finditer(log_text))
    if not loop_matches:
        raise ValueError("missing loop timing")

    loop_match = loop_matches[-1]
    metrics: dict[str, float] = {
        "loop_seconds": float(loop_match.group(1)),
        "logged_mpi_ranks": float(loop_match.group(2)),
        "logged_run_steps": float(loop_match.group(3)),
        "logged_atom_count": float(loop_match.group(4)),
    }

    section_map = {
        "Pair": "pair_seconds",
        "Kspace": "kspace_seconds",
        "Neigh": "neigh_seconds",
        "Comm": "comm_seconds",
    }

    for raw_line in log_text.splitlines():
        if "|" not in raw_line:
            continue
        parts = [part.strip() for part in raw_line.split("|")]
        if len(parts) < 3:
            continue
        section = parts[0]
        metric_name = section_map.get(section)
        if metric_name is None:
            continue
        avg_field = parts[2]
        if avg_field:
            metrics[metric_name] = float(avg_field)

    missing = [name for name in ("pair_seconds", "kspace_seconds", "neigh_seconds", "comm_seconds") if name not in metrics]
    if missing:
        raise ValueError(f"missing timer sections: {', '.join(missing)}")

    return metrics


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


def _extract_peak_rss_mb(stderr_text: str) -> float:
    match = TIME_RSS_RE.search(stderr_text)
    if match:
        return float(match.group(1)) / 1024.0
    return float(resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss) / 1024.0


def _tail_text(text: str, limit: int = 2000) -> str:
    stripped = text.strip()
    if not stripped:
        return ""
    if len(stripped) <= limit:
        return stripped
    return stripped[-limit:]


def _prepare_wrapper_input(
    artifact_dir: Path,
    case_id: str,
    input_script: str,
    force_dump_path: Path,
) -> Path:
    wrapper_path = artifact_dir / f"{case_id}.wrapped.in"
    wrapper_lines = [
        "timer full",
        f"include {input_script}",
        f"write_dump all custom {force_dump_path} id type fx fy fz modify sort id",
    ]
    wrapper_path.write_text("\n".join(wrapper_lines) + "\n", encoding="utf-8")
    return wrapper_path


def _weighted_median(values: list[tuple[float, float]]) -> float:
    filtered = sorted((value, weight) for value, weight in values if weight > 0.0)
    if not filtered:
        return FAILURE_PENALTY

    total_weight = sum(weight for _, weight in filtered)
    midpoint = total_weight / 2.0
    cumulative = 0.0
    for value, weight in filtered:
        cumulative += weight
        if cumulative >= midpoint:
            return value
    return filtered[-1][0]


def _run_case(
    repo_root: Path,
    benchmark_path: Path,
    benchmark_payload: dict[str, Any],
    case: dict[str, Any],
) -> tuple[dict[str, Any], float]:
    benchmark_dir = benchmark_path.resolve().parent
    artifact_root = benchmark_dir / "runner_artifacts"
    artifact_root.mkdir(parents=True, exist_ok=True)

    case_id = str(case.get("id") or "case")
    artifact_dir = Path(tempfile.mkdtemp(prefix=f"{case_id}-", dir=str(artifact_root)))
    input_root = _resolve_input_root(repo_root, benchmark_dir)
    input_script = str(case.get("input_script") or "")
    data_file = str(case.get("data_file") or "")
    expected_atoms = _as_int(case.get("expected_atoms"), 0)
    mpi_ranks = _as_int(case.get("mpi_ranks"), 1)
    omp_num_threads = _as_int(case.get("omp_num_threads"), _as_int(os.environ.get("OMP_NUM_THREADS"), 1))
    requested_steps = _as_int(case.get("run_steps"), 0)

    if not input_script:
        return (
            {
                "id": case_id,
                "converged": False,
                "wall_seconds": 0.0,
                "total_seconds": 0.0,
                "error": "missing case.input_script",
            },
            0.0,
        )

    script_path = input_root / input_script
    data_path = input_root / data_file if data_file else None
    if not script_path.is_file():
        return (
            {
                "id": case_id,
                "converged": False,
                "wall_seconds": 0.0,
                "total_seconds": 0.0,
                "error": f"missing input script: {script_path}",
            },
            0.0,
        )
    if data_path is not None and not data_path.is_file():
        return (
            {
                "id": case_id,
                "converged": False,
                "wall_seconds": 0.0,
                "total_seconds": 0.0,
                "error": f"missing data file: {data_path}",
            },
            0.0,
        )

    lammps_binary = _resolve_lammps_binary(repo_root)
    if not lammps_binary.is_file():
        return (
            {
                "id": case_id,
                "converged": False,
                "wall_seconds": 0.0,
                "total_seconds": 0.0,
                "error": f"missing LAMMPS binary: {lammps_binary}",
            },
            0.0,
        )

    log_path = artifact_dir / f"{case_id}.log"
    force_dump_path = artifact_dir / f"{case_id}.forces.dump"
    wrapper_path = _prepare_wrapper_input(artifact_dir, case_id, input_script, force_dump_path)

    base_command = [
        "mpirun",
        "-np",
        str(mpi_ranks),
        str(lammps_binary),
        "-in",
        str(wrapper_path),
        "-log",
        str(log_path),
        "-screen",
        "none",
    ]

    time_binary = shutil.which("time") or "/usr/bin/time"
    use_external_time = Path(time_binary).is_file() and os.path.basename(time_binary) == "time"
    command = (
        [time_binary, "-f", "FERMILINK_MAXRSS_KB=%M", *base_command]
        if use_external_time
        else base_command
    )

    exec_env = os.environ.copy()
    exec_env["OMP_NUM_THREADS"] = str(omp_num_threads)
    timeout_seconds = _as_float(
        case.get("timeout_seconds"),
        _as_float(
            benchmark_payload.get("controller", {}).get("timeout_seconds"),
            0.0,
        ),
    )
    if timeout_seconds <= 0.0:
        timeout_seconds = None

    started = time.perf_counter()
    try:
        completed = subprocess.run(
            command,
            cwd=str(input_root),
            env=exec_env,
            text=True,
            capture_output=True,
            timeout=timeout_seconds,
            check=False,
        )
    except subprocess.TimeoutExpired:
        elapsed = max(0.0, time.perf_counter() - started)
        return (
            {
                "id": case_id,
                "converged": False,
                "wall_seconds": elapsed,
                "total_seconds": elapsed,
                "error": "timeout",
            },
            0.0,
        )
    except OSError as exc:
        elapsed = max(0.0, time.perf_counter() - started)
        return (
            {
                "id": case_id,
                "converged": False,
                "wall_seconds": elapsed,
                "total_seconds": elapsed,
                "error": str(exc),
            },
            0.0,
        )

    elapsed = max(0.0, time.perf_counter() - started)
    peak_rss_mb = _extract_peak_rss_mb(str(completed.stderr or ""))
    log_text = log_path.read_text(encoding="utf-8", errors="replace") if log_path.is_file() else ""

    if completed.returncode != 0:
        error_text = _tail_text("\n".join(part for part in (completed.stderr, completed.stdout, log_text) if part))
        return (
            {
                "id": case_id,
                "converged": False,
                "wall_seconds": elapsed,
                "total_seconds": elapsed,
                "error": error_text or f"subprocess failed with return code {completed.returncode}",
                "return_code": int(completed.returncode),
            },
            peak_rss_mb,
        )

    try:
        thermo_metrics = _parse_thermo_table(log_text)
        timer_metrics = _parse_loop_and_timing(log_text)
        force_sample_atom_ids, forces_xyz, dumped_atoms = _parse_force_dump(force_dump_path)
    except Exception as exc:
        error_text = _tail_text(f"{exc}\n{log_text}")
        return (
            {
                "id": case_id,
                "converged": False,
                "wall_seconds": elapsed,
                "total_seconds": elapsed,
                "error": error_text or str(exc),
            },
            peak_rss_mb,
        )

    atom_count = int(thermo_metrics["atom_count"])
    if dumped_atoms != atom_count:
        return (
            {
                "id": case_id,
                "converged": False,
                "wall_seconds": elapsed,
                "total_seconds": elapsed,
                "error": f"atom count mismatch between thermo ({atom_count}) and force dump ({dumped_atoms})",
            },
            peak_rss_mb,
        )
    if expected_atoms and atom_count != expected_atoms:
        return (
            {
                "id": case_id,
                "converged": False,
                "wall_seconds": elapsed,
                "total_seconds": elapsed,
                "error": f"unexpected atom count: expected {expected_atoms}, observed {atom_count}",
            },
            peak_rss_mb,
        )

    loop_seconds = float(timer_metrics["loop_seconds"])
    if loop_seconds <= 0.0:
        return (
            {
                "id": case_id,
                "converged": False,
                "wall_seconds": elapsed,
                "total_seconds": elapsed,
                "error": "non-positive loop_seconds",
            },
            peak_rss_mb,
        )

    run_steps = requested_steps or int(timer_metrics["logged_run_steps"])
    pair_plus_kspace_seconds = float(timer_metrics["pair_seconds"]) + float(timer_metrics["kspace_seconds"])
    steps_per_second = float(run_steps) / loop_seconds if run_steps > 0 else 0.0

    result = {
        "id": case_id,
        "converged": True,
        "wall_seconds": elapsed,
        "total_seconds": elapsed,
        "loop_seconds": loop_seconds,
        "pair_seconds": float(timer_metrics["pair_seconds"]),
        "kspace_seconds": float(timer_metrics["kspace_seconds"]),
        "pair_plus_kspace_seconds": pair_plus_kspace_seconds,
        "neigh_seconds": float(timer_metrics["neigh_seconds"]),
        "comm_seconds": float(timer_metrics["comm_seconds"]),
        "steps_per_second": steps_per_second,
        "thermo_steps": thermo_metrics["thermo_steps"],
        "atom_count": atom_count,
        "etotal_series": thermo_metrics["etotal_series"],
        "energy_drift_per_atom_per_step": thermo_metrics["energy_drift_per_atom_per_step"],
        "pe_series": thermo_metrics["pe_series"],
        "ke_series": thermo_metrics["ke_series"],
        "temp_series": thermo_metrics["temp_series"],
        "press_series": thermo_metrics["press_series"],
        "density_series": thermo_metrics["density_series"],
        "force_sample_atom_ids": force_sample_atom_ids,
        "forces_xyz": forces_xyz,
        "error": "",
    }
    return result, peak_rss_mb


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--benchmark", required=True)
    parser.add_argument("--emit-json", action="store_true")
    args = parser.parse_args()

    benchmark_path = Path(args.benchmark).expanduser().resolve()
    primary_metric_name = DEFAULT_PRIMARY_METRIC

    try:
        benchmark_payload = _load_yaml(benchmark_path)
        repo_root = _find_repo_root(benchmark_path)
        primary_metric_name = str(
            benchmark_payload.get("controller", {})
            .get("objective", {})
            .get("primary_metric", DEFAULT_PRIMARY_METRIC)
        )
        cases_raw = benchmark_payload.get("cases")
        if not isinstance(cases_raw, list) or not cases_raw:
            raise ValueError("benchmark file requires a non-empty cases list")
        cases = [item for item in cases_raw if isinstance(item, dict)]
        if not cases:
            raise ValueError("benchmark file requires dict entries in cases")

        case_results: list[dict[str, Any]] = []
        weighted_values: list[tuple[float, float]] = []
        peak_rss_mb = 0.0
        total_failures = 0

        for case in cases:
            result, case_peak_rss = _run_case(repo_root, benchmark_path, benchmark_payload, case)
            case_results.append(result)
            peak_rss_mb = max(peak_rss_mb, case_peak_rss)
            if result.get("converged"):
                weighted_values.append(
                    (
                        float(result["pair_plus_kspace_seconds"]),
                        _as_float(case.get("weight"), 1.0),
                    )
                )
            else:
                total_failures += 1

        primary_metric_value = (
            _weighted_median(weighted_values)
            if total_failures == 0 and len(weighted_values) == len(cases)
            else FAILURE_PENALTY
        )

        output = {
            "benchmark_id": str(benchmark_payload.get("benchmark_id") or "benchmark"),
            "correctness_ok": total_failures == 0,
            "summary_metrics": {
                primary_metric_name: float(primary_metric_value),
                "peak_rss_mb": float(peak_rss_mb),
                "total_failures": int(total_failures),
            },
            "cases": case_results,
        }
    except Exception as exc:
        output = {
            "benchmark_id": benchmark_path.stem,
            "correctness_ok": False,
            "summary_metrics": {
                primary_metric_name: FAILURE_PENALTY,
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

    print(json.dumps(output, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())