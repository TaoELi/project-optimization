#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import math
import os
import re
import shutil
import statistics
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

import yaml


FAILURE_PENALTY_METRIC = 1.0e9


def _as_float(value: Any) -> float | None:
    try:
        parsed = float(value)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(parsed):
        return None
    return parsed


def _as_int(value: Any) -> int | None:
    try:
        return int(str(value).strip())
    except (TypeError, ValueError):
        return None


def _is_number(token: str) -> bool:
    try:
        float(token)
        return True
    except ValueError:
        return False


def _find_project_root(start: Path) -> Path:
    for candidate in [start, *start.parents]:
        has_git = (candidate / ".git").exists()
        has_src = (candidate / "src").is_dir()
        if has_git and has_src:
            return candidate
    return start.resolve()


def _resolve_input_root(project_root: Path) -> Path:
    env_value = os.environ.get("FERMILINK_GOAL_INPUT_ROOT", "").strip()
    if env_value:
        env_path = Path(env_value).expanduser()
        if not env_path.is_absolute():
            env_path = (project_root / env_path).resolve()
        return env_path
    fallback = (project_root / ".fermilink-optimize" / "inputs" / "all").resolve()
    if fallback.is_dir():
        return fallback
    return project_root.resolve()


def _resolve_case_file(raw_path: str, input_root: Path, project_root: Path) -> Path:
    path = Path(raw_path).expanduser()
    if path.is_absolute():
        return path
    in_input = (input_root / path).resolve()
    if in_input.exists():
        return in_input
    return (project_root / path).resolve()


def _command_path(path: Path, cwd: Path) -> str:
    try:
        rel = path.resolve().relative_to(cwd.resolve())
        return str(rel)
    except ValueError:
        return str(path.resolve())


def _choose_mpi_launcher() -> str | None:
    for candidate in ("mpiexec", "mpirun"):
        if shutil.which(candidate):
            return candidate
    return None


def _choose_time_wrapper() -> tuple[list[str], str]:
    time_bin = Path("/usr/bin/time")
    if not time_bin.exists():
        return ([], "")
    if sys.platform == "darwin":
        return ([str(time_bin), "-l"], "bsd")
    return ([str(time_bin), "-v"], "gnu")


def _parse_max_rss_mb(stderr_text: str, mode: str) -> float | None:
    if not stderr_text:
        return None
    if mode == "gnu":
        match = re.search(r"Maximum resident set size \(kbytes\):\s*([0-9]+)", stderr_text)
        if match:
            return float(int(match.group(1)) / 1024.0)
    if mode == "bsd":
        match = re.search(r"^\s*([0-9]+)\s+maximum resident set size", stderr_text, flags=re.IGNORECASE | re.MULTILINE)
        if match:
            value = int(match.group(1))
            if value > 10_000_000:
                return float(value / (1024.0 * 1024.0))
            return float(value / 1024.0)
    generic = re.search(r"resident set size[^0-9]*([0-9]+)", stderr_text, flags=re.IGNORECASE)
    if generic:
        return float(int(generic.group(1)) / 1024.0)
    return None


def _parse_input_script_metadata(script_path: Path) -> dict[str, float | int]:
    metadata: dict[str, float | int] = {}
    try:
        text = script_path.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return metadata
    for line in text.splitlines():
        body = line.split("#", 1)[0].strip()
        if not body:
            continue
        tokens = body.split()
        if not tokens:
            continue
        cmd = tokens[0].lower()
        if cmd == "timestep" and len(tokens) >= 2:
            timestep = _as_float(tokens[1])
            if timestep is not None:
                metadata["timestep_fs"] = timestep
        elif cmd == "run" and len(tokens) >= 2:
            run_steps = _as_int(tokens[1])
            if run_steps is not None:
                metadata["run_steps"] = run_steps
    return metadata


def _extract_last_thermo_row(text: str) -> dict[str, float]:
    selected_header: list[str] | None = None
    selected_row: list[str] | None = None
    lines = text.splitlines()
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        tokens = line.split()
        if tokens and tokens[0].lower() == "step" and len(tokens) >= 2:
            header = tokens
            j = i + 1
            last_row: list[str] | None = None
            while j < len(lines):
                row = lines[j].strip()
                if not row:
                    if last_row is not None:
                        break
                    j += 1
                    continue
                if row.startswith("Loop time of") or row.startswith("ERROR"):
                    break
                row_tokens = row.split()
                if len(row_tokens) == len(header) and all(_is_number(token) for token in row_tokens):
                    last_row = row_tokens
                    j += 1
                    continue
                if last_row is None:
                    j += 1
                    continue
                break
            if last_row is not None:
                selected_header = header
                selected_row = last_row
            i = j
            continue
        i += 1

    if selected_header is None or selected_row is None:
        return {}

    values: dict[str, float] = {}
    for key, token in zip(selected_header, selected_row):
        parsed = _as_float(token)
        if parsed is not None:
            values[key.lower()] = parsed
    return values


def _parse_loop_metrics(text: str) -> dict[str, float | int]:
    output: dict[str, float | int] = {}
    loop_matches = list(
        re.finditer(
            r"Loop time of\s+([0-9eE+\-.]+)\s+on\s+\d+\s+procs?\s+for\s+(\d+)\s+steps?\s+with\s+(\d+)\s+atoms",
            text,
        )
    )
    if loop_matches:
        last = loop_matches[-1]
        output["loop_time_seconds"] = float(last.group(1))
        output["steps"] = int(last.group(2))
        output["loop_atoms"] = int(last.group(3))

    perf_matches = list(
        re.finditer(r"Performance:\s+([0-9eE+\-.]+)\s+ns/day.*?([0-9eE+\-.]+)\s+timesteps/s", text)
    )
    if perf_matches:
        perf = perf_matches[-1]
        output["ns_per_day"] = float(perf.group(1))
        output["timesteps_per_second"] = float(perf.group(2))

    return output


def _pick_first(mapping: dict[str, float], keys: list[str]) -> float | None:
    for key in keys:
        if key in mapping:
            return mapping[key]
    return None


def _count_nan_inf(text: str, numeric_values: list[float]) -> int:
    count = len(re.findall(r"(?<![A-Za-z0-9_])[+-]?(?:nan|inf|infinity)(?![A-Za-z0-9_])", text, flags=re.IGNORECASE))
    for value in numeric_values:
        if not math.isfinite(value):
            count += 1
    return count


def _parse_case_outputs(full_text: str, input_meta: dict[str, float | int]) -> dict[str, Any]:
    thermo = _extract_last_thermo_row(full_text)
    loop = _parse_loop_metrics(full_text)

    outputs: dict[str, Any] = {}
    etotal = _pick_first(thermo, ["etotal", "toteng", "e_total", "etot"])
    pe = _pick_first(thermo, ["pe", "poteng"])
    ke = _pick_first(thermo, ["ke", "kineng"])
    press = _pick_first(thermo, ["press", "pressure"])
    atoms = _pick_first(thermo, ["atoms"])

    if etotal is not None:
        outputs["thermo.etotal"] = float(etotal)
    if pe is not None:
        outputs["thermo.pe"] = float(pe)
    if ke is not None:
        outputs["thermo.ke"] = float(ke)
    if press is not None:
        outputs["thermo.press"] = float(press)

    if atoms is None and "loop_atoms" in loop:
        atoms = float(loop["loop_atoms"])
    if atoms is not None:
        outputs["atom_count"] = int(round(atoms))

    if "loop_time_seconds" in loop:
        outputs["loop_time_seconds"] = float(loop["loop_time_seconds"])
    if "steps" in loop:
        outputs["steps"] = int(loop["steps"])
    elif "run_steps" in input_meta:
        outputs["steps"] = int(input_meta["run_steps"])

    if "timesteps_per_second" in loop:
        outputs["timesteps_per_second"] = float(loop["timesteps_per_second"])
    elif "steps" in outputs and "loop_time_seconds" in outputs and float(outputs["loop_time_seconds"]) > 0.0:
        outputs["timesteps_per_second"] = float(outputs["steps"]) / float(outputs["loop_time_seconds"])

    if "ns_per_day" in loop:
        outputs["ns_per_day"] = float(loop["ns_per_day"])
    else:
        timestep_fs = _as_float(input_meta.get("timestep_fs"))
        tps = _as_float(outputs.get("timesteps_per_second"))
        if timestep_fs is not None and tps is not None:
            outputs["ns_per_day"] = tps * timestep_fs * 1.0e-6 * 86400.0

    numeric_vals = [value for value in outputs.values() if isinstance(value, (float, int))]
    outputs["nan_inf_count"] = int(_count_nan_inf(full_text, [float(item) for item in numeric_vals]))
    return outputs


def _run_lammps_once(
    case_id: str,
    repeat_idx: int,
    phase: str,
    timeout_seconds: float | None,
    input_root: Path,
    project_root: Path,
    case: dict[str, Any],
    mpi_launcher: str,
) -> dict[str, Any]:
    input_script_raw = str(case.get("input_script") or "")
    data_file_raw = str(case.get("data_file") or "")
    if not input_script_raw:
        return {
            "success": False,
            "elapsed": 0.0,
            "error": "missing case.input_script",
            "return_code": 2,
            "outputs": {},
            "rss_mb": 0.0,
        }

    script_path = _resolve_case_file(input_script_raw, input_root, project_root)
    if not script_path.exists():
        return {
            "success": False,
            "elapsed": 0.0,
            "error": f"input script not found: {input_script_raw}",
            "return_code": 2,
            "outputs": {},
            "rss_mb": 0.0,
        }

    if data_file_raw:
        data_file_path = _resolve_case_file(data_file_raw, input_root, project_root)
        if not data_file_path.exists():
            return {
                "success": False,
                "elapsed": 0.0,
                "error": f"data file not found: {data_file_raw}",
                "return_code": 2,
                "outputs": {},
                "rss_mb": 0.0,
            }

    input_meta = _parse_input_script_metadata(script_path)

    lmp_exec_raw = str(case.get("lmp_executable") or "build/lmp")
    lmp_exec = Path(lmp_exec_raw).expanduser()
    if not lmp_exec.is_absolute():
        lmp_exec = (project_root / lmp_exec).resolve()

    if not lmp_exec.exists():
        return {
            "success": False,
            "elapsed": 0.0,
            "error": f"LAMMPS executable not found: {lmp_exec}",
            "return_code": 2,
            "outputs": {},
            "rss_mb": 0.0,
        }

    mpi_ranks = _as_int(case.get("mpi_ranks")) or 1
    omp_threads = _as_int(case.get("omp_num_threads")) or 1

    runs_dir = (project_root / ".fermilink-optimize" / "runs" / "autogen").resolve()
    runs_dir.mkdir(parents=True, exist_ok=True)
    log_path = runs_dir / f"{case_id}.{phase}.{repeat_idx}.log"

    script_arg = _command_path(script_path, input_root)
    command: list[str] = [
        mpi_launcher,
        "-n",
        str(mpi_ranks),
        str(lmp_exec),
        "-in",
        script_arg,
        "-log",
        str(log_path),
    ]

    time_prefix, time_mode = _choose_time_wrapper()
    timed_command = [*time_prefix, *command]

    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(omp_threads)
    env["OPENBLAS_NUM_THREADS"] = "1"
    env["MKL_NUM_THREADS"] = "1"
    env["NUMEXPR_NUM_THREADS"] = "1"

    started = time.perf_counter()
    try:
        completed = subprocess.run(
            timed_command,
            cwd=str(input_root),
            text=True,
            capture_output=True,
            timeout=timeout_seconds,
            check=False,
            env=env,
        )
        elapsed = max(0.0, time.perf_counter() - started)
    except subprocess.TimeoutExpired:
        elapsed = max(0.0, time.perf_counter() - started)
        return {
            "success": False,
            "elapsed": elapsed,
            "error": "timeout",
            "return_code": 124,
            "outputs": {},
            "rss_mb": 0.0,
        }
    except OSError as exc:
        elapsed = max(0.0, time.perf_counter() - started)
        return {
            "success": False,
            "elapsed": elapsed,
            "error": str(exc),
            "return_code": 1,
            "outputs": {},
            "rss_mb": 0.0,
        }

    log_text = ""
    if log_path.exists():
        try:
            log_text = log_path.read_text(encoding="utf-8", errors="replace")
        except OSError:
            log_text = ""

    stdout_text = str(completed.stdout or "")
    stderr_text = str(completed.stderr or "")
    combined = "\n".join([stdout_text, stderr_text, log_text])

    outputs = _parse_case_outputs(combined, input_meta)
    rss_mb = _parse_max_rss_mb(stderr_text, time_mode) or 0.0

    success = completed.returncode == 0
    error = ""
    if not success:
        error = stderr_text.strip() or stdout_text.strip() or f"return code {completed.returncode}"

    return {
        "success": success,
        "elapsed": elapsed,
        "error": error,
        "return_code": int(completed.returncode),
        "outputs": outputs,
        "rss_mb": float(rss_mb),
    }


def _weighted_median(values_and_weights: list[tuple[float, float]]) -> float:
    cleaned = [(v, w) for v, w in values_and_weights if math.isfinite(v) and math.isfinite(w) and w > 0.0]
    if not cleaned:
        return 0.0
    cleaned.sort(key=lambda item: item[0])
    total_weight = sum(weight for _, weight in cleaned)
    cutoff = total_weight / 2.0
    cumulative = 0.0
    for value, weight in cleaned:
        cumulative += weight
        if cumulative >= cutoff:
            return float(value)
    return float(cleaned[-1][0])


def _median(values: list[float]) -> float:
    if not values:
        return 0.0
    return float(statistics.median(values))


def _has_tolerance_threshold(item: dict[str, Any]) -> bool:
    for key in ("abs_delta", "relative_delta", "rms_delta"):
        if _as_float(item.get(key)) is not None:
            return True
    legacy_mode = item.get("mode")
    legacy_value = _as_float(item.get("value"))
    if isinstance(legacy_mode, str) and legacy_mode in {"abs_delta", "relative_delta", "rms_delta"}:
        return legacy_value is not None
    return False


def _run_case(
    project_root: Path,
    input_root: Path,
    case: dict[str, Any],
    required_fields: list[str],
    controller_timeout: float | None,
    mpi_launcher: str,
) -> dict[str, Any]:
    case_id = str(case.get("id") or "case")
    warmups = max(0, _as_int(case.get("warmup_repeats")) or 0)
    measured = max(1, _as_int(case.get("repeats")) or 1)
    case_timeout = _as_float(case.get("timeout_seconds"))
    timeout_seconds = case_timeout if case_timeout is not None and case_timeout > 0.0 else controller_timeout
    weight = _as_float(case.get("weight")) or 1.0

    max_rss_mb = 0.0
    measured_records: list[dict[str, Any]] = []
    errors: list[str] = []

    for idx in range(warmups):
        outcome = _run_lammps_once(
            case_id=case_id,
            repeat_idx=idx,
            phase="warmup",
            timeout_seconds=timeout_seconds,
            input_root=input_root,
            project_root=project_root,
            case=case,
            mpi_launcher=mpi_launcher,
        )
        max_rss_mb = max(max_rss_mb, float(outcome.get("rss_mb") or 0.0))
        if not bool(outcome.get("success")):
            errors.append(f"warmup_{idx}: {outcome.get('error') or 'failed'}")

    for idx in range(measured):
        outcome = _run_lammps_once(
            case_id=case_id,
            repeat_idx=idx,
            phase="measured",
            timeout_seconds=timeout_seconds,
            input_root=input_root,
            project_root=project_root,
            case=case,
            mpi_launcher=mpi_launcher,
        )
        max_rss_mb = max(max_rss_mb, float(outcome.get("rss_mb") or 0.0))
        if bool(outcome.get("success")):
            measured_records.append(outcome)
        else:
            errors.append(f"measured_{idx}: {outcome.get('error') or 'failed'}")

    elapsed_values = [float(item.get("elapsed") or 0.0) for item in measured_records]
    wall_seconds = _median(elapsed_values)
    total_seconds = float(sum(elapsed_values)) if elapsed_values else 0.0

    steps_values: list[float] = []
    tps_values: list[float] = []
    ns_day_values: list[float] = []
    loop_values: list[float] = []
    for item in measured_records:
        outputs = item.get("outputs") if isinstance(item.get("outputs"), dict) else {}
        steps = _as_float(outputs.get("steps"))
        if steps is not None:
            steps_values.append(steps)
        tps = _as_float(outputs.get("timesteps_per_second"))
        if tps is not None:
            tps_values.append(tps)
        ns_day = _as_float(outputs.get("ns_per_day"))
        if ns_day is not None:
            ns_day_values.append(ns_day)
        loop_time = _as_float(outputs.get("loop_time_seconds"))
        if loop_time is not None:
            loop_values.append(loop_time)

    steps = int(round(_median(steps_values))) if steps_values else (_as_int(case.get("run_steps")) or 0)
    wall_per_100 = wall_seconds * (100.0 / steps) if steps > 0 else wall_seconds

    representative_outputs: dict[str, Any] = {}
    if measured_records:
        maybe_outputs = measured_records[-1].get("outputs")
        if isinstance(maybe_outputs, dict):
            representative_outputs = dict(maybe_outputs)

    result: dict[str, Any] = {
        "id": case_id,
        "weight": float(weight),
        "converged": len(errors) == 0 and len(measured_records) == measured,
        "wall_seconds": float(wall_seconds),
        "total_seconds": float(total_seconds),
        "loop_time_seconds": float(_median(loop_values)) if loop_values else 0.0,
        "steps": int(steps),
        "wall_seconds_per_100_steps": float(wall_per_100),
        "timesteps_per_second": float(_median(tps_values)) if tps_values else 0.0,
        "ns_per_day": float(_median(ns_day_values)) if ns_day_values else 0.0,
        "peak_rss_mb": float(max_rss_mb),
        "error": "; ".join(errors),
    }

    for key, value in representative_outputs.items():
        result[key] = value

    for field in required_fields:
        if field not in result:
            result["converged"] = False
            if result["error"]:
                result["error"] += "; "
            result["error"] += f"missing output field: {field}"

    nan_inf = _as_int(result.get("nan_inf_count"))
    if nan_inf is not None and nan_inf > 0:
        result["converged"] = False
        if result["error"]:
            result["error"] += "; "
        result["error"] += "nan/inf detected"

    if not result["converged"] and not result["error"]:
        result["error"] = "case failed"

    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--benchmark", required=True)
    parser.add_argument("--emit-json", action="store_true")
    args = parser.parse_args()

    benchmark_path = Path(args.benchmark).expanduser().resolve()
    payload_raw = yaml.safe_load(benchmark_path.read_text(encoding="utf-8"))
    if not isinstance(payload_raw, dict):
        raise SystemExit("benchmark YAML must be a mapping")

    payload = dict(payload_raw)
    cases_raw = payload.get("cases")
    if not isinstance(cases_raw, list) or not cases_raw:
        raise SystemExit("benchmark YAML must define non-empty cases")
    cases: list[dict[str, Any]] = [dict(item) for item in cases_raw if isinstance(item, dict)]
    if not cases:
        raise SystemExit("benchmark YAML must define at least one object case")

    project_root = _find_project_root(benchmark_path.parent)
    input_root = _resolve_input_root(project_root)
    if not input_root.exists() or not input_root.is_dir():
        raise SystemExit(f"input root not found or not a directory: {input_root}")

    mpi_launcher = _choose_mpi_launcher()
    if mpi_launcher is None:
        raise SystemExit("unable to find MPI launcher: expected mpiexec or mpirun in PATH")

    objective = payload.get("controller", {}).get("objective", {}) if isinstance(payload.get("controller"), dict) else {}
    primary_metric = str(objective.get("primary_metric") or "weighted_median_wall_seconds_per_100_steps")

    correctness = payload.get("correctness")
    required_fields: list[str] = []
    if isinstance(correctness, dict):
        field_tolerances = correctness.get("field_tolerances")
        if isinstance(field_tolerances, list):
            for item in field_tolerances:
                if isinstance(item, dict):
                    field = item.get("field")
                    if isinstance(field, str) and field and _has_tolerance_threshold(item):
                        required_fields.append(field)

    controller_timeout = None
    controller = payload.get("controller")
    if isinstance(controller, dict):
        timeout_candidate = _as_float(controller.get("timeout_seconds"))
        if timeout_candidate is not None and timeout_candidate > 0.0:
            controller_timeout = timeout_candidate

    case_results: list[dict[str, Any]] = []
    for case in cases:
        case_results.append(
            _run_case(
                project_root=project_root,
                input_root=input_root,
                case=case,
                required_fields=required_fields,
                controller_timeout=controller_timeout,
                mpi_launcher=mpi_launcher,
            )
        )

    weighted_values: list[tuple[float, float]] = []
    peak_rss_mb = 0.0
    correctness_ok = True

    for case in case_results:
        weight = _as_float(case.get("weight")) or 1.0
        metric_value = _as_float(case.get("wall_seconds_per_100_steps"))
        converged = bool(case.get("converged"))
        nan_inf = _as_int(case.get("nan_inf_count")) or 0

        if metric_value is None:
            metric_value = FAILURE_PENALTY_METRIC
        if not converged:
            metric_value = max(metric_value, FAILURE_PENALTY_METRIC)
        weighted_values.append((float(metric_value), float(weight)))

        peak_rss_mb = max(peak_rss_mb, _as_float(case.get("peak_rss_mb")) or 0.0)
        if not converged or nan_inf != 0:
            correctness_ok = False

    primary_value = _weighted_median(weighted_values)

    output = {
        "benchmark_id": str(payload.get("benchmark_id") or benchmark_path.stem),
        "correctness_ok": bool(correctness_ok),
        "summary_metrics": {
            primary_metric: float(primary_value),
            "peak_rss_mb": float(peak_rss_mb),
        },
        "cases": case_results,
    }

    print(json.dumps(output, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())