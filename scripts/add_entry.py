"""Add or refresh a Project Optimization entry.

Thin wrapper around ``skills/optimize-report/assets/build_report.py`` that
standardizes the output path inside this Sphinx project:

    source/entries/<package>/<task>/

Usage:

    python project-optimization/scripts/add_entry.py <path-to-project>

The script accepts either a project root containing ``.fermilink-optimize/``
or the ``.fermilink-optimize`` directory itself. ``--package`` and ``--task``
remain available as overrides, but when omitted they are inferred from the
target repo's current branch name. The default branch pattern is expected to be:

    fermilink-optimize/<package>-<task>

Optional flags (passed through to the skill): ``--title``, ``--metric-label``,
``--direction {lower,higher}``.
"""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
SKILL_BUILDER = REPO_ROOT / "skills" / "optimize-report" / "assets" / "build_report.py"
PROJECT_ROOT = Path(__file__).resolve().parents[1]
ENTRIES_ROOT = PROJECT_ROOT / "source" / "entries"

SLUG_OK = re.compile(r"^[a-z0-9][a-z0-9\-_]*$")
ENTRY_BRANCH = re.compile(r"^fermilink-optimize/(?P<spec>[a-z0-9][a-z0-9\-_]*)$")
WORKER_BRANCH = re.compile(
    r"^fermilink-optimize-worker/fermilink-optimize-"
    r"(?P<spec>[a-z0-9][a-z0-9\-_]*?)(?:-[0-9a-f]{12})?$"
)


def _validate_slug(name: str, label: str) -> str:
    value = name.strip().lower()
    if not SLUG_OK.match(value):
        raise SystemExit(
            f"invalid --{label} '{name}': use lowercase letters, digits, '-', '_' "
            "(must start with a letter or digit)"
        )
    return value


def _resolve_optimize_dir(path: Path) -> Path:
    candidate = path.resolve()
    if (candidate / "results.tsv").exists():
        return candidate

    nested = candidate / ".fermilink-optimize"
    if (nested / "results.tsv").exists():
        return nested

    raise SystemExit(
        f"no results.tsv under {candidate} or {nested} — pass a project directory "
        "containing .fermilink-optimize/ or the .fermilink-optimize directory itself"
    )


def _get_current_branch(path: Path) -> str | None:
    result = subprocess.run(
        ["git", "-C", str(path), "branch", "--show-current"],
        check=False,
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        return None
    branch = result.stdout.strip()
    return branch or None


def _parse_package_task_from_branch(branch: str) -> tuple[str, str] | None:
    for pattern in (ENTRY_BRANCH, WORKER_BRANCH):
        match = pattern.match(branch.strip())
        if not match:
            continue
        spec = match.group("spec")
        package, sep, task = spec.partition("-")
        if sep and package and task:
            return package, task
    return None


def _infer_package_task(path: Path) -> tuple[str, str, str]:
    branch = _get_current_branch(path)
    if not branch:
        raise SystemExit(
            "could not determine the current git branch from "
            f"{path}; pass --package and --task explicitly"
        )

    parsed = _parse_package_task_from_branch(branch)
    if parsed is None:
        raise SystemExit(
            "could not infer --package/--task from branch "
            f"'{branch}'; expected 'fermilink-optimize/<package>-<task>' "
            "or pass --package and --task explicitly"
        )

    package, task = parsed
    return branch, package, task


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "optimize_dir",
        type=Path,
        help="path to a project directory or completed .fermilink-optimize directory",
    )
    ap.add_argument(
        "--package",
        default=None,
        help="package slug (e.g. pyscf, lammps); defaults to the current branch",
    )
    ap.add_argument(
        "--task",
        default=None,
        help="task slug (e.g. diis-scf); defaults to the current branch",
    )
    ap.add_argument("--title", default=None, help="report title (optional)")
    ap.add_argument("--metric-label", default=None)
    ap.add_argument("--direction", choices=["lower", "higher"], default=None)
    args = ap.parse_args()

    if not SKILL_BUILDER.exists():
        raise SystemExit(f"builder not found: {SKILL_BUILDER}")
    optimize_dir = _resolve_optimize_dir(args.optimize_dir)

    inferred_branch = None
    inferred_package = None
    inferred_task = None
    if args.package is None or args.task is None:
        inferred_branch, inferred_package, inferred_task = _infer_package_task(optimize_dir)

    package = _validate_slug(args.package or inferred_package, "package")
    task = _validate_slug(args.task or inferred_task, "task")
    out_dir = ENTRIES_ROOT / package / task

    cmd: list[str] = [
        sys.executable,
        str(SKILL_BUILDER),
        str(optimize_dir),
        "--out",
        str(out_dir),
    ]
    if args.title:
        cmd += ["--title", args.title]
    if args.metric_label:
        cmd += ["--metric-label", args.metric_label]
    if args.direction:
        cmd += ["--direction", args.direction]

    print(f"[add_entry] invoking {SKILL_BUILDER.name}")
    print(f"[add_entry] source : {optimize_dir}")
    if inferred_branch:
        print(f"[add_entry] branch : {inferred_branch}")
        print(f"[add_entry] entry  : {package}/{task}")
    print(f"[add_entry] target : {out_dir}")
    result = subprocess.run(cmd, check=False)
    if result.returncode != 0:
        raise SystemExit(result.returncode)

    print(f"[add_entry] done. Rebuild the site with: cd {PROJECT_ROOT.name} && make doc")


if __name__ == "__main__":
    main()
