"""Add or refresh a Project Optimization entry.

Thin wrapper around ``skills/optimize-report/assets/build_report.py`` that
standardizes the output path inside this Sphinx project:

    source/entries/<package>/<task>/

Usage:

    python project-optimization/scripts/add_entry.py \\
        <path-to>/.fermilink-optimize \\
        --package pyscf --task diis-scf

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


def _validate_slug(name: str, label: str) -> str:
    value = name.strip().lower()
    if not SLUG_OK.match(value):
        raise SystemExit(
            f"invalid --{label} '{name}': use lowercase letters, digits, '-', '_' "
            "(must start with a letter or digit)"
        )
    return value


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "optimize_dir",
        type=Path,
        help="path to a completed .fermilink-optimize directory",
    )
    ap.add_argument("--package", required=True, help="package slug (e.g. pyscf, lammps)")
    ap.add_argument("--task", required=True, help="task slug (e.g. diis-scf)")
    ap.add_argument("--title", default=None, help="report title (optional)")
    ap.add_argument("--metric-label", default=None)
    ap.add_argument("--direction", choices=["lower", "higher"], default=None)
    args = ap.parse_args()

    if not SKILL_BUILDER.exists():
        raise SystemExit(f"builder not found: {SKILL_BUILDER}")
    optimize_dir = args.optimize_dir.resolve()
    if not (optimize_dir / "results.tsv").exists():
        raise SystemExit(
            f"no results.tsv under {optimize_dir} — not a valid "
            ".fermilink-optimize workspace"
        )

    package = _validate_slug(args.package, "package")
    task = _validate_slug(args.task, "task")
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
    print(f"[add_entry] target : {out_dir}")
    result = subprocess.run(cmd, check=False)
    if result.returncode != 0:
        raise SystemExit(result.returncode)

    print(f"[add_entry] done. Rebuild the site with: cd {PROJECT_ROOT.name} && make doc")


if __name__ == "__main__":
    main()
