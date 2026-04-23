"""Sphinx config for the FermiLink Project Optimization gallery.

This site aggregates per-task report bundles produced by
``skills/optimize-report/assets/build_report.py`` (driven via
``project-optimization/scripts/add_entry.py``). Entries live under
``source/entries/<package>/<task>/`` and are treated as immutable inputs.

On every build, ``_regenerate_package_pages`` scans ``source/entries/`` for
``data/summary.json`` manifests and writes one RST page per package under
``source/packages/<package>.rst`` with a summary table + toctree pointing at
the task bundles. ``source/packages/`` is gitignored because it is fully
derived from the entries.
"""

from __future__ import annotations

import json
from datetime import datetime
from pathlib import Path

SOURCE = Path(__file__).resolve().parent
ROOT = SOURCE.parent
ENTRIES_DIR = SOURCE / "entries"
PACKAGES_DIR = SOURCE / "packages"

# --- Project info ---
project = "FermiLink Project Optimization"
author = "Tao E. Li"
copyright = "TEL Research Group 2026"
release = version = datetime.now().strftime("%Y.%m.%d")
html_title = "FermiLink · Project Optimization"
html_short_title = "Project Optimization"

# --- Extensions ---
extensions = [
    "sphinx.ext.mathjax",
    "myst_parser",
]

# --- Theme ---
html_theme = "furo"
html_theme_options = {
    "sidebar_hide_name": False,
    "navigation_with_keys": False,
    "light_logo": "img/mark-light.svg",
    "dark_logo": "img/mark.svg",
    "light_css_variables": {
        "color-brand-primary": "#1264a3",
        "color-brand-content": "#0d2a4d",
        "color-sidebar-background": "#f6f9ff",
        "color-admonition-background": "rgba(18, 100, 163, 0.08)",
    },
    "dark_css_variables": {
        "color-brand-primary": "#66c7ff",
        "color-brand-content": "#d6ecff",
        "color-sidebar-background": "#0d1829",
        "color-admonition-background": "rgba(102, 199, 255, 0.12)",
    },
}

templates_path = ["_templates"]
exclude_patterns = [
    "_build",
    # Contract sidecars are reference text inside each entry, not site pages.
    "entries/**/contract/**",
]
html_static_path = ["_static"]
html_css_files = ["css/custom.css"]

html_meta = {"referrer": "strict-origin-when-cross-origin"}


# ---------------------------------------------------------------------------
# Per-package page regeneration
# ---------------------------------------------------------------------------


def _discover_entries() -> dict[str, list[dict]]:
    """Return {package_id: [entry_record, ...]} sorted by task name."""
    if not ENTRIES_DIR.exists():
        return {}
    by_package: dict[str, list[dict]] = {}
    for package_dir in sorted(p for p in ENTRIES_DIR.iterdir() if p.is_dir()):
        tasks: list[dict] = []
        for task_dir in sorted(t for t in package_dir.iterdir() if t.is_dir()):
            index_rst = task_dir / "index.rst"
            summary_json = task_dir / "data" / "summary.json"
            if not index_rst.exists():
                continue
            record: dict = {
                "task_id": task_dir.name,
                "doc_ref": f"../entries/{package_dir.name}/{task_dir.name}/index",
            }
            if summary_json.exists():
                try:
                    record["summary"] = json.loads(summary_json.read_text())
                except json.JSONDecodeError:
                    record["summary"] = {}
            else:
                record["summary"] = {}
            tasks.append(record)
        if tasks:
            by_package[package_dir.name] = tasks
    return by_package


def _format_metric(value) -> str:
    if value is None:
        return "—"
    try:
        return f"{float(value):.6g}"
    except (TypeError, ValueError):
        return str(value)


def _render_package_page(package_id: str, tasks: list[dict]) -> str:
    title = package_id
    pretty = package_id.replace("-", " ").replace("_", " ").title()
    header = (
        f"{pretty} — optimization runs\n"
        f"{'=' * (len(pretty) + len(' — optimization runs'))}\n"
    )

    lines: list[str] = [header, ""]
    lines.append(
        f"{len(tasks)} optimize run{'s' if len(tasks) != 1 else ''} "
        f"recorded for ``{package_id}``."
    )
    lines.append("")

    lines.append(".. list-table:: Runs")
    lines.append("   :header-rows: 1")
    lines.append("   :widths: 28 18 18 12 12 12")
    lines.append("")
    lines.append("   * - Task")
    lines.append("     - Metric")
    lines.append("     - Direction")
    lines.append("     - Baseline")
    lines.append("     - Best")
    lines.append("     - Δ vs baseline")
    for t in tasks:
        s = t.get("summary", {}) or {}
        baseline = (s.get("baseline") or {}).get("metric")
        best = (s.get("best") or {}).get("metric")
        pct = (s.get("best") or {}).get("pct_vs_baseline")
        metric_label = s.get("metric_label") or "—"
        direction = s.get("direction") or "—"
        pct_str = "—" if pct is None else f"{float(pct):+.2f}%"
        lines.append(f"   * - :doc:`{t['task_id'] } <{t['doc_ref']}>`")
        lines.append(f"     - ``{metric_label}``")
        lines.append(f"     - {direction}")
        lines.append(f"     - {_format_metric(baseline)}")
        lines.append(f"     - {_format_metric(best)}")
        lines.append(f"     - {pct_str}")
    lines.append("")

    lines.append(".. toctree::")
    lines.append("   :maxdepth: 1")
    lines.append("   :hidden:")
    lines.append("")
    for t in tasks:
        lines.append(f"   {t['doc_ref']}")
    lines.append("")

    return "\n".join(lines)


def _regenerate_package_pages() -> None:
    by_package = _discover_entries()
    PACKAGES_DIR.mkdir(parents=True, exist_ok=True)
    # Wipe stale generated pages (they are gitignored).
    for existing in PACKAGES_DIR.glob("*.rst"):
        existing.unlink()
    for package_id, tasks in by_package.items():
        out = PACKAGES_DIR / f"{package_id}.rst"
        out.write_text(_render_package_page(package_id, tasks), encoding="utf-8")


def _build_landing_packages() -> list[dict]:
    """Summarize discovered packages for the landing-page card grid."""
    cards: list[dict] = []
    for package_id, tasks in _discover_entries().items():
        pretty = package_id.replace("-", " ").replace("_", " ").title()
        best_delta_pct = None
        for t in tasks:
            pct = ((t.get("summary") or {}).get("best") or {}).get("pct_vs_baseline")
            if pct is None:
                continue
            try:
                pct_f = float(pct)
            except (TypeError, ValueError):
                continue
            # pct_vs_baseline is already oriented as "improvement"
            # (positive = better) by optimize-report, so pick the max.
            if best_delta_pct is None or pct_f > best_delta_pct:
                best_delta_pct = pct_f
        cards.append({
            "id": package_id,
            "title": pretty,
            "docname": f"packages/{package_id}",
            "runs_count": len(tasks),
            "best_delta_pct": (
                None if best_delta_pct is None else f"{best_delta_pct:+.1f}%"
            ),
        })
    cards.sort(key=lambda c: c["id"])
    return cards


html_context = {"po_packages": _build_landing_packages()}


def _on_builder_inited(app: object) -> None:
    del app
    _regenerate_package_pages()


def setup(app: object) -> dict[str, bool]:
    app.connect("builder-inited", _on_builder_inited)
    return {"parallel_read_safe": True, "parallel_write_safe": True}
