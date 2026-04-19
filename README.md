# FermiLink Project Optimization

Aggregated gallery of `fermilink optimize` runs across scientific packages.
Each entry is a self-contained Sphinx bundle produced by
`skills/optimize-report/assets/build_report.py` from a completed
`.fermilink-optimize` workspace.

The built site is intended to be served at
`https://project-optimization.fermilink.org` and linked from the main
FermiLink documentation hero page.

## Add a new entry

```bash
python project-optimization/scripts/add_entry.py <path-to-project>
```

The script accepts either the project root or the `.fermilink-optimize`
directory itself. By default it infers `<package>` and `<task>` from the
target repo's current branch, which is expected to look like
`fermilink-optimize/<package>-<task>`. You can still override either value with
`--package` and `--task` when needed.

This writes a per-task bundle to
`project-optimization/source/entries/<package>/<task>/`. Per-package summary
pages under `source/packages/` are regenerated automatically on every build.

## Build the site

```bash
cd project-optimization
make doc     # build HTML into build/html
make html    # open the built index in a browser
```

## Layout

```
project-optimization/
  Makefile
  scripts/add_entry.py            # wraps optimize-report skill into this site
  source/
    conf.py                       # Furo theme + auto-aggregation hook
    index.rst                     # landing page (toctree globs over packages/)
    _static/                      # copied from docs/source/_static
    packages/<pkg>.rst            # auto-generated per-package summary (gitignored)
    entries/<pkg>/<task>/         # drop target for bundles (checked in)
  build/html                      # build output (gitignored)
```

`source/packages/` is regenerated on every build from
`source/entries/*/*/data/summary.json`, so it is gitignored. To add a new
package grouping, just drop an entry under `source/entries/<new-pkg>/`.
# project-optimization
