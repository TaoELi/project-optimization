#!/usr/bin/env bash
set -euo pipefail
cd /anvil/scratch/x-tli22/fermilink_optimize/project_pyscf/pyscf-davidson
if [ ! -d .venv ]; then
  python -m venv .venv
fi
. .venv/bin/activate
python -m pip install -U pip
python -m pip install -e .
