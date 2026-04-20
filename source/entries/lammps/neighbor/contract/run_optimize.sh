#!/usr/bin/env bash
set -euo pipefail

fermilink optimize \
  lammps \
  /anvil/scratch/x-tli22/fermilink_optimize/project_lammps/lammps-neighbor \
  --benchmark /anvil/scratch/x-tli22/fermilink_optimize/project_lammps/lammps-neighbor/.fermilink-optimize/autogen/benchmark.yaml \
  --skills-source existing \
  --max-iterations 30 \
  --stop-on-consecutive-rejections 8 \
  --worker-max-iterations 8 \
  --worker-wait-seconds 1 \
  --worker-max-wait-seconds 900 \
  --worker-pid-stall-seconds 300 \
  "$@"
