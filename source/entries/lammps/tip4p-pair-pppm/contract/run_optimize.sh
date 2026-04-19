#!/usr/bin/env bash
set -euo pipefail

fermilink optimize \
  lammps \
  /Users/taoli/.fermilink/workspaces/test_optimize/lammps-optimize-tip4p \
  --benchmark /Users/taoli/.fermilink/workspaces/test_optimize/lammps-optimize-tip4p/.fermilink-optimize/autogen/benchmark.yaml \
  --skills-source existing \
  --max-iterations 30 \
  --stop-on-consecutive-rejections 8 \
  --worker-max-iterations 8 \
  --worker-wait-seconds 1 \
  --worker-max-wait-seconds 900 \
  --worker-pid-stall-seconds 300 \
  "$@"
