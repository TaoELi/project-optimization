Optimization Report — lammps-neighbor
=====================================


Primary metric: ``Weighted median neigh time (s)`` (lower is better).

Summary
-------


- baseline (`e7c0ed95a333 <summary-baseline-e7c0ed95a333_>`_): ``0.36158``
- best accepted (`11db74893452 <summary-best-11db74893452_>`_): ``0.28487`` (+21.22% vs baseline)
- published GitHub branch: `fermilink-optimize/lammps-neighbor <https://github.com/skilled-scipkg/lammps/tree/fermilink-optimize%2Flammps-neighbor>`_
- iterations: 28 total | 6 accepted | 20 rejected | 0 correctness failure

Optimization Trajectory
-----------------------


.. image:: img/metric_vs_iter.svg
   :width: 100%
   :alt: metric vs iteration

.. image:: img/improvement_cumulative.svg
   :width: 100%
   :alt: running incumbent

All iterations
--------------


+------+-------------------------------------------------+-------------------+---------+--------------------------------------------------------------------------------------------------------+
| iter | commit                                          | status            | metric  | summary                                                                                                |
+======+=================================================+===================+=========+========================================================================================================+
| 0    | `e7c0ed95a333 <iter-0000-table-e7c0ed95a333_>`_ | baseline          | 0.36158 | baseline                                                                                               |
+------+-------------------------------------------------+-------------------+---------+--------------------------------------------------------------------------------------------------------+
| 1    | `2d1eda35b944 <iter-0001-table-2d1eda35b944_>`_ | rejected          | 0.47712 | Split the orthogonal half-bin/newton self-bin traversal in \`npair_bin.cpp\` and cache per-atom cut…   |
+------+-------------------------------------------------+-------------------+---------+--------------------------------------------------------------------------------------------------------+
| 2    | `32b495fdb66b <iter-0002-table-32b495fdb66b_>`_ | accepted          | 0.33535 | Specialize the orthogonal half/bin/newton neighbor build in \`src/npair_bin.cpp\` to remove generic…   |
+------+-------------------------------------------------+-------------------+---------+--------------------------------------------------------------------------------------------------------+
| 3    | `380bf2b0146b <iter-0003-table-380bf2b0146b_>`_ | accepted          | 0.32842 | Add a benchmark-specific fast path in \`NPairBin<1,1,0,0,0>::build()\` for standard molecular syste…   |
+------+-------------------------------------------------+-------------------+---------+--------------------------------------------------------------------------------------------------------+
| 4    | `3bd7cbc7ee03 <iter-0004-table-3bd7cbc7ee03_>`_ | rejected          | 0.32287 | Skip \`find_special()\` and \`minimum_image_check()\` for different-molecule pairs inside the incumbe… |
+------+-------------------------------------------------+-------------------+---------+--------------------------------------------------------------------------------------------------------+
| 5    | `1a1935b59c00 <iter-0005-table-1a1935b59c00_>`_ | rejected          | 0.32384 | Add a \`maxspecial <= 2\` molecular/no-exclusion fast path in \`NPairBin<1,1,0,0,0>::build()\` that b… |
+------+-------------------------------------------------+-------------------+---------+--------------------------------------------------------------------------------------------------------+
| 6    | `c4128d3970b2 <iter-0006-table-c4128d3970b2_>`_ | accepted          | 0.31969 | Refine the incumbent \`NPairBin<1,1,0,0,0>::build()\` fast path by early-accepting different-molecu…   |
+------+-------------------------------------------------+-------------------+---------+--------------------------------------------------------------------------------------------------------+
| 7    | `d89e6f820474 <iter-0007-table-d89e6f820474_>`_ | rejected          | 0.32426 | Replace the remaining same-molecule \`find_special()\` scan in the incumbent \`NPairBin<1,1,0,0,0>::…  |
+------+-------------------------------------------------+-------------------+---------+--------------------------------------------------------------------------------------------------------+
| 8    | `d381c38751ec <iter-0008-table-d381c38751ec_>`_ | rejected          | 0.33329 | Split the benchmark-hot molecular/no-exclusion neighbor build into encoded-special and generic va…     |
+------+-------------------------------------------------+-------------------+---------+--------------------------------------------------------------------------------------------------------+
| 9    | `2894990122f9 <iter-0009-table-2894990122f9_>`_ | rejected          | 0.3187  | Cache per-neighbor ownership/tag/molecule values and use an encoded-only special lookup inside th…     |
+------+-------------------------------------------------+-------------------+---------+--------------------------------------------------------------------------------------------------------+
| 10   | `b21ee6319e89 <iter-0010-table-b21ee6319e89_>`_ | rejected          | 0.3221  | Split the incumbent molecular/no-exclusion self-bin walk in \`NPairBin<1,1,0,0,0>::build()\` into o…   |
+------+-------------------------------------------------+-------------------+---------+--------------------------------------------------------------------------------------------------------+
| 11   | `8c29b124a2c6 <iter-0011-table-8c29b124a2c6_>`_ | accepted          | 0.3115  | Use contiguous \`_noalias\` aliases for hot coordinate/type/tag/molecule arrays in the incumbent \`N…  |
+------+-------------------------------------------------+-------------------+---------+--------------------------------------------------------------------------------------------------------+
| 12   | `e5d6d4548fe4 <iter-0012-table-e5d6d4548fe4_>`_ | rejected          | 0.31003 | Extend the incumbent \`NPairBin<1,1,0,0,0>::build()\` hot path with additional neighbor-array alias…   |
+------+-------------------------------------------------+-------------------+---------+--------------------------------------------------------------------------------------------------------+
| 13   | `7adf3bfc72a8 <iter-0013-table-7adf3bfc72a8_>`_ | rejected          | 0.31711 | Replace the incumbent molecular/no-exclusion \`NPairBin<1,1,0,0,0>::build()\` special handling with…   |
+------+-------------------------------------------------+-------------------+---------+--------------------------------------------------------------------------------------------------------+
| 14   | `e536920bfbe5 <iter-0014-table-e536920bfbe5_>`_ | rejected          | 0.30666 | Add explicit \`jnext\` traversal and light next-neighbor coordinate prefetching to the incumbent \`N…  |
+------+-------------------------------------------------+-------------------+---------+--------------------------------------------------------------------------------------------------------+
| 15   | `e5eef089d2d0 <iter-0015-table-e5eef089d2d0_>`_ | rejected          | 0.30618 | Combine explicit jnext traversal and next-neighbor coordinate prefetching with cached periodic ha…     |
+------+-------------------------------------------------+-------------------+---------+--------------------------------------------------------------------------------------------------------+
| 16   | `d5b0238dc967 <iter-0016-table-d5b0238dc967_>`_ | rejected          | 0.30551 | Refine the single-path molecular half/bin/newton builder in src/npair_bin.cpp with explicit jnext…     |
+------+-------------------------------------------------+-------------------+---------+--------------------------------------------------------------------------------------------------------+
| 17   | `4a9dda87cf45 <iter-0017-table-4a9dda87cf45_>`_ | rejected          | 0.30695 | Combine iteration-16 style \`NPairBin<1,1,0,0,0>\` jnext/prefetch/half-box fast-path refinements wi…   |
+------+-------------------------------------------------+-------------------+---------+--------------------------------------------------------------------------------------------------------+
| 18   | `08f5ab4a4e02 <iter-0018-table-08f5ab4a4e02_>`_ | accepted          | 0.29094 | Use explicit \`jnext\` linked-list traversal with bin/list \`_noalias\` aliases and delayed generic-o… |
+------+-------------------------------------------------+-------------------+---------+--------------------------------------------------------------------------------------------------------+
| 19   | `11db74893452 <iter-0019-table-11db74893452_>`_ | accepted          | 0.28487 | Optimize Neighbor movement tracking by scanning contiguous x/xhold arrays in check_distance() and…     |
+------+-------------------------------------------------+-------------------+---------+--------------------------------------------------------------------------------------------------------+
| 20   | `11db74893452 <iter-0020-table-11db74893452_>`_ | worker_incomplete | nan     | optimize iteration 20                                                                                  |
+------+-------------------------------------------------+-------------------+---------+--------------------------------------------------------------------------------------------------------+
| 21   | `e35c58b86e07 <iter-0021-table-e35c58b86e07_>`_ | rejected          | 0.28995 | Widen Neighbor::check_distance() to 8-atom unrolled blocks and compare one blockwise max displace…     |
+------+-------------------------------------------------+-------------------+---------+--------------------------------------------------------------------------------------------------------+
| 22   | `708f9b2d16b1 <iter-0022-table-708f9b2d16b1_>`_ | rejected          | 0.28503 | Specialize the incumbent \`NPairBin<1,1,0,0,0>::build()\` fast path for all-zero special-bond weigh…   |
+------+-------------------------------------------------+-------------------+---------+--------------------------------------------------------------------------------------------------------+
| 23   | `a3f1052997d8 <iter-0023-table-a3f1052997d8_>`_ | rejected          | 0.28469 | Add a runtime uniform-cutoff fast path to the specialized \`NPairBin<1,1,0,0,0>\` molecular/no-excl…   |
+------+-------------------------------------------------+-------------------+---------+--------------------------------------------------------------------------------------------------------+
| 24   | `84babbbb36e3 <iter-0024-table-84babbbb36e3_>`_ | rejected          | 0.49372 | Add an axis-aligned early reject in the specialized \`NPairBin<1,1,0,0,0>\` fast path so pairs exce…   |
+------+-------------------------------------------------+-------------------+---------+--------------------------------------------------------------------------------------------------------+
| 25   | `b78acabc90a5 <iter-0025-table-b78acabc90a5_>`_ | rejected          | 0.28592 | Isolate the iteration-23 Neighbor cadence fast-check by caching the benchmark-common every=1/dela…     |
+------+-------------------------------------------------+-------------------+---------+--------------------------------------------------------------------------------------------------------+
| 26   | `ba18979d7807 <iter-0026-table-ba18979d7807_>`_ | rejected          | 0.28444 | Add a benchmark-common \`NPairBin<1,1,0,0,0>\` fast path for the uniform-cutoff plus all-special-en…   |
+------+-------------------------------------------------+-------------------+---------+--------------------------------------------------------------------------------------------------------+
| 27   | `321fb756a836 <iter-0027-table-321fb756a836_>`_ | rejected          | 0.28559 | Cache j-side neighbor coords/tag/molecule once and inline orthogonal half-box minimum-image check…     |
+------+-------------------------------------------------+-------------------+---------+--------------------------------------------------------------------------------------------------------+

Accepted Commits
----------------


Accepted candidate detail pages and current manual-review status:

+-----------------------------------------------------+----------------------------------------+
| accepted commit                                     | Human verification                     |
+=====================================================+========================================+
| :doc:`32b495fdb66b <iterations/iter_0002_accepted>` | not verified                           |
+-----------------------------------------------------+----------------------------------------+
| :doc:`380bf2b0146b <iterations/iter_0003_accepted>` | not verified                           |
+-----------------------------------------------------+----------------------------------------+
| :doc:`c4128d3970b2 <iterations/iter_0006_accepted>` | not verified                           |
+-----------------------------------------------------+----------------------------------------+
| :doc:`8c29b124a2c6 <iterations/iter_0011_accepted>` | not verified                           |
+-----------------------------------------------------+----------------------------------------+
| :doc:`08f5ab4a4e02 <iterations/iter_0018_accepted>` | not verified                           |
+-----------------------------------------------------+----------------------------------------+
| :doc:`11db74893452 <iterations/iter_0019_accepted>` | not verified                           |
+-----------------------------------------------------+----------------------------------------+

.. toctree::
   :maxdepth: 1
   :hidden:

   iterations/iter_0002_accepted
   iterations/iter_0003_accepted
   iterations/iter_0006_accepted
   iterations/iter_0011_accepted
   iterations/iter_0018_accepted
   iterations/iter_0019_accepted

Benchmark Contracts
-------------------


Benchmark contract and runner used for this optimization:

- :download:`benchmark.yaml <contract/benchmark.yaml>`
- :download:`benchmark_runner.py <contract/benchmark_runner.py>`
- :download:`goal.md <contract/goal.md>`

Input files for Benchmarks
--------------------------


Copied auxiliary benchmark inputs from ``.fermilink-optimize/inputs/all/``:

- :download:`in.tip4p_nve <inputs/all/in.tip4p_nve>`
- :download:`in.tip4p_nve_long <inputs/all/in.tip4p_nve_long>`
- :download:`water_216_data.lmp <inputs/all/water_216_data.lmp>`

Runtime Data
------------


- :download:`results.tsv <data/results.tsv>`
- :download:`summary.json <data/summary.json>`

Rerun Guide
-----------


Agent provider ``codex``; model ``gpt-5.4-xhigh``

Use the bundled contract files from this report to recreate the optimization against a fresh upstream checkout.

- default upstream clone: ``git@github.com:skilled-scipkg/lammps.git``
- confirm the upstream default branch before creating the worktree: `develop on GitHub <https://github.com/skilled-scipkg/lammps/tree/develop>`_
- detected package language: ``cpp``; use ``fermilink-optimize-cpp`` for goal-mode reruns
- ``contract/run_optimize.sh`` and ``contract/setup_env.sh`` record the original campaign, but they can contain site-specific absolute paths
- if :download:`goal_inputs.json <contract/goal_inputs.json>` is present, restage the listed auxiliary workload files before rerunning
- copied benchmark input files are bundled under ``inputs/all/`` and should be restored into ``.fermilink-optimize/inputs/all/`` for deterministic reruns

.. code-block:: bash

   git clone git@github.com:skilled-scipkg/lammps.git
   cd lammps
   git worktree add -b fermilink-optimize/lammps-<modified-feature> ../lammps-<modified-feature> develop

Path 1: rerun from the bundled :download:`goal.md <contract/goal.md>`.

Run this from the cloned main repo so the launcher can create or reuse the sibling worktree:

.. code-block:: bash

   fermilink-optimize-cpp \
     --project-root "$PWD" \
     --goal /path/to/report/contract/goal.md \
     --branch fermilink-optimize/lammps-<modified-feature> \
     --worktree-root .. \
     --worktree-name lammps-<modified-feature>

Path 2: rerun more deterministically from the copied :download:`benchmark.yaml <contract/benchmark.yaml>` and :download:`benchmark_runner.py <contract/benchmark_runner.py>`.

This avoids regenerating the benchmark contract from ``goal.md`` before the campaign starts:

.. code-block:: bash

   cd ../lammps-<modified-feature>
   mkdir -p .fermilink-optimize/autogen .fermilink-optimize/inputs/all
   cp /path/to/report/contract/benchmark.yaml .fermilink-optimize/autogen/benchmark.yaml
   cp /path/to/report/contract/benchmark_runner.py .fermilink-optimize/autogen/benchmark_runner.py
   cp -R /path/to/report/inputs/all/. .fermilink-optimize/inputs/all/
   printf '%s\n' '.fermilink-optimize/' >> .git/info/exclude
   fermilink optimize lammps "$PWD" \
     --benchmark "$PWD/.fermilink-optimize/autogen/benchmark.yaml" \
     --skills-source existing

Building environment
~~~~~~~~~~~~~~~~~~~~


These commands come from the copied ``## Build`` block in :download:`goal.md <contract/goal.md>` and are rerun before benchmarks through the benchmark configuration.

Check that they work in your local environment before launching a long run. If they do not, update the ``## Build`` section in ``goal.md`` or the corresponding ``runtime.pre_commands`` setting in :download:`benchmark.yaml <contract/benchmark.yaml>`.

.. code-block:: bash

   mkdir -p build
   cd build
   cmake -C ../cmake/presets/most.cmake -C ../cmake/presets/nolib.cmake -D PKG_GPU=off ../cmake
   cmake --build . -j4

Benchmark Examples
------------------


Worker iterations run the ``train-*`` benchmark cases below while searching for candidate changes:

.. code-block:: yaml

   cases:
     - id: train-16r-long
       weight: 1.0
       description: Longer-run 16-rank training case with repeated rebuilds and lower communication noise.
       input_script: in.tip4p_nve_long
       data_file: water_216_data.lmp
       mpi_ranks: 16
       omp_threads: 1
     - id: train-32r-long
       weight: 1.0
       description: Longer-run 32-rank training case for a second domain decomposition.
       input_script: in.tip4p_nve_long
       data_file: water_216_data.lmp
       mpi_ranks: 32
       omp_threads: 1
     - id: train-32r-short
       weight: 1.0
       description: Shorter-turnaround 32-rank training case with the same neighbor settings.
       input_script: in.tip4p_nve
       data_file: water_216_data.lmp
       mpi_ranks: 32
       omp_threads: 1

Controller reviews run the ``test-*`` benchmark cases below to validate accepted candidates:

.. code-block:: yaml

   cases:
     - id: test-16r-short
       weight: 1.0
       description: Held-out lower-rank validation case.
       input_script: in.tip4p_nve
       data_file: water_216_data.lmp
       mpi_ranks: 16
       omp_threads: 1
     - id: test-64r-short
       weight: 1.0
       description: Held-out scaling-sensitive validation case.
       input_script: in.tip4p_nve
       data_file: water_216_data.lmp
       mpi_ranks: 64
       omp_threads: 1


.. _summary-baseline-e7c0ed95a333: https://github.com/skilled-scipkg/lammps/commit/e7c0ed95a333
.. _summary-best-11db74893452: https://github.com/skilled-scipkg/lammps/commit/11db74893452
.. _iter-0000-table-e7c0ed95a333: https://github.com/skilled-scipkg/lammps/commit/e7c0ed95a333
.. _iter-0001-table-2d1eda35b944: https://github.com/skilled-scipkg/lammps/commit/2d1eda35b944
.. _iter-0002-table-32b495fdb66b: https://github.com/skilled-scipkg/lammps/commit/32b495fdb66b
.. _iter-0003-table-380bf2b0146b: https://github.com/skilled-scipkg/lammps/commit/380bf2b0146b
.. _iter-0004-table-3bd7cbc7ee03: https://github.com/skilled-scipkg/lammps/commit/3bd7cbc7ee03
.. _iter-0005-table-1a1935b59c00: https://github.com/skilled-scipkg/lammps/commit/1a1935b59c00
.. _iter-0006-table-c4128d3970b2: https://github.com/skilled-scipkg/lammps/commit/c4128d3970b2
.. _iter-0007-table-d89e6f820474: https://github.com/skilled-scipkg/lammps/commit/d89e6f820474
.. _iter-0008-table-d381c38751ec: https://github.com/skilled-scipkg/lammps/commit/d381c38751ec
.. _iter-0009-table-2894990122f9: https://github.com/skilled-scipkg/lammps/commit/2894990122f9
.. _iter-0010-table-b21ee6319e89: https://github.com/skilled-scipkg/lammps/commit/b21ee6319e89
.. _iter-0011-table-8c29b124a2c6: https://github.com/skilled-scipkg/lammps/commit/8c29b124a2c6
.. _iter-0012-table-e5d6d4548fe4: https://github.com/skilled-scipkg/lammps/commit/e5d6d4548fe4
.. _iter-0013-table-7adf3bfc72a8: https://github.com/skilled-scipkg/lammps/commit/7adf3bfc72a8
.. _iter-0014-table-e536920bfbe5: https://github.com/skilled-scipkg/lammps/commit/e536920bfbe5
.. _iter-0015-table-e5eef089d2d0: https://github.com/skilled-scipkg/lammps/commit/e5eef089d2d0
.. _iter-0016-table-d5b0238dc967: https://github.com/skilled-scipkg/lammps/commit/d5b0238dc967
.. _iter-0017-table-4a9dda87cf45: https://github.com/skilled-scipkg/lammps/commit/4a9dda87cf45
.. _iter-0018-table-08f5ab4a4e02: https://github.com/skilled-scipkg/lammps/commit/08f5ab4a4e02
.. _iter-0019-table-11db74893452: https://github.com/skilled-scipkg/lammps/commit/11db74893452
.. _iter-0020-table-11db74893452: https://github.com/skilled-scipkg/lammps/commit/11db74893452
.. _iter-0021-table-e35c58b86e07: https://github.com/skilled-scipkg/lammps/commit/e35c58b86e07
.. _iter-0022-table-708f9b2d16b1: https://github.com/skilled-scipkg/lammps/commit/708f9b2d16b1
.. _iter-0023-table-a3f1052997d8: https://github.com/skilled-scipkg/lammps/commit/a3f1052997d8
.. _iter-0024-table-84babbbb36e3: https://github.com/skilled-scipkg/lammps/commit/84babbbb36e3
.. _iter-0025-table-b78acabc90a5: https://github.com/skilled-scipkg/lammps/commit/b78acabc90a5
.. _iter-0026-table-ba18979d7807: https://github.com/skilled-scipkg/lammps/commit/ba18979d7807
.. _iter-0027-table-321fb756a836: https://github.com/skilled-scipkg/lammps/commit/321fb756a836