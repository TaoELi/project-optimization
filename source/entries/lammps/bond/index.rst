Optimization Report — lammps-bond
=================================


Primary metric: ``Weighted median bond time (s)`` (lower is better).

Summary
-------


- baseline (`e7c0ed95a333 <summary-baseline-e7c0ed95a333_>`_): ``0.32627``
- best accepted (`bf07f28c4118 <summary-best-bf07f28c4118_>`_): ``0.21875`` (+32.95% vs baseline)
- published GitHub branch: `fermilink-optimize/lammps-bond <https://github.com/skilled-scipkg/lammps/tree/fermilink-optimize%2Flammps-bond>`_
- iterations: 27 total | 7 accepted | 16 rejected | 3 correctness failure

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


+------+-------------------------------------------------+---------------------+---------+-------------------------------------------------------------------------------------------------------------+
| iter | commit                                          | status              | metric  | summary                                                                                                     |
+======+=================================================+=====================+=========+=============================================================================================================+
| 0    | `e7c0ed95a333 <iter-0000-table-e7c0ed95a333_>`_ | baseline            | 0.32627 | baseline                                                                                                    |
+------+-------------------------------------------------+---------------------+---------+-------------------------------------------------------------------------------------------------------------+
| 1    | `911f06a50100 <iter-0001-table-911f06a50100_>`_ | accepted            | 0.28156 | Hoist \`evflag\`/\`eflag\`/\`newton_bond\` branches out of the serial \`bond_class2\` and \`angle_harmonic… |
+------+-------------------------------------------------+---------------------+---------+-------------------------------------------------------------------------------------------------------------+
| 2    | `25f839bdfc42 <iter-0002-table-25f839bdfc42_>`_ | correctness_failure | 0.28546 | Reuse per-type coefficients and reciprocals, evaluate the class2 polynomial in a lower-overhead f…          |
+------+-------------------------------------------------+---------------------+---------+-------------------------------------------------------------------------------------------------------------+
| 3    | `4a79c588ea08 <iter-0003-table-4a79c588ea08_>`_ | accepted            | 0.27301 | Add a \`bond_class2\` single-bond-type fast path that dispatches to a dedicated \`eval_one_type\` hel…      |
+------+-------------------------------------------------+---------------------+---------+-------------------------------------------------------------------------------------------------------------+
| 4    | `e5c1ddf090fa <iter-0004-table-e5c1ddf090fa_>`_ | rejected            | 0.27052 | Add an \`angle_harmonic\` single-angle-type fast path that dispatches to \`eval_one_type\` and hoists…      |
+------+-------------------------------------------------+---------------------+---------+-------------------------------------------------------------------------------------------------------------+
| 5    | `224da0125669 <iter-0005-table-224da0125669_>`_ | rejected            | 0.27041 | Cache per-bond force components in bond_class2 and restore the single-angle-type fast path in ang…          |
+------+-------------------------------------------------+---------------------+---------+-------------------------------------------------------------------------------------------------------------+
| 6    | `e639b3949f97 <iter-0006-table-e639b3949f97_>`_ | rejected            | 0.27144 | Accumulate consecutive same-\`i1\` bonds in the single-type \`bond_class2\` fast path to reuse the le…      |
+------+-------------------------------------------------+---------------------+---------+-------------------------------------------------------------------------------------------------------------+
| 7    | `47d88a349494 <iter-0007-table-47d88a349494_>`_ | rejected            | 0.27624 | Restore the \`angle_harmonic\` single-angle-type fast path and add an explicit same-\`i1\` two-bond b…      |
+------+-------------------------------------------------+---------------------+---------+-------------------------------------------------------------------------------------------------------------+
| 8    | `e3c2ae9ee96a <iter-0008-table-e3c2ae9ee96a_>`_ | rejected            | 0.27111 | Restore \`angle_harmonic\` single-angle-type dispatch and tighten one-type bonded hot paths by only…        |
+------+-------------------------------------------------+---------------------+---------+-------------------------------------------------------------------------------------------------------------+
| 9    | `6fbda0abc94b <iter-0009-table-6fbda0abc94b_>`_ | rejected            | 0.26892 | Add a dedicated two-way unrolled \`bond_class2\` one-type hot path for the dominant no-\`evflag\`, Ne…      |
+------+-------------------------------------------------+---------------------+---------+-------------------------------------------------------------------------------------------------------------+
| 10   | `cf4ba8e7204d <iter-0010-table-cf4ba8e7204d_>`_ | accepted            | 0.26621 | Combine the correctness-safe \`bond_class2\` one-type two-way unrolled hot kernel with the correctn…        |
+------+-------------------------------------------------+---------------------+---------+-------------------------------------------------------------------------------------------------------------+
| 11   | `672cbe901cae <iter-0011-table-672cbe901cae_>`_ | accepted            | 0.24606 | Add a dedicated \`angle_harmonic\` one-type no-\`evflag\`, Newton-on hot path with scalar force tempo…      |
+------+-------------------------------------------------+---------------------+---------+-------------------------------------------------------------------------------------------------------------+
| 12   | `0826ffd39444 <iter-0012-table-0826ffd39444_>`_ | accepted            | 0.23436 | Cache one-type bond_class2 bond geometry and reuse paired bond data in angle_harmonic's dominant …          |
+------+-------------------------------------------------+---------------------+---------+-------------------------------------------------------------------------------------------------------------+
| 13   | `30584dbe2ab0 <iter-0013-table-30584dbe2ab0_>`_ | accepted            | 0.22529 | Replace the one-type bond-to-angle reuse with an angle-oriented cache built in \`bond_class2\`, the…        |
+------+-------------------------------------------------+---------------------+---------+-------------------------------------------------------------------------------------------------------------+
| 14   | `0890ed81e853 <iter-0014-table-0890ed81e853_>`_ | rejected            | 0.26506 | Replace the one-type bond-to-angle geometry cache with a prepared per-angle force cache built in …          |
+------+-------------------------------------------------+---------------------+---------+-------------------------------------------------------------------------------------------------------------+
| 15   | `fdcd68522a83 <iter-0015-table-fdcd68522a83_>`_ | correctness_failure | 0.21063 | Compact the angle hot cache to bond deltas plus inverse bond lengths and consume them in the cach…          |
+------+-------------------------------------------------+---------------------+---------+-------------------------------------------------------------------------------------------------------------+
| 16   | `1f058192fe33 <iter-0016-table-1f058192fe33_>`_ | rejected            | 0.23245 | Compact the angle hot cache to bond deltas plus exact bond lengths, recomputing cached-path rsq v…          |
+------+-------------------------------------------------+---------------------+---------+-------------------------------------------------------------------------------------------------------------+
| 17   | `738ffb9e3a30 <iter-0017-table-738ffb9e3a30_>`_ | rejected            | 0.22321 | Precompute the one-type bond/angle pairing layout once per topology build and reuse it in bond_cl…          |
+------+-------------------------------------------------+---------------------+---------+-------------------------------------------------------------------------------------------------------------+
| 18   | `bed99cfa6486 <iter-0018-table-bed99cfa6486_>`_ | rejected            | 0.22343 | Prepare and reuse a prevalidated one-type bond-angle orientation layout in \`bond_class2\` hot cach…        |
+------+-------------------------------------------------+---------------------+---------+-------------------------------------------------------------------------------------------------------------+
| 19   | `acd350270fd3 <iter-0019-table-acd350270fd3_>`_ | rejected            | 0.22126 | Compact the cached angle hot entry to store the exact bond-length product \`r12 = r1*r2\` instead o…        |
+------+-------------------------------------------------+---------------------+---------+-------------------------------------------------------------------------------------------------------------+
| 20   | `c0a62a11cfaf <iter-0020-table-c0a62a11cfaf_>`_ | rejected            | 0.22182 | Combine prevalidated bond-angle orientation layout reuse in \`bond_class2\` with exact \`r12\` hot-ca…      |
+------+-------------------------------------------------+---------------------+---------+-------------------------------------------------------------------------------------------------------------+
| 21   | `e4103f7dab04 <iter-0021-table-e4103f7dab04_>`_ | correctness_failure | 0.21864 | Store exact \`r12 = r1*r2\` in the hot angle cache and reuse a shared local \`inv_r12\` in the cached…      |
+------+-------------------------------------------------+---------------------+---------+-------------------------------------------------------------------------------------------------------------+
| 22   | `efcb4c8f585c <iter-0022-table-efcb4c8f585c_>`_ | rejected            | 0.22213 | Compact the one-type hot angle cache to exact \`r12 = r1*r2\` and add prepared hot-angle index/layo…        |
+------+-------------------------------------------------+---------------------+---------+-------------------------------------------------------------------------------------------------------------+
| 23   | `567a79cfb769 <iter-0023-table-567a79cfb769_>`_ | rejected            | 0.22195 | Widen the cached \`angle_harmonic\` hot consumer to four-way exact unrolling while preserving the i…        |
+------+-------------------------------------------------+---------------------+---------+-------------------------------------------------------------------------------------------------------------+
| 24   | `2c69753dfc48 <iter-0024-table-2c69753dfc48_>`_ | rejected            | 0.22158 | Combine exact \`r12\` hot-angle cache compaction in \`bond_class2\` with a four-way cached \`angle_har…     |
+------+-------------------------------------------------+---------------------+---------+-------------------------------------------------------------------------------------------------------------+
| 25   | `bf07f28c4118 <iter-0025-table-bf07f28c4118_>`_ | accepted            | 0.21875 | Compact the hot angle cache to exact r12 payloads and widen the one-type cached bond_class2 produ…          |
+------+-------------------------------------------------+---------------------+---------+-------------------------------------------------------------------------------------------------------------+
| 26   | `a0800c8480fc <iter-0026-table-a0800c8480fc_>`_ | rejected            | 0.22115 | Restore prepared hot-angle layout/index reuse in bond_class2 and a four-way cached angle_harmonic…          |
+------+-------------------------------------------------+---------------------+---------+-------------------------------------------------------------------------------------------------------------+

Accepted Commits
----------------


Accepted candidate detail pages and current manual-review status:

+-----------------------------------------------------+----------------------------------------+
| accepted commit                                     | Human verification                     |
+=====================================================+========================================+
| :doc:`911f06a50100 <iterations/iter_0001_accepted>` | not verified                           |
+-----------------------------------------------------+----------------------------------------+
| :doc:`4a79c588ea08 <iterations/iter_0003_accepted>` | not verified                           |
+-----------------------------------------------------+----------------------------------------+
| :doc:`cf4ba8e7204d <iterations/iter_0010_accepted>` | not verified                           |
+-----------------------------------------------------+----------------------------------------+
| :doc:`672cbe901cae <iterations/iter_0011_accepted>` | not verified                           |
+-----------------------------------------------------+----------------------------------------+
| :doc:`0826ffd39444 <iterations/iter_0012_accepted>` | not verified                           |
+-----------------------------------------------------+----------------------------------------+
| :doc:`30584dbe2ab0 <iterations/iter_0013_accepted>` | not verified                           |
+-----------------------------------------------------+----------------------------------------+
| :doc:`bf07f28c4118 <iterations/iter_0025_accepted>` | not verified                           |
+-----------------------------------------------------+----------------------------------------+

.. toctree::
   :maxdepth: 1
   :hidden:

   iterations/iter_0001_accepted
   iterations/iter_0003_accepted
   iterations/iter_0010_accepted
   iterations/iter_0011_accepted
   iterations/iter_0012_accepted
   iterations/iter_0013_accepted
   iterations/iter_0025_accepted

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
     - id: train-16r-short
       weight: 1.0
       input_script: in.tip4p_nve
       data_file: water_216_data.lmp
       mpi_ranks: 16
       omp_num_threads: 1
       expected_atoms: 41472
       expected_bonds: 27648
       expected_angles: 13824
       run_steps: 1200
       thermo_every: 100
       timer_mode: full
       timestep: 0.5
       pair_style: lj/cut/tip4p/long 1 2 1 1 0.278072379 17.007
       bond_style: class2
       angle_style: harmonic
       kspace_style: pppm/tip4p 0.0001
       neighbor: 2.0 bin
       replicate:
         - 4
         - 4
         - 4
     - id: train-32r-short
       weight: 1.0
       input_script: in.tip4p_nve
       data_file: water_216_data.lmp
       mpi_ranks: 32
       omp_num_threads: 1
       expected_atoms: 41472
       expected_bonds: 27648
       expected_angles: 13824
       run_steps: 1200
       thermo_every: 100
       timer_mode: full
       timestep: 0.5
       pair_style: lj/cut/tip4p/long 1 2 1 1 0.278072379 17.007
       bond_style: class2
       angle_style: harmonic
       kspace_style: pppm/tip4p 0.0001
       neighbor: 2.0 bin
       replicate:
         - 4
         - 4
         - 4
     - id: train-16r-long
       weight: 1.0
       input_script: in.tip4p_nve_long
       data_file: water_216_data.lmp
       mpi_ranks: 16
       omp_num_threads: 1
       expected_atoms: 41472
       expected_bonds: 27648
       expected_angles: 13824
       run_steps: 10000
       thermo_every: 100
       timer_mode: full
       timestep: 0.5
       pair_style: lj/cut/tip4p/long 1 2 1 1 0.278072379 17.007
       bond_style: class2
       angle_style: harmonic
       kspace_style: pppm/tip4p 0.0001
       neighbor: 2.0 bin
       replicate:
         - 4
         - 4
         - 4

Controller reviews run the ``test-*`` benchmark cases below to validate accepted candidates:

.. code-block:: yaml

   cases:
     - id: test-32r-long
       weight: 0.5
       input_script: in.tip4p_nve_long
       data_file: water_216_data.lmp
       mpi_ranks: 32
       omp_num_threads: 1
       expected_atoms: 41472
       expected_bonds: 27648
       expected_angles: 13824
       run_steps: 10000
       thermo_every: 100
       timer_mode: full
       timestep: 0.5
       pair_style: lj/cut/tip4p/long 1 2 1 1 0.278072379 17.007
       bond_style: class2
       angle_style: harmonic
       kspace_style: pppm/tip4p 0.0001
       neighbor: 2.0 bin
       replicate:
         - 4
         - 4
         - 4
     - id: test-64r-short
       weight: 0.25
       input_script: in.tip4p_nve
       data_file: water_216_data.lmp
       mpi_ranks: 64
       omp_num_threads: 1
       expected_atoms: 41472
       expected_bonds: 27648
       expected_angles: 13824
       run_steps: 1200
       thermo_every: 100
       timer_mode: full
       timestep: 0.5
       pair_style: lj/cut/tip4p/long 1 2 1 1 0.278072379 17.007
       bond_style: class2
       angle_style: harmonic
       kspace_style: pppm/tip4p 0.0001
       neighbor: 2.0 bin
       replicate:
         - 4
         - 4
         - 4


.. _summary-baseline-e7c0ed95a333: https://github.com/skilled-scipkg/lammps/commit/e7c0ed95a333
.. _summary-best-bf07f28c4118: https://github.com/skilled-scipkg/lammps/commit/bf07f28c4118
.. _iter-0000-table-e7c0ed95a333: https://github.com/skilled-scipkg/lammps/commit/e7c0ed95a333
.. _iter-0001-table-911f06a50100: https://github.com/skilled-scipkg/lammps/commit/911f06a50100
.. _iter-0002-table-25f839bdfc42: https://github.com/skilled-scipkg/lammps/commit/25f839bdfc42
.. _iter-0003-table-4a79c588ea08: https://github.com/skilled-scipkg/lammps/commit/4a79c588ea08
.. _iter-0004-table-e5c1ddf090fa: https://github.com/skilled-scipkg/lammps/commit/e5c1ddf090fa
.. _iter-0005-table-224da0125669: https://github.com/skilled-scipkg/lammps/commit/224da0125669
.. _iter-0006-table-e639b3949f97: https://github.com/skilled-scipkg/lammps/commit/e639b3949f97
.. _iter-0007-table-47d88a349494: https://github.com/skilled-scipkg/lammps/commit/47d88a349494
.. _iter-0008-table-e3c2ae9ee96a: https://github.com/skilled-scipkg/lammps/commit/e3c2ae9ee96a
.. _iter-0009-table-6fbda0abc94b: https://github.com/skilled-scipkg/lammps/commit/6fbda0abc94b
.. _iter-0010-table-cf4ba8e7204d: https://github.com/skilled-scipkg/lammps/commit/cf4ba8e7204d
.. _iter-0011-table-672cbe901cae: https://github.com/skilled-scipkg/lammps/commit/672cbe901cae
.. _iter-0012-table-0826ffd39444: https://github.com/skilled-scipkg/lammps/commit/0826ffd39444
.. _iter-0013-table-30584dbe2ab0: https://github.com/skilled-scipkg/lammps/commit/30584dbe2ab0
.. _iter-0014-table-0890ed81e853: https://github.com/skilled-scipkg/lammps/commit/0890ed81e853
.. _iter-0015-table-fdcd68522a83: https://github.com/skilled-scipkg/lammps/commit/fdcd68522a83
.. _iter-0016-table-1f058192fe33: https://github.com/skilled-scipkg/lammps/commit/1f058192fe33
.. _iter-0017-table-738ffb9e3a30: https://github.com/skilled-scipkg/lammps/commit/738ffb9e3a30
.. _iter-0018-table-bed99cfa6486: https://github.com/skilled-scipkg/lammps/commit/bed99cfa6486
.. _iter-0019-table-acd350270fd3: https://github.com/skilled-scipkg/lammps/commit/acd350270fd3
.. _iter-0020-table-c0a62a11cfaf: https://github.com/skilled-scipkg/lammps/commit/c0a62a11cfaf
.. _iter-0021-table-e4103f7dab04: https://github.com/skilled-scipkg/lammps/commit/e4103f7dab04
.. _iter-0022-table-efcb4c8f585c: https://github.com/skilled-scipkg/lammps/commit/efcb4c8f585c
.. _iter-0023-table-567a79cfb769: https://github.com/skilled-scipkg/lammps/commit/567a79cfb769
.. _iter-0024-table-2c69753dfc48: https://github.com/skilled-scipkg/lammps/commit/2c69753dfc48
.. _iter-0025-table-bf07f28c4118: https://github.com/skilled-scipkg/lammps/commit/bf07f28c4118
.. _iter-0026-table-a0800c8480fc: https://github.com/skilled-scipkg/lammps/commit/a0800c8480fc