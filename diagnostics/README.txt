Production MCMC diagnostics for bat_MPI
======================================

This directory contains the generic analyzer, its Python dependencies and
a generic JSON schema example. Model-specific parameter names, partitions
and mode definitions belong in the corresponding model repository. Pass
that repository's JSON file to the analyzer with --regions.

New sampler tracking
--------------------
PT clients can call WriteParallelTemperingTransport(prefix) after production
sampling to save:
  <prefix>pt_roundtrips.csv
  <prefix>pt_swaps_by_ladder.csv
The existing WriteParallelTemperingDiagnostics(filename) saves aggregate
adjacent-temperature swap statistics. A client must call the writers to
export the accumulated statistics. MCMC ROOT trees also contain a Replica
branch identifying configurations arriving at cold temperature slots.

Physical replica IDs move with configurations through accepted exchanges;
RNGs, proposal matrices and iteration counters stay with temperature slots.
Statistics reset at production start, excluding adaptation and discarded
loaded-proposal warm-up. A hot-origin replica must visit cold before it can
start a counted cold -> hottest -> cold round trip.

Trip durations are full sampler sweeps, including time waiting at cold.
They are not per-parameter Iteration counts in factorized runs. Missing
durations are nan. journey_stage is 0 before the first cold visit, 1 after
cold but before subsequently reaching hottest, and 2 after hottest until
returning cold. Replica IDs are local to each production run. The full hot
path is not saved; round trips are accumulated inside the sampler and
cannot be reconstructed retrospectively from cold samples alone.

Install
-------
Use Python 3.10 or newer. From the bat_MPI repository:

  python3 -m venv diagnostics-env
  source diagnostics-env/bin/activate
  python -m pip install -r diagnostics/requirements.txt

Python is not a runtime dependency of the C++ sampler. Reports may be
produced on a workstation after copying the cluster's ROOT/CSV files.
The implementation uses ArviZ's maintained rank-normalized split/folded
R-hat, bulk/tail ESS and mean ESS/MCSE estimators.

Analyze
-------
  python diagnostics/analyze_mcmc.py /path/to/run_cold_chains.root \
    --regions /path/to/model/diagnostics/regions.json \
    --output /path/to/report

If the filename ends in cold_chains.root, the preceding prefix is used to
find PT transport CSVs. Otherwise provide --transport-prefix /path/to/run_.
Use --tree NAME when there is more than one MCMC tree in a ROOT file.

Multiple independent runs may be compared by giving several ROOT files.
Their chains keep separate run:chain identities. Do not treat successive
pieces of one continuing chain as independent runs. Repeated explicit
--transport-prefix arguments must follow the same order as the input files.

Only Phase=1 production samples are analyzed. Rejections remain in the
sample sequence. Chain ordering is preserved; unequal chain lengths,
duplicate/nonmonotonic/irregular iterations and explicitly heated input
are rejected. All recorded draws are used for statistical calculations;
only plotted/cumulative traces are decimated for large outputs.

All saved numerical parameter, observable and log-density branches are
considered by default. Use --var NAME repeatedly to restrict this set.
A histogram alone cannot supply ESS or MCSE because sample order is lost.
Variables ending in _arg additionally get sin/cos diagnostics; use
--angular NAME for other radian variables. Check circular transforms when
a phase branch cut makes raw-angle diagnostics misleading.

Output
------
  summary.txt / summary.json              Issues, thresholds and provenance
  diagnostics.csv                        R-hat, bulk/tail/mean ESS, MCSE
  region_weights.csv                     Pooled weights and indicator ESS/MCSE
  region_weights_by_chain.csv            Per-chain occupancy and uncertainty
  region_weights_by_block.csv            Occupancy in successive time blocks
  region_weights_cumulative.csv          Cumulative occupancy per chain
  swaps.csv / round_trips.csv             Available transport data by run
  traces.png / rank_plots.png             Parameter traces and rank distributions
  region_traces.png                      Region indicators and cumulative weights
  swap_acceptance.png / round_trips.png   Replica transport plots

Plots show up to eight variables/region entries; tables cover all requested
quantities. Flagged, nonconstant variables get priority in trace/rank plots.
Use a separate report directory for each analysis: generated files in the
chosen output directory are replaced. --no-plots produces tables only.

Model-owned region configurations
---------------------------------
A generic example is provided at diagnostics/regions.example.json. Replace
the example parameter x and its zero threshold with definitions appropriate
to your model; the tool does not automatically discover physical modes.

An explicitly disjoint regions list defines joint mode bins:
{
  "regions": [
    {"name": "negative_x", "conditions": [
      {"parameter": "x", "op": "<", "value": 0}
    ]},
    {"name": "nonnegative_x", "conditions": [
      {"parameter": "x", "op": ">=", "value": 0}
    ]}
  ]
}

Conditions within a region are ANDed. Overlap is rejected and uncovered
samples are labeled other. Explicitly declared unvisited regions remain
in the report. Supported transforms are identity, sin and cos; comparisons
are <, <=, >, >= and ==. No Python expressions are evaluated.

Alternatively, a partitions list defines separate binary diagnostic views:
{
  "partitions": [
    {"name": "real_part_sign", "conditions": [
      {"parameter": "phase", "transform": "cos", "op": ">=", "value": 0}
    ], "true_label": "positive", "false_label": "negative", "optional": true}
  ]
}

Each partition's two probabilities sum to one. Different partitions overlap:
do not add their probabilities or interpret them as joint mode weights.
Optional partitions with absent parameters are skipped and listed. An
unvisited bin outside the selected prior's support may be impossible by
construction; interpret flags with the actual sampling domain in mind.

Interpretation
--------------
Configurable defaults:
  --rhat-limit 1.01
  --ess-per-chain 100
  --region-mcse 0.01         Absolute standard error of a region weight
  --blocks 10

Region indicator ESS/MCSE refer to the mean of the raw 0/1 indicator, not
to the rank-normalized bulk ESS. A never-visited region or a chain without
region transitions is not assigned a misleading zero MCSE. Fixed parameter
metadata distinguish known constants from potentially stuck variables.
CSV nan denotes unavailable results; JSON uses null for nonfinite numbers.

Report status:
  attention          Some computed checks are flagged.
  incomplete         No computed flags, but region/transport checks missing.
  screening_passed   Configured, available checks passed the selected limits.

None establishes that every relevant mode was discovered. No universal
swap acceptance threshold is imposed. Zero swaps and missing round trips
are flagged; examine all rates, trip counts/durations, mode visits and
weight stability. A few trips alone do not establish mixing.

The tool normally exits 0 after successfully writing a report. With
--fail-on-warning, attention/incomplete exits 2. Invalid input exits 1.
Older cold ROOT files support statistical/region checks, and aggregate
swap CSVs are read as a fallback. Missing historical round trips are
explicitly reported as unavailable.

Tests
-----
  python test/test_mcmc_diagnostics.py

This runs 23 tests, including healthy/shifted/stuck/autocorrelated chains,
ROOT ordering/production filtering, fixed metadata, angular transforms,
region boundaries/overlap, transport provenance and report generation.
The C++ PTTransport_MPI_TEST covers configuration IDs, accepted/rejected
swaps, cold-hot-cold cycles, hot-origin exclusion, durations and reset.

To include the Python suite in CTest:
  cmake -S . -B build-pt -DBAT_BUILD_PT_TESTS=ON \
    -DBAT_DIAGNOSTICS_PYTHON=/absolute/path/to/diagnostics-env/bin/python
  cmake --build build-pt -j 4
  ctest --test-dir build-pt --output-on-failure

References
----------
Modern R-hat and ESS: https://arxiv.org/abs/1903.08008
PT round trips: https://arxiv.org/abs/cond-mat/0602085
ArviZ: https://python.arviz.org/
