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
adjacent-temperature swap statistics. These explicit writers remain available. Streaming output, described below,
automatically refreshes the production snapshots during sampling. MCMC ROOT trees also contain a Replica
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
returning cold. Replica IDs are local to each production run. Without streaming output the full hot
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

This runs statistical and streaming regression tests, including healthy/shifted/stuck/autocorrelated chains,
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

All-temperature streaming diagnostics
--------------------------------------
Rebuild and install BAT, then rebuild every application against those same
headers and libraries. The engine layout and public API have changed; an
old mainGM/HEPfit binary cannot be mixed with the new library. CMake exports
BAT_HAS_PT_STREAMING_DIAGNOSTICS; bat-config --features includes pt-streaming.

Before initialization, after SetParallelTempering, clients enable output:
  model.SetPTDiagnosticsOutput("results/run_", 1000, 100);
The first interval is the number of full sweeps between saved counter blocks.
The second is the state trace stride; 0 disables traces but retains counters.
An empty prefix disables tracking. Parent directories must already exist.
Do not share a prefix between simultaneous runs. Reusing a prefix replaces
its previous diagnostics. Output settings are intentionally not copied when
a model is copied. Only rank 0 writes files.

Saved counters cover every temperature and independent ladder:
  pt_local_moves.csv          Proposal attempts, accepts, invalid proposals,
                              nonfinite evaluations; factorized/joint separate
  pt_region_occupancy.csv     Each partition/region's post-sweep occupancy
  pt_region_transitions.csv  All region changes, cause=local or cause=swap
  pt_swaps_by_block.csv       Adjacent-link attempts/accepts in each time block
  pt_roundtrip_regions.csv   Region at cold departure and cold return, by walker
  pt_partitions.csv          Model-provided partition and region label schema
  pt_states.csv              Decimated slot, walker ID, T, log L, log prior,
                              and region indices (not full parameter vectors)
  pt_checkpoints.csv         Commit ledger for complete counter blocks
  pt_progress_PHASE.csv      Latest saved sweep and phase completion status

PHASE is prerun (adaptive tuning), warmup (discarded loaded-proposal sweeps)
or production. Sweep numbering and block counters restart in each phase.
An ordinary run has prerun and production; a loaded-proposal run has warmup
and production. Never combine prerun/warmup with posterior draws. Trace
files include sweep 0 and the final sweep when enabled. Occupancy counts
start at sweep 1, after local updates and any scheduled swaps.

A factorized attempt is one coordinate proposal; a multivariate attempt is
one joint proposal. A sweep contains all coordinate proposals or one joint
proposal per slot. Local crossings count every accepted proposal's label
change, including brief excursions invisible in decimated traces. A rejected
proposal cannot change a region. Accepted swaps can change both slots'
labels and are counted separately. A physical walker carries its identity
through swaps; proposals and RNGs stay with slots. Same-label exchanges
still count as accepted swaps but cause no region transition.

Round-trip region pairs use the actual last cold departure before visiting
the hottest level and returning to cold. Hot-origin walkers must first visit
cold. Different endpoint labels establish a change during the excursion;
they do not locate the change in temperature or establish that the hottest
local kernel crossed. Same-region returns are also retained. Production
journey counters reset independently of prerun/warmup.

Each saved block is flushed and closed before its ledger entry is appended.
The analyzer ignores uncommitted trailing rows after an interruption. At
most the unfinished checkpoint window is lost from the exact counters.
The legacy production swap/round-trip snapshots are also refreshed at these
checkpoints using atomic file replacement. This protects against process
interruption; it is not a power-loss durability or sampler-restart guarantee.
Fixed-size sampler counters do not grow with run length. Disk usage grows
with saved blocks and decimated traces. Full-run ROOT statistical analysis
still loads cold samples into memory; --pt-only streams the CSV counters.

Model-independent hooks
-----------------------
BAT has no model-specific parameter names or cuts. A model may override:
  DefinePTDiagnosticPartitions() const
  ClassifyPTDiagnosticRegions(parameters, labels) const
The first returns unique partition names with unique nonempty label lists.
The classifier fills exactly one valid unsigned label index per partition,
in the same order. It must be deterministic and must not alter model state,
RNGs, likelihood caches or observables. Cache parameter indices when defining
the schema. With no partitions, generic proposal/swap/trajectory diagnostics
still work; regional mixing cannot be evaluated. Label counts are capped at
1024 per partition and transition counters scale quadratically in that count.

Analyze streaming diagnostics
-----------------------------
The normal cold ROOT analysis automatically discovers these files from its
transport prefix. It adds streaming_pt to summary.json and flags transport
problems alongside the existing cold-chain checks. Explicit --regions JSON
still controls posterior-region ESS/MCSE; sampler labels describe transport.
For an ongoing/crashed run, or without a readable cold ROOT file:

  python diagnostics/analyze_mcmc.py --pt-only \
    --transport-prefix results/run_ --output results/run_transport

This reads only committed windows. A partial phase is reported explicitly.
A completed transport-only report has status transport_checked, never
screening_passed: cold-chain R-hat/ESS/posterior weights were not checked.
--fail-on-warning exits 2 for transport-only reports as those checks are
incomplete. --no-plots gives tables only. Older runs without the streaming
files remain supported by the ordinary cold-chain analysis, with
streaming_pt.available=false. The new events cannot be recovered afterwards.

New report tables (all phases, runs, ladders and partitions):
  pt_local_acceptance.csv    Cumulative proposal counters/rates by kernel/T
  pt_swap_acceptance.csv     Cumulative link counters/rates per ladder
  pt_occupancy.csv           Occupancy by temperature and independent ladder
  pt_region_crossings.csv    Exact local versus swap-induced label changes
  pt_region_returns.csv      Round-trip endpoint pairs and changed_region
  pt_*_blocks.csv            Committed source counter windows, with run ID
Plots: pt_local_acceptance.png, pt_swap_history.png, pt_regions_runN.png,
pt_paths_runN.png. Region plots show at most eight partitions; trajectory
plots show at most eight walkers from the first ladder. Tables retain all.
Trajectory and block-history plot points are bounded/decimated; event counts
are exact through the last committed sweep. The hottest/cold diagnostic views
prefer production, falling back to prerun or warmup while a run is ongoing.

Interpret the findings together:
- No hot local crossings plus disagreeing hot-ladder occupancies suggests
  inadequate exploration at the top of the ladder. No hot crossings alone
  does not establish trapping: one region may dominate that tempered target.
- Hot crossings with separated, noncrossing cold ladders suggests a transport
  problem. Inspect swap histories, intermediate occupancies and replica paths.
- Round trips that repeatedly return to the same region can coexist with
  good temperature movement and poor global mixing. Changed-region returns
  are useful evidence, but still require stable cold weights and R-hat/ESS.
Heated occupancies belong to different tempered distributions. Never pool
hot samples or interpret their weights as cold posterior probabilities.

Regression coverage
--------------------
PT_streaming tests exact local/swap attribution, brief crossings, region
pairs on round trips, interrupted checkpoints, phase isolation, disabled
traces, bounded counters, and identical fixed-seed sample paths with tracking
on/off for factorized and multivariate kernels. Proposal_pt covers the
separate loaded warmup phase. The Python suite includes streaming reports,
committed-window validation, partial trailing output, phase/run identities,
transport-only analysis and plotting, alongside the statistical tests.
