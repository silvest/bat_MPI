BAT - Bayesian Analysis Toolkit - MPI fork
==========================================

The Bayesian Analysis Toolkit is a software analysis package which allows
to address the main goals of a typical data analysis:

 - compare model predictions with data,
 - draw conclusions on the validity of the model as a representation
   of the data and
 - extract the possible values of parameters within the context of
   a model.

The BAT is based on Bayes' Theorem and is realized with the use of Markov
Chain Monte Carlo. This gives access to the full posterior probability
distribution and enables straightforward parameter estimation, limit
setting and uncertainty propagation.

The BAT is implemented in C++ on top of CERN's ROOT with emphasis put
on flexibility and modularity in defining models while keeping in mind
the reliability and speed requirements of the numerical operations.

Authors
--------

See `doc/CREDITS` file for list of contributors to BAT. This fork was adapted for HEPfit by the HEPfit contributors.

Availability
-------------

The latest code and official bat_MPI releases are published online at
https://github.com/silvest/bat_MPI. BAT is open-source software under the LGPLv3 or later.
See `doc/COPYING` and `doc/LICENSE` for licensing details.

Documentation
-------------

CMake is the primary build system. See [INSTALL.md](INSTALL.md) for dependencies,
installation, exported CMake targets, optional features, tests and cluster setup.
Autotools is deprecated and retained during the migration.

```sh
cmake -S . -B build-cmake -DCMAKE_INSTALL_PREFIX="$HOME/opt/bat-pt"
cmake --build build-cmake -j 4
cmake --install build-cmake
```

Contents
---------

Directory `bat/`:

* `doc/`: documentation about BAT
* `examples/`: well commented example programs, see `examples/README.md`
* `BAT/`: BAT header files
* `src/`: source files of the BAT core library
* `models`: models for specific data-analysis problems: fast fitters
  in `base/` and the multi-template fitter in `mtf`
* `benchmarks`: performance test suite with html output
* `test`: unit tests
* `tools`: template script to generate a BAT project and build utilities

Other files distributed with BAT are part of the configuration and
build system.

## Production MCMC diagnostics

PT runs can record configuration round trips and per-ladder swap counts.
`diagnostics/analyze_mcmc.py` analyzes cold ROOT samples with ArviZ for modern
R-hat, ESS, MCSE, and configurable region occupancy. See
[diagnostics/README.txt](diagnostics/README.txt) for installation, examples,
interpretation, and regression tests. Python dependencies are optional and
are not needed by the sampler itself.

The optional staged warmup policy bounds proposal adaptation and adds discarded
fixed-proposal settling before production. Legacy stopping remains the default.
See [diagnostic documentation](diagnostics/README.txt) for the API, status output
and limitations; completing warmup does not establish posterior convergence.

Smooth scale adaptation is independently selectable with `SetScaleAdaptation`;
legacy scale updates remain the default. See the diagnostic documentation for
the diminishing log-scale rule and recorded comparison settings.
