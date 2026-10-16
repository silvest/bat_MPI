# Build and install bat_MPI with CMake

CMake is the primary build system for this MPI fork. Autotools is deprecated
and retained temporarily; its historical instructions are in
[doc/INSTALL_AUTOTOOLS.md](doc/INSTALL_AUTOTOOLS.md).

## Dependencies

Use CMake 3.16 or newer, a C++ compiler compatible with your ROOT build, ROOT 6.08
or newer with MathMore/Minuit, and an MPI implementation such as OpenMPI. ROOT's
required C++ standard is used by default. Use the same compiler, C++ runtime,
ROOT and MPI installation for BAT and downstream applications. Python is optional
and only needed for offline MCMC diagnostics and their tests.

## Build

```sh
git clone https://github.com/silvest/bat_MPI.git
cd bat_MPI
cmake -S . -B build-cmake -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX="$HOME/opt/bat-pt"
cmake --build build-cmake -j 4
cmake --install build-cmake
```

For a nonstandard ROOT installation, pass `-DROOT_DIR=/path/to/root/lib/cmake/ROOT`
or add its prefix to `CMAKE_PREFIX_PATH`. Select MPI with CMake's `MPI_CXX_COMPILER`
or your module environment. A plain C++ compiler works: the exported BAT targets
carry MPI include and link dependencies. Use separate build directories when
changing ROOT, MPI, compiler, optional features or install layout. In-source builds
are intentionally rejected so an old autotools `config.h` cannot affect CMake.

The default build installs shared and static `BAT`, `BATmodels` and `BATmtf`
libraries, public headers, ROOT dictionaries, the CMake package, compatibility
`bat-config`/`bat.pc`, the `bat-project` generator and offline diagnostic scripts.
GNUInstallDirs controls `lib` versus `lib64`. To use a consistent cluster prefix,
configure explicitly with `-DCMAKE_INSTALL_LIBDIR=lib` if desired.

## Options

| Option | Default | Purpose |
| --- | --- | --- |
| `BAT_BUILD_SHARED_LIBS` | ON | Shared libraries. |
| `BAT_BUILD_STATIC_LIBS` | ON | Static libraries; at least one library form is required. |
| `BAT_BUILD_MODELS` | ON | General fitters in `BATmodels`. |
| `BAT_BUILD_MTF` | ON | Template fitters in `BATmtf`. |
| `BAT_ENABLE_ROOSTATS` | OFF | RooFit/RooStats adapter in BATmodels; requires compatible ROOT RooFit/RooStats components. |
| `BAT_ENABLE_CUBA` | OFF | Existing Cuba 3.3/4.x; pass `Cuba_ROOT` or `Cuba_INCLUDE_DIR`/`Cuba_LIBRARY`. No automatic downloads. |
| `BAT_ENABLE_OPENMP` | OFF | Legacy OpenMP thread bookkeeping; MPI remains the sampler backend. |
| `BAT_BUILD_PT_TESTS` | OFF | PT, transport and proposal-cache regression tests. |
| `BAT_BUILD_INSTALL_TESTS` | OFF | Relocated installation, exported targets, bat-config and generated project tests. |
| `BAT_BUILD_LEGACY_TESTS` | OFF | Upstream unit tests, running one MPI process. |
| `BAT_BUILD_EXAMPLES` | OFF | Standalone C++ examples, running one MPI process. |
| `BAT_BUILD_DOCS` | OFF | Enable `BAT_reference` Doxygen target. |

The optional Cuba interface is detected by compilation; a requested feature with
missing dependencies fails configuration. Unlike MPI, enabling the legacy OpenMP
hooks does not parallelize likelihood evaluation across additional workers.
Historical ROOT `.C` example macros and the LaTeX manual remain source resources;
they are not automatically executed by the CMake build.

## CMake consumers

```cmake
find_package(BAT 1 CONFIG REQUIRED COMPONENTS PT)
add_executable(my_analysis main.cpp)
target_link_libraries(my_analysis PRIVATE BAT::BAT)
# Optional: BAT::BATmodels and BAT::BATmtf.
set_property(TARGET my_analysis APPEND PROPERTY INSTALL_RPATH "${BAT_LIBRARY_DIR}")
```

Configure the consumer with `-DCMAKE_PREFIX_PATH="$HOME/opt/bat-pt"` or the exact
`BAT_DIR` printed by installation (`lib/cmake/BAT` or `lib64/cmake/BAT`). Components
`models`, `mtf`, `RooStats` and `PT` can be required explicitly. When both forms are
built, the default targets select shared libraries; `BAT::BAT_static`,
`BAT::BATmodels_static` and `BAT::BATmtf_static` select the static alternatives.
With a static-only build, the default targets select static libraries. ROOT and
MPI can still be shared dependencies in a static BAT build.

Applications must call `MPI_Init` before constructing a BAT engine and destroy
engines before `MPI_Finalize`. Existing upstream examples use one process and do
not implement a distributed likelihood worker loop. A multi-process application
must implement that protocol, as HEPfit and GoldenModes do.

`BAT_HAS_PARALLEL_TEMPERING` and `BAT_HAS_TRANSPORT_DIAGNOSTICS` describe this fork's
APIs. Headers and libraries must come from the same installation.

## Existing consumers and cluster jobs

`bat-config` supports `--prefix`, `--version`, `--cflags`, `--libs`, `--bindir`,
`--incdir`, `--libdir`, and `--features`. Compile legacy consumers with `mpicxx`;
`--libs` includes an explicit BAT runtime search path. `bat.pc` supports
pkg-config consumers. The preferred CMake package propagates dependencies directly.

GoldenModes now prefers `find_package(BAT CONFIG)` and preserves the chosen BAT
runtime path in installed `mainGM`. Its explicit `BAT_SOURCE_DIR`/`BAT_LIBRARY`
mode and the bat-config fallback remain available. Use a fresh build directory:

```sh
cmake -S /path/to/ThErrorsInGoldenModes -B /path/to/ThErrorsInGoldenModes/build-cmake \
  -DCMAKE_PREFIX_PATH="$HOME/opt/bat-pt" \
  -DCMAKE_INSTALL_PREFIX=/path/to/ThErrorsInGoldenModes/dist
cmake --build /path/to/ThErrorsInGoldenModes/build-cmake -j 4
cmake --install /path/to/ThErrorsInGoldenModes/build-cmake
cmake -S /path/to/HEPfit -B /path/to/HEPfit/build-pt \
  -DMPIBAT=ON -DBAT_INSTALL=OFF -DBAT_INSTALL_DIR="$HOME/opt/bat-pt" \
  -DRGESOLVER_CONFIG_DIR=/path/to/RGESolver/bin
cmake --build /path/to/HEPfit/build-pt --target analysis -j 4
```

On Linux, an old BAT directory in `LD_LIBRARY_PATH` can override executable
RUNPATH. Set the correct environment in the submitted job script on every node,
not only in the login shell, and check the installed executable:

```sh
export PATH="$HOME/opt/bat-pt/bin:$PATH"
export LD_LIBRARY_PATH="$HOME/opt/bat-pt/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
bat-config --prefix --features
ldd /path/to/ThErrorsInGoldenModes/dist/bin/mainGM | grep BAT
nm -D "$HOME/opt/bat-pt/lib/libBAT.so" | c++filt | grep SetParallelTempering
```

Use `lib64` if that is your configured library directory. The installation is
relocatable; consumers must be reconfigured when its prefix changes. No build
system can prevent an explicitly conflicting loader environment from selecting
an older ABI-compatible-named library.

## Tests and examples

```sh
cmake -S . -B build-cmake -DBAT_BUILD_PT_TESTS=ON \
  -DBAT_BUILD_INSTALL_TESTS=ON -DBAT_BUILD_LEGACY_TESTS=ON \
  -DBAT_BUILD_EXAMPLES=ON
cmake --build build-cmake -j 4
ctest --test-dir build-cmake --output-on-failure
```

Pass `BAT_DIAGNOSTICS_PYTHON=/path/to/venv/bin/python` to include the statistical
analyzer's tests; install its dependencies from `diagnostics/requirements.txt`.
Launcher-specific arguments go in the semicolon-separated `MPIEXEC_PREFLAGS`.
For example, local OpenMPI on some macOS machines needs
`--map-by;slot;--bind-to;none;--host;localhost:8`.

On the local ROOT 6.40 installation, upstream `legacy_BCHistogram` (an exact
floating-point boundary assertion) and `legacy_BCPrior` (numeric versus histogram
mode comparison) fail. Both failures were reproduced against the pre-migration
BAT installation. They remain visible in CTest; their assertions were not relaxed.
The PT, proposal-cache, diagnostics and installation tests are independent.

The installation test moves a staging prefix, builds and runs shared/static CMake
consumers and a bat-config consumer, runs installed executables, and generates,
builds and runs new projects. It does not modify your configured install prefix.

With `BAT_BUILD_DOCS=ON`, build `BAT_reference` explicitly to generate the HTML
reference. To create a project from an installed BAT, run `bat-project MyAnalysis
MyModel`, then configure it with the BAT prefix as above. A Makefile is retained
for compatibility.

## Migration validation

The macOS/ROOT 6.40 checks cover shared and static libraries, a static-only
`lib64` installation, relocation to a prefix containing spaces, CMake and
bat-config/pkg-config consumers, generated projects, compiled examples,
RooStats/OpenMP builds, and the Doxygen reference target. GoldenModes' prior
checks and HEPfit's MPI/PT integration tests pass against the installed package.
A dangling reference in an upstream integration test was corrected without
changing its numerical assertion. The two upstream failures described above
also occur with the old BAT installation.

Cuba is not installed on the validation machine: its missing-dependency error
was checked, but a real Cuba-enabled build still needs validation. The Linux
cluster build/run must also be checked before retiring autotools entirely.
