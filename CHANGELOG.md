## v0.9

Initial public version.

## v1.0

### New features
- Added `--pad` option to the `segment` CLI command
- Added `total_cn` ploidy output to `main_ploidy`
- Arms labelling funcion `add_arms` added to segments API.

### Infrastructure
- Added GitHub Actions CI with automated test suite
- Test badge added to README

## v1.1.0

### New features
- Added `remove_overlaps` for removing overlaps between consecutive copy-number segments. Overlapping bases are removed equally from both segments, with an odd remaining base removed from the first segment.

### Fixes
- `main_align` now removes segment overlaps before merging neighbouring segments.
- Overlap removal now raises `ValueError` if trimming would completely remove either segment.

## v1.1.1

### Infrastructure
- Relaxed the dependency requirements. `pandas>=2.2` was never a real minimum and is now
  `pandas>=1.5` (verified down to 1.3.5), `numpy>=1.23`, and `numba>=0.57` (0.56.4 fails
  to compile `np.round` of a 2-D array in the clustering kernel).
- matplotlib is now an optional dependency, installed via the `plot` extra
  (`pip install "CNSistent[plot]"`). It is only needed by the plotting API, not by any CLI
  command, and it accounted for ~83 MB of the install. `cns.analyze.plot` is imported on
  first access, so `cns.fig_lines` and friends keep working unchanged. Note that plotting
  names are no longer star-exported, `from cns import *` no longer brings them into scope.
- numba is still a default dependency, but is no longer required to import the package.
  When it is absent the JIT decorators become no-ops and a warning is logged. Aggregation
  is then about 2x slower.
- Added tests for the plotting API and CI now installs the `plot` extra so they run.
  Since the module is imported lazily, nothing else in the suite loads it, so without
  these tests a breakage in `cns.analyze.plot` would only surface for users.
