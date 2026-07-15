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
