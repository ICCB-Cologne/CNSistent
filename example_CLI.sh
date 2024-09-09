#/usr/bin/env bash

# %%
# Version and help
cns -v
cns -h

mkdir -p ./tests/out

# %%
# Fill in gaps in the CNS data
cns fill ./tests/in/test_cns_source.tsv --sample ./tests/in/test_sample_source.tsv --out ./tests/out/test_cns_fill.tsv --verbose

# %%
# Impute the filled regions in the CNS data
cns impute ./tests/out/test_cns_fill.tsv --out ./tests/out/test_cns_imp.tsv --verbose

# %%
# Calculate portion of the filled regions in the CNS data
cns coverage ./tests/out/test_cns_fill.tsv --out ./tests/out/test_sample_cover.tsv --verbose

# %%
# Calculate the aneuploidy of the filled-in CNS data
cns ploidy ./tests/out/test_cns_imp.tsv --samples ./tests/in/test_sample_source.tsv --out ./tests/out/test_sample_ploidy.tsv --verbose

# %%
# Obtain regions by clustering neighboring breakpoints
cns cluster ./tests/out/test_cns_fill.tsv --dist 100000 --out ./tests/out/mcs_regions.tsv --verbose 

# %%
# Calculate arm-level bins
cns bin --select arms ./tests/out/test_cns_fill.tsv --out ./tests/out/test_cns_arms.tsv --verbose

# %%
# Calculate band-level bins
cns bin --select bands ./tests/out/test_cns_fill.tsv --out ./tests/out/test_cns_bands.tsv --verbose

# %% 
# Calculate 1MB bins
cns bin --bins 1000000 ./tests/out/test_cns_fill.tsv --out ./tests/out/test_cns_1MB.tsv --verbose

# %% 
# Calculate 1MB bins with gaps
cns bin --bins 1000000 ./tests/out/test_cns_fill.tsv --out ./tests/out/test_cns_1MB_gaps.tsv --remove gaps --filter 500000 --verbose

# %% 
# Create only segments from arms where gaps bigger than 100KB are removed
cns bin --select arms ./tests/out/test_cns_fill.tsv --out ./tests/out/test_segs_arms_gaps.tsv --remove gaps --filter 100000 --segfile

# %%
# Calculate consisten segment bins
cns bin ./tests/out/test_cns_fill.tsv --select ./tests/out/mcs_regions.tsv --out ./tests/out/test_cns_mcs.tsv --verbose
