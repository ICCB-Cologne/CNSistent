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

# Calculate the aneuploidy of the filled-in CNS data
cns signatures ./tests/out/test_cns_imp.tsv --samples ./tests/in/test_sample_source.tsv --out ./tests/out/test_sample_signatures.tsv --verbose

# %%
# Obtain regions by clustering neighboring breakpoints
cns segment ./tests/out/test_cns_fill.tsv --merge 100000 --out ./tests/out/mcs_regions.bed --verbose 

# %%
# Calculate arm-level segments
cns segment dummy --select arms --out ./tests/out/test_segs_arms.bed --verbose

# %%
# Calculate band-level segments
cns segment dummy --select bands ./tests/out/test_cns_fill.tsv --out ./tests/out/test_segs_bands.bed --verbose

# %% 
# Calculate 1MB segments
cns segment dummy --split 1000000 ./tests/out/test_cns_fill.tsv --out ./tests/out/test_segs_1MB.bed --verbose

# %% 
# Calculate 1MB segments with gaps
cns segment dummy --split 1000000 ./tests/out/test_cns_fill.tsv --out ./tests/out/test_segs_1MB_gaps.bed --remove gaps --filter 500000 --verbose

# %% 
# Create only segments from arms where gaps bigger than 100KB are removed
cns segment dummy --select arms ./tests/out/test_cns_fill.tsv --out ./tests/out/test_segs_arms_gaps.bed --remove gaps --filter 100000

# %% 
# Create segments for COSMIC
cns segment dummy --select ./data/COSMIC_consensus_genes.bed --out ./tests/out/test_COSMIC_gaps.bed --remove gaps --filter 100000

# %%
# Calculate consisten segment segments
cns aggregate ./tests/out/test_cns_fill.tsv --segments ./tests/out/test_segs_1MB.tsv --out ./tests/out/test_cns_1MB.tsv --verbose
