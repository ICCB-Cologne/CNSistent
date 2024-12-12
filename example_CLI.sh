#/usr/bin/env bash

# %%
# Version and help
cns -v
cns -h
cns fill -h
set -x

mkdir -p ./tests/out

# %%
# Fill in gaps in the CNS data
cns fill ./tests/in/test_cns_double.tsv --samples ./tests/in/test_samples.tsv --out ./tests/out/test_cns_fill.tsv --verbose
cns fill ./tests/in/test_cns_single.tsv --out ./tests/out/test_cns_single_fill.tsv --verbose

# %%
# Impute the filled regions in the CNS data
cns impute ./tests/out/test_cns_fill.tsv --out ./tests/out/test_cns_imp.tsv --verbose
cns impute ./tests/out/test_cns_single_fill.tsv --out ./tests/out/test_cns_single_imp.tsv --verbose

# %%
# Calculate portion of the filled regions in the CNS data
cns coverage ./tests/out/test_cns_fill.tsv --out ./tests/out/test_sample_cover.tsv --verbose
cns coverage ./tests/out/test_cns_single_fill.tsv --out ./tests/out/test_sample_single.tsv --verbose

# %%
# Calculate the aneuploidy of the filled-in CNS data
cns ploidy ./tests/out/test_cns_imp.tsv --samples ./tests/in/test_samples.tsv --out ./tests/out/test_sample_ploidy.tsv --verbose
cns ploidy ./tests/out/test_cns_single_imp.tsv --samples ./tests/out/test_sample_single.tsv --out ./tests/out/test_sample_single.tsv --verbose

# Calculate the aneuploidy of the filled-in CNS data
cns breakage ./tests/out/test_cns_imp.tsv --samples ./tests/in/test_samples.tsv --out ./tests/out/test_sample_breakage.tsv --verbose
cns breakage ./tests/out/test_cns_single_imp.tsv --samples ./tests/out/test_sample_single.tsv --out ./tests/out/test_sample_single.tsv --verbose

# %%
# Try different segmentations
cns segment ./tests/out/test_cns_fill.tsv --merge 100000 --out ./tests/out/mcs_regions.bed --verbose 
cns segment arms --out ./tests/out/test_segs_arms.bed --verbose
cns segment bands --out ./tests/out/test_segs_bands.bed --verbose
cns segment ./tests/out/test_cns_fill.tsv --split 1000000 --out ./tests/out/test_segs_1MB.bed --verbose
cns segment ./tests/out/test_cns_fill.tsv --split 1000000 --out ./tests/out/test_segs_1MB_gaps.bed --remove gaps --filter 500000 --verbose
cns segment arms --out ./tests/out/test_segs_arms_gaps.bed --remove gaps --filter 100000
cns segment ./data/COSMIC_consensus_genes.bed --out ./tests/out/test_COSMIC_gaps.bed --remove gaps --filter 100000

# %%
# Calculate consisten segmentation 
cns aggregate ./tests/out/test_cns_fill.tsv --segments ./tests/out/test_segs_1MB.bed --out ./tests/out/test_cns_1MB.tsv --verbose
# Calculate gene selection
cns aggregate ./tests/out/test_cns_imp.tsv --segments ./tests/out/test_COSMIC_gaps.bed --out ./tests/out/test_cns_COSMIC.tsv --verbose --how min
# Calculate single-column segmentation 
cns aggregate ./tests/in/test_cns_single.tsv --segments ./tests/out/test_segs_arms_gaps.bed --out ./tests/out/test_cns_arms.tsv --verbose
