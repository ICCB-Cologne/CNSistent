#/usr/bin/env bash
set -x
trap read debug

# %%
# Version and help
python cns.py -v
python cns.py -h

# %%
# Fill in gaps in the CNA data
python cns.py fill ./tests/in/test_cna_source.tsv --sample ./tests/in/test_sample_source.tsv --out ./tests/out/test_cna_fill.tsv --verbose

# %%
# Impute the filled regions in the CNA data
python cns.py impute ./tests/out/test_cna_fill.tsv --out ./tests/out/test_cna_imp.tsv --verbose

# %%
# Calculate portion of the filled regions in the CNA data
python cns.py coverage ./tests/out/test_cna_fill.tsv --out ./tests/out/test_sample_cover.tsv --verbose

# %%
# Calculate the aneuploidy of the filled-in CNA data
python cns.py ploidy ./tests/out/test_cna_imp.tsv --samples ./tests/in/test_sample_source.tsv --out ./tests/out/test_sample_ploidy.tsv --verbose

# %%
# Obtain regions by calculating MCS
python cns.py mcs ./tests/out/test_cna_fill.tsv --dist 100000 --out ./tests/out/mcs_regions.tsv --verbose 

# %%
# Calculate arm-level bins
python cns.py bin --select arms ./tests/out/test_cna_fill.tsv --out ./tests/out/test_cna_arms.tsv --verbose

# %%
# Calculate band-level bins
python cns.py bin --select bands ./tests/out/test_cna_fill.tsv --out ./tests/out/test_cna_bands.tsv --verbose

# %% 
# Calculate 1MB bins
python cns.py bin --bins 1000000 ./tests/out/test_cna_fill.tsv --out ./tests/out/test_cna_1MB.tsv --verbose

# %% 
# Calculate 1MB bins with gaps
python cns.py bin --bins 1000000 ./tests/out/test_cna_fill.tsv --out ./tests/out/test_cna_1MB_gaps.tsv --remove gaps --filter 500000 --verbose

# %% 
# Create only segments from arms where gaps bigger than 100KB are removed
python cns.py bin --select arms ./tests/out/test_cna_fill.tsv --out ./tests/out/test_segs_arms_gaps.tsv --remove gaps --filter 100000 --onlybins

# %%
# Calculate consisten segment bins
python cns.py bin ./tests/out/test_cna_fill.tsv --select ./tests/out/mcs_regions.tsv --out ./tests/out/test_cna_mcs.tsv --verbose
