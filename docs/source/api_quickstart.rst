API Quickstart
==============

The repository contains a set of utility functions for loading the available data in `cns.data_utils`.

* `main_load_data`: The entry point function. It loads and filters samples from the datasets (PCAWG, TCGA, TRACERx) together with the imputed CNS data. If bins value is provided, instead of CNS data, binned data is loaded for the perspective bin size.
* `load_COSMIC`: Loads the gene set from COSMIC with columns `gene, chrom, start, end`.
* `load_ENSEMBL`: Loads the gene set from Ensembl with columns `gene, chrom, start, end`.

### Data subfunctions:

* `filter_samples`: Filters samples based on the provided criteria.
* `load_all_samples`: Loads samples only, with possible filtering.
* `load_merged_samples`: Same as `load_all_samples`, but merges the samples from the datasets into a dataframe.
* `load_merged_bins`: Only loads bins for provided list of samples.
* `load_merged_cns`: Only loads CNS for provided list of samples.

### Utility functions

* `get_cns_for_type`: Obtains CNS for samples for a given cancer type.

### File access:
* `out_path`: Returns the default output path where output of processing/binning files are stored.
* `data_path`: Returns the default data path where the input files (bundled together with the repository) are stored.
* `load_cns_out`: Loads CNS file from the output folder.
* `load_samples_out`: Loads samples file from the output folder.
* `load_bins`: Loads bins from the output folder.