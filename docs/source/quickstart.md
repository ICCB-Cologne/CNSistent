# Quickstart

## Processing bundled data

The `cnsistent` repository contains data from the PCAWG, TCGA, and TRACERx studies, which can be preprocessed and segmented using bundled scripts under the following steps:

1. Start a `bash`-compatible shell.	
2. Create the conda environment: `conda env create -f ./cnsistent.yml`
3. Activate the environment: `conda activate cns`
4. Install the package from location: `pip install -e .`
5. Process data: `bash ./scripts/data_process.sh`
6. Optional: Aggregate processed data: `bash ./scripts/data_aggregate.sh`

**NOTES:**

* The aggregation step will take several hours and produce ~40GB of data.
* By default, 16 threads are used, if that causes problems (crashes), reduce the number of threads in the `data_process.sh` and `data_aggregate.sh` scripts.
* You can also install the package with `pip install .`, however there is a set of utility functions for loading data in `cns.data_utils.py` that will not be accesible then.

## Accessing bundled data 
TODO

## Basic tool usage

`cnsistent` SCNA profiles as `.tsv` files. Have an example file `data.tsv`:

```data.tsv
sample_id   	chrom   start   end     total_cn
sample1         chr1    100     200     1       
...
```

To preprocess the segments:
```bash
cns fill data.tsv --out filled.tsv
cns impute filled.tsv --out imputed.csv
```

To create statistics:
```bash
cns coverage data.tsv --out samples.tsv
cns ploidy imputed.tsv --samples samples.tsv --out samples.tsv
cns signatures imputed.tsv --samples samples.tsv --out samples.tsv
```

To calculate the mean ploidy per chromosome arm:
```bash
cns segment arms --out arms.bed
cns aggregate imputed.tsv --segments --out a_bins.tsv
```

To conduct breakpoint clustering with 1 mb distance:
```bash
cns cluster imputed.tsv --merge 1000000 --out clust.bed
cns aggregate imputed.tsv  --segments clust.bed --out c_bins.tsv
```

To conduct segmentation using 5 mb bins:
```bash
cns segment whole --step 5000000 --out clust.bed
cns aggregate data.tsv  --segments clust.bed --out c_bins.tsv
```

## Basic library usage
TODO