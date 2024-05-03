# CNSistent

CNSistent is a Python tool for processing and analyzing copy number data. It is designed to work with data from a variety of sources, including TCGA, PCAWG, TRACERx, and COSMIC. The tool is designed to be easy to use, and to provide a comprehensive set of analyses and visualizations.

# Usage

CNSistent is available both as a [stand-alone](#cli) tool and as a [Python library](#library).

## Quickstart (from Bash)

1. Create the conda environment: `conda env create -f cns.yml`
2. Activate the environment: `conda activate cns`
3. Process data: `bash ./data_process.sh`
4. Bin processed data: `bash ./data_bin.sh`
5. Run example analysis: `python example_analysis.py`


NOTES:

* By default, 15 threads are used, if that causes problems (crashes), reduce the number of threads in the `data_process.sh` and `data_bin.sh` scripts.
* The `example_analysis.py` is split into cells that can be run individually in an IDE.

# Input

The tool expects an unprocessed copy number dataset in the form of a `TSV` file with the following column scheme: `sample_id, chrom, start, end, major_cn, minor_cn`. The names of the columns can be different, but the order must be the same. 

* The `sample_id` is the identifier of the sample, 
* `chrom` is the name of the chromosome, 
* `start` and `end` are the start and end positions of the segment, 
* `major_cn` and `minor_cn` are the copy numbers of the segment.

E.g.:

```tsv
sample_id   	chrom   start   end     major_cn    minor_cn
sample1         chr1    100     200     1           0
sample1         chr1    550     1000    2           0
...
```
![Example Data](./doc/example_data.png)

> **NOTE:** For conformation with the standard, the start and end positions are 1-based, and the end position is inclusive. 
> However, for the sake of sanity of the author, these are converted to 0-based, and the end position is exclusive. See details on usage: TODO: link.

## BED File Input

BED files do not have a sample identifier, so the `sample_id` column is not present. If you aim to process just a single sample, you can format it for input using the following command:

```bash
awk 'BEGIN{FS="[ \t]+";OFS="\t"} {print "sample1", $1, $2, $3}' yourfile.txt | sed '1isample_id\tchrom\tstart\tend' > modified_file.tsv
```

# CLI

The command line interface uses the following pattern:

`cnsistent [command] cna_file_path [options]`

The following commands are available (see details below)

* `fill`: Adds NaN regions to the CNA data to match the assembly.
* `impute`: Imputes missing values in the CNA data.
* `coverage`: Calculates coverage for filled (but not imputed) CNA data.
* `ploidy`: Calculates aneuploidy for CNA data (NaNs are ignored).
* `bin`: Creates bins for CNA data.
* `cluster`: Creates a clustering of breakpoints.

The `cna_file_path` must point to a CNA file as described in the [Input](#input) section.

The following additional optional arguments are shared:

* `--samples string`: Path to the file with labelling of the samples. The file should be a TSV with the following columns: `sample_id, sex`. This is an optional argument. 
    * If `sex` is not provided, it is derived from the presence of `chrY` in the dataset. 
    * If `--samples` is used, only the samples in the file are processed.
    * If the samples file contains additional columns, and the program outputs a `samples` file (e.g. `coverage`), these columns are copied to the output file.
* `--out string`: Output file. This is an optional argument. If not provided, it defaults to `./cns.out.tsv`.
* `--assembly [hg19, hg38]`: Assembly to use.  If not provided, it defaults to `hg19`.
* `--threads int`: Number of threads to use, defaults to `1`. (Note: `cluster` is not parallelizeable).	
* `--verbose`: If provided, progress will be printed.


> For examples of using the CLI, you can refer to  the file `example_CLI.sh`.

## `fill`

Fills any gaps in *CNS* file with Nan values. The following steps are performed: 

1. Added NaN segments to the telomeres.
2. Fill gaps in the data with NaN values.
3. Add missing chromosomes, if they are missing compared to the reference.
4. Merge neighbouring segments with the same copy numbers (or NaNs). Both minor and major must match.

## `impute`

Replaces any NaNs in the *CNS* file with the values of the closest neighbouring region that is not NaN. The following steps are performed:

1. Assign telomers the values of the closest neighbouring region is not NaN.
2. Split the gaps and to each side, assign the values of the closes neighbouring region that is not NaN, in the direction from the center towards the side (see example below). 
3. If a whole chromosome is missing, or declared as NaN, its assigned to 0 for its whole length.
4. Merge neighbouring segments with the same copy numbers (or NaNs). Both minor and major CN values must match to be merged.

| ![Example Imputation](./doc/example_impute.png) | ![Example Merging](./doc/example_merge.png) | 
| -- | -- |


## `coverage`

Calculates the coverage of the *CNS* file. The coverage is calculated as the fraction of the genome that has a CN value assigned.

> coverage should ber run on a filled, but **not** imputed dataset.

The following statistics are calculated and stored in a *samples* file:

* `sex`: `xy` for male, `xx` for female. If this information is not specified, `xy` is used if and only if `chrY` is present in the sample.
* `chrom_count`: the number of autosomes that had any CN values assigned
* `chrom_missing`: the list of chromosomes that have no CN values assigned
* `bases_aut`: total number of bases found to have an assigned copy number on autosomes
* `bases_sex`: total number of bases found to have an assigned copy number on sex chromosomes
* `bases_tot`: sum of the two above
* `frac_aut`: the number of autosomal bases that had an assigned CN value divided by the total number of autosomal bases
* `frac_sex`: the number of sex bases that had an assigned CN value divided by the total number
* `frac_tot`: the number of chromosomes that have no CN values assigned

## `ploidy`

Calculates the portions of the genome that are aneuploid, or for absent in case of male sex chromosomes.

> ploidy should ber run on an imputed dataset.

The following statistics are calculated and stored in a *samples* file:

* `breaks_aut`: the number of autosomal breakpoints, a healthy genome should have 0
* `breaks_sex`: the number of breakpoints on sex chromosomes
* `breaks_tot`: the sum of the two above
* `ane_major_cn_aut`: the fraction of major copy number on autosomes that is not 1	
* `ane_minor_cn_aut`: the fraction of minor copy number on autosomes that is not 1
* `ane_total_cn_aut`: the fraction of total copy number on autosomes that is not 2 - note that this does not prevent e.g. 0 + 2 from being counted as diploid
* `ane_major_cn_frac_aut`, `ane_minor_cn_frac_aut`, `ane_total_cn_frac_aut`: the fraction of the genome that has a CN different from 1, 1, and 2, respectively, on autosomes.
* `ane_major_cn_sex`, `ane_minor_cn_sex`, `ane_total_cn_sex`, `ane_major_cn_frac_sex`, `ane_minor_cn_frac_sex`, `ane_total_cn_frac_sex`: the same as above, but for the sex chromosomes

## `bin`

Conducts the binning - for selected segments, the aggregate CN value for selected samples are calculated. 

> binning should ber run on an imputed dataset.

Binning can be done on the whole genome, or on selected segments. Additionally, segments can be removed from the dataset before binning. The following steps are performed:

1. If `--select` is provided, only the selected segments are used for binning. The segments are selected based on the `chrom`, `start`, and `end` columns. The segments can be selected by chromosome, chromosome arm, or chromosome band. If select is not provided, the whole genome is used, based on the assembly.
    * If `--filter` is provided, segments that are strictly smaller than the value are removed.
2. If `--remove` is provided, these segments are subtracted from the selection. The segments are removed based on the `chrom`, `start`, and `end` columns. The segments can be removed by chromosome, chromosome arm, or chromosome band. 
    * If `--filter` is provided, segments that are strictly smaller than the value are removed both before and after the subtraction process, i.e. a if a remove segment is smaller than the filter value, it is not used in subtraction. If the subtraction results in a segments smaller than the filter, it is likewise not used for binning.
4. If `--bins` is provided, the data is binned into segments of the given size. The segments are created by aggregating the CN values of the selected segments. 
    * The aggregation can be done using one of the following `--aggregate` functions: `mean`, `min`, `max`. For `mean`, the aggregation is weighted by the segment length.
5. The binned data is stored in a TSV file with the following columns: `sample_id, chrom, start, end, major_cn, minor_cn`. The names of the columns can be different, but the order must be the same.
    * If `--onlybins` is provided, only the segemnents are created, with the columns `chrom, start, end`.


![Example Binning](./doc/example_bins.png)

Additional arguments:

* `--select string`: Path to the selection file or one of [`arms`, `bands`]. The file should be a TSV with the following columns: `chrom, start, end`. This is an optional argument. Default is whole chromosomes.
* `--remove string`: Path to the removal file or `gaps`. The file should be a TSV with the following columns: `chrom, start, end`. This is an optional argument. Default is no removal.
* `--bins int`: Size of the bins. This is an optional argument. If not provided, the default value is `0` for no bins.
* `--filter int`: Minimum size of the segment. This is an optional argument. If not provided, the default value is `0` for no filtering.
* `--aggregate [mean, min, max]`: Aggregation function to use. This is an optional argument. If not provided, the default value is `mean`.
* `--onlybins`: If provided, only the bins are created, without the sample_id column.


An example of use:

```bash	
# Create bins of 1Mb size, exclude gaps, and filter out small segments
python cns.py bin cns.tsv --bins 200000 --remove gaps --filter 
# Sample copy numbers of genes, take the minimum CN value across the gene to exclude fractional segments
python cns.py bin cns.tsv --select genes.tsv --aggregate min
```

## `cluster`

Calculates consistent regions using breakpoints clustering.

Additional arguments:

* `--dist`: Maximum distance between breakpoints to be considered a part of the same segment. If not provided, the default value is `0`, which will create the exact segements of the dataset.

> NOTE: If first or last breakpoint is less than `--dist` from the end of the chromosome, the segment is extended to the end of the chromosome.

# API (Python library)

The library is split into three blocks:
* `cns.process`: The methods for processing the data, including filling, imputing, and calculating coverage.
* `cns.analyze`: The methods for deriving sample metrics from the CNA data and plotting.
* `cns.utils`: Miscellaneous utilities, mainly for converstion and file IO.

> For examples of using the Library, you can check the file `example_API.py`.

# REFERENCE

## Data

TCGA data obtained from ASCATv3 at: https://github.com/VanLoo-lab/ascat/tree/master/ReleasedData    
Cite: https://www.pnas.org/doi/full/10.1073/pnas.1009843107    

PCAWG data obtained from: https://dcc.icgc.org/releases/PCAWG/consensus_cnv    
Cite: https://www.nature.com/articles/s41587-019-0055-9    

TRACERx data obtained from: https://zenodo.org/records/7649257    
Cite: https://www.nature.com/articles/s41586-023-05729-x

COSMIC cancer set obtained from: https://cancer.sanger.ac.uk/census   
Cite: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6450507

Human genome gene set obtained using PyENSEMBL (2023).
Cite: https://academic.oup.com/nar/article/51/D1/D933/6786199

Cytoband, Gap data obtained from: https://genome.ucsc.edu
Cite: https://www.nature.com/articles/35057062

## Please cite

...

## License

Copyright © 2023 Dr. Adam Streck, adam.streck@gmail.com

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the “Software”), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
