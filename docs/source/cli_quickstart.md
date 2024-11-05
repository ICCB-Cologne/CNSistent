Tool usage
==========

The command line interface uses the following pattern:

`cns [command] cns_file_path [options]`

The following commands are available (see details below)

* `fill`: Adds NaN regions to the CNS data to match the assembly.
* `impute`: Imputes missing values in the CNS data.
* `coverage`: Calculates coverage for filled (but not imputed) CNS data.
* `ploidy`: Calculates aneuploidy for CNS data (NaNs are ignored).
* `bin`: Creates bins for CNS data.
* `cluster`: Creates a clustering of breakpoints.

The `cns_file_path` must point to a CNS file as described in the [Data](#Data) section.

## Data

The tool expects an unprocessed copy number dataset in the form of a `TSV` file with the following column scheme: `sample_id, chrom, start, end, [*_cn]`. 

* The `sample_id` is the identifier of the sample, 
* `chrom` is the name of the chromosome, 
* `start` and `end` are the start and end positions of the segment, 
* `[*_cn]` is typically one or two copy number segments.

E.g.:

```tsv
sample_id   	chrom   start   end     major_cn    minor_cn
sample1         chr1    100     200     1           0
sample1         chr1    550     1000    2           0
...
```

| ![Example Data](../example_data.png) | -- | 
| -- | -- |


> **NOTE:** For conformation with the standard, the start and end positions are 1-based, and the end position is inclusive. 
> However, for the sake of sanity of the author, these are converted to 0-based, and the end position is exclusive. See details on usage: TODO: link.

### BED File Input

BED files do not have a sample identifier, so the `sample_id` column is not present. If you aim to process just a single sample, you can format it for input using the following command:

```bash
awk 'BEGIN{FS="[ \t]+";OFS="\t"} {print "sample1", $1, $2, $3}' yourfile.txt | sed '1isample_id\tchrom\tstart\tend' > modified_file.tsv
```
## Commands

The following additional optional arguments are shared:

* `--cols int`: Number of copy number columns. If specified the columns on position `5:5+cols` are used as copy numbers. If not specified, the parser will consider columns that start or end with `"CN"` or `"cn"`.
* `--samples string`: Path to the file with labelling of the samples. The file should be a TSV with the following columns: `sample_id, sex`. This is an optional argument. 
    * If `sex` is not provided, it is derived from the presence of `chrY` in the dataset. 
    * If `--samples` is used, only the samples in the file are processed.
    * If the samples file contains additional columns, and the program outputs a `samples` file (e.g. `coverage`), these columns are copied to the output file.
* `--out string`: Output file. This is an optional argument. If not provided, it defaults to `./cns.out.tsv`.
* `--noheader`: If provided, the header is not expected in the input file.
* `--nosample`: If provided, the sample_id column is not expected in the input file.
* `--assembly [hg19, hg38]`: Assembly to use.  If not provided, it defaults to `hg19`.
* `--threads int`: Number of threads to use, defaults to `1`. (Note: `cluster` is not parallelizable).	
* `--subsplit int`: Number of segments to split the data into for parallel processing. If not provided, the default value is `1` for no splitting. This helps to limit memory consumption for large datasets.
* `--verbose`: If provided, progress will be printed.

> For examples of using the CLI, you can refer to  the file `example_CLI.sh`.

### `fill`

Fills any gaps in *CNS* file with Nan values. The following steps are performed: 

1. Added NaN segments to the telomeres.
2. Fill gaps in the data with NaN values.
3. Add missing chromosomes, if they are missing compared to the reference.
4. Merge neighbouring segments with the same copy numbers (or NaNs). Both minor and major must match.

### `impute`

Replaces any NaNs in the *CNS* file with the values of the closest neighbouring region that is not NaN. The following steps are performed:

1. Assign telomeres the values of the closest neighbouring region is not NaN.
2. Split the gaps and to each side, assign the values of the closes neighbouring region that is not NaN, in the direction from the center towards the side (see example below). 
3. If a whole chromosome is missing, or declared as NaN, its assigned to 0 for its whole length.
4. Merge neighbouring segments with the same copy numbers (or NaNs). Both minor and major CN values must match to be merged.

| ![Example Imputation](../example_impute.png) | ![Example Merging](../example_merge.png) | 
| -- | -- |


### `coverage`

> For all sample statistics, the values are calculated for autosomes, sex chromosomes, and the total genome, with the values being suffixed with `_aut`, `_sex`, `_tot`, respectively. If sex chromosomes are missing from data altogether, only `_aut` values are calculated.

Calculates the coverage of the *CNS* file. The coverage is calculated as the fraction of the genome that has a CN value assigned.

> coverage should ber run on a filled, but **not** imputed dataset.

The following statistics are calculated and stored in a *samples* file:

* `sex`: `xy` for male, `xx` for female. If this information is not specified, `xy` is used if and only if `chrY` is present in the sample.
* `chrom_count`: the number of autosomes that had any CN values assigned
* `chrom_missing`: the list of chromosomes that have no CN values assigned
* `bases_*`: total number of bases with CNS values assigned. `*` is for `aut`, `sex`, `tot`, referring to autosomes, sex chromosomes, and the sum of both, respectively
* `frac_*`: fraction of bases with CNS values assigned over the total number of bases

### `ploidy`

Calculates the portions of the genome that are aneuploid, or for absent in case of male sex chromosomes.

> ploidy should ber run on an imputed dataset.

The following statistics are calculated and stored in a *samples* file:

* `breaks_*`: the number of breakpoints, a healthy genome should have 0
* `ane_+_*`: the number of bases that have a CN different from normal (so 2, or 1 for male sex chromosomes), `+` expands into CN columns names
* `ane_+_frac_*`: the fraction of aneuploid bases over the total number of bases

### `signatures`

TODO

### `segment`

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
    * If `--segfile` is provided, only the segments are created, with the columns `chrom, start, end`.


![Example Binning](../example_bins.png)

Additional arguments:

* `--select string`: Path to the selection file or one of [`arms`, `bands`]. The file should be a TSV with the following columns: `chrom, start, end`. This is an optional argument. Default is whole chromosomes.
* `--remove string`: Path to the removal file or `gaps`. The file should be a TSV with the following columns: `chrom, start, end`. This is an optional argument. Default is no removal.
* `--split int`: Size of the bins. This is an optional argument. If not provided, the default value is `0` for no bins.
* `--merge`: Maximum distance between breakpoints to be considered a part of the same segment. If not provided, the default value is `0`, which will create the exact segments of the dataset.
* `--filter int`: Minimum size of the segment. This is an optional argument. If not provided, the default value is `0` for no filtering.
* `--aggregate [mean, min, max]`: Aggregation function to use. This is an optional argument. If not provided, the default value is `mean`.


An example of use:

```bash	
## Create bins of 1Mb size, exclude gaps, and filter out small segments
cns segment cns.tsv --split 200000 --remove gaps --filter 
## Sample copy numbers of genes, take the minimum CN value across the gene to exclude fractional segments
cns segment cns.tsv --select genes.tsv --aggregate min
```

### `aggregate`

TODO

```bash	
cns aggregate cns.tsv --segments segments.tsv
```