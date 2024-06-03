# CNSistent

CNSistent is a Python tool for processing and analyzing copy number data. It is designed to work with data from a variety of sources, including TCGA, PCAWG, TRACERx, and COSMIC. The tool is designed to be easy to use, and to provide a comprehensive set of analyses and visualizations.

## Usage

CNSistent is available both as a [stand-alone](#cli) tool and as a [Python library](#library).

Requires:
* Python 3.8+
* Pip 21.3+

### Quickstart (from Bash)

1. Install the package from location: `pip install -e .`
2. Create the conda environment: `conda env create -f ./cnsistent.yml`
3. Activate the environment: `conda activate cns`
4. Process data: `bash ./scripts/data_process.sh`
5. Bin processed data: `bash ./scripts/data_bin.sh`
6. Run example analysis: `python ./example_API.py`


NOTES:

* By default, 16 threads are used, if that causes problems (crashes), reduce the number of threads in the `data_process.sh` and `data_bin.sh` scripts.
* The `example_API.py` is split into cells that can be run individually in an IDE.
* You can also install the package with `pip install .`, however there is a set of utility functions for loading data in `cns.data_utils.py` that will not be accesible then.

## Input

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
![Example Data](./docs/example_data.png)

> **NOTE:** For conformation with the standard, the start and end positions are 1-based, and the end position is inclusive. 
> However, for the sake of sanity of the author, these are converted to 0-based, and the end position is exclusive. See details on usage: TODO: link.

### BED File Input

BED files do not have a sample identifier, so the `sample_id` column is not present. If you aim to process just a single sample, you can format it for input using the following command:

```bash
awk 'BEGIN{FS="[ \t]+";OFS="\t"} {print "sample1", $1, $2, $3}' yourfile.txt | sed '1isample_id\tchrom\tstart\tend' > modified_file.tsv
```

## CLI

The command line interface uses the following pattern:

`cns [command] cns_file_path [options]`

The following commands are available (see details below)

* `fill`: Adds NaN regions to the CNS data to match the assembly.
* `impute`: Imputes missing values in the CNS data.
* `coverage`: Calculates coverage for filled (but not imputed) CNS data.
* `ploidy`: Calculates aneuploidy for CNS data (NaNs are ignored).
* `bin`: Creates bins for CNS data.
* `cluster`: Creates a clustering of breakpoints.

The `cns_file_path` must point to a CNS file as described in the [Input](#input) section.

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

| ![Example Imputation](./docs/example_impute.png) | ![Example Merging](./docs/example_merge.png) | 
| -- | -- |


### `coverage`

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

### `bin`

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
    * If `--onlybins` is provided, only the segments are created, with the columns `chrom, start, end`.


![Example Binning](./docs/example_bins.png)

Additional arguments:

* `--select string`: Path to the selection file or one of [`arms`, `bands`]. The file should be a TSV with the following columns: `chrom, start, end`. This is an optional argument. Default is whole chromosomes.
* `--remove string`: Path to the removal file or `gaps`. The file should be a TSV with the following columns: `chrom, start, end`. This is an optional argument. Default is no removal.
* `--bins int`: Size of the bins. This is an optional argument. If not provided, the default value is `0` for no bins.
* `--filter int`: Minimum size of the segment. This is an optional argument. If not provided, the default value is `0` for no filtering.
* `--aggregate [mean, min, max]`: Aggregation function to use. This is an optional argument. If not provided, the default value is `mean`.
* `--onlybins`: If provided, only the bins are created, without the sample_id column.


An example of use:

```bash	
## Create bins of 1Mb size, exclude gaps, and filter out small segments
cns bin cns.tsv --bins 200000 --remove gaps --filter 
## Sample copy numbers of genes, take the minimum CN value across the gene to exclude fractional segments
cns bin cns.tsv --select genes.tsv --aggregate min
```

### `cluster`

Calculates consistent regions using breakpoints clustering.

Additional arguments:

* `--dist`: Maximum distance between breakpoints to be considered a part of the same segment. If not provided, the default value is `0`, which will create the exact segments of the dataset.

> NOTE: If first or last breakpoint is less than `--dist` from the end of the chromosome, the segment is extended to the end of the chromosome.

## API (Python library)

The library is split into three blocks:
* `cns.process`: The methods for processing the data, including filling, imputing, and calculating coverage.
* `cns.analyze`: Statistics and summary metrics.
* `cns.display`: Plotting of binned data.
* `cns.utils`: Conversion, selection, file IO.

> Note that many functions have assembly as an optional parameter, if not provided, the default assembly is `hg19`.

### Data types:

Functions for conversions between data types can be found in `cns.utils.conversions`.

#### CNS DataFrame

The CNS DataFrame is a pandas DataFrame with the following columns: `sample_id, chrom, start, end` followed by copy number column or columns. Typically, the columns are `major_cn, minor_cn` or `total_cn`. CNSiststen considers any column whose name starts or ends with `CN` or `cn` as a copy number column. 

```tsv
sample_id   	chrom   start   end     major_cn    minor_cn
sample1         chr1    100     200     1           0
sample1         chr1    550     1000    2           0
...
```

This is the main data type used in the library.

#### Breakpoints

Breakpoints are stored in a dictionary mapping a set of breakpoints to a chromosome. Boundaries are normally included, e.g.
```
{
    'chr1' : [0, 125000000, 249250621],
    ...
}
```
would be the arm breakpoints for hg19 chromosome 1.

#### Segments

Segments are a list of triples: `[(chrom, start, end), ...]`, used to map into chromosomes. The range is start-inclusive, end-exclusive.

#### Assemblies

Assembly a class that provides information about the species. CNSistent currently supports `hg19` and `hg38` assemblies, if you want to use it with a different assembly, you need to create a new object.

Note that sex chromosomes are always expected to be named `chrX` and `chrY`.	

### Process

#### Pipelines

The commands that are available from CMD are executed via a main function for each command, e.g. `main_impute`.

The file also contains high-level functions for region manipulation, in particular `get_genome_segments` has the following procedure:

1. Select regions from the assembly (emptystring for the whole assembly), filter those below `filter_size`.
2. If `remove` is specified, remove regions from the selection, and filter the remaining regions.
3. If `bin_size` is specified, bin the regions into equidistant segments of the given size.

#### Segment

Functions for working with segments. Segments are tuples `(chrom, start, end)`, where the start is inclusive, and the end is exclusive.

Note that you can pass longer tuples, but the result will discard the 4th and further elements.

Notable functions:

* `merge_segments`: Will merge overlapping segments, merging is possible if `end==start` for two consecutive segments on the same chromosome. Note that if the segments are not sorted, you need to set `sort=True` to sort them first.
* `split_segments`: Will split into equidistant chunks based on specified size (useful for binning).
* `segment_difference`: Will remove regions from a list of segments found in another list of segments.
* `filter_min_size`: Will remove segments strictly smaller than the specified size.

#### Imputation

Functions for adding missing segments and values in the CNS data. The process is to first add missing regions with NaN values and then impute the missing values.

There are separate functions to fill the telomeres, fill the gaps, and add missing chromosomes. 

If guessing values in imputation is not desired, the `fill_nans_with_zeros` function can be used to simply fill with 0 instead.


#### Breakpoints

Creating of breakpoints (see Breakpoint data type above). The function `make_breaks` will create de-novo breakpoints of a certain size, whereas `get_breaks` will return the breakpoints for a given `CNS_df`.

#### Binning

Binning will produce segments of a certain size, aggregating the copy number values of the segments within the bin. 

There are the following aggregate functions: `mean`, `min`, `max`, and `none`. The `none` function will just split exiting bins, without additional aggregation. This is useful if you want to introduce additional breakpoints into the data.

Binning can be done either using explicit segments, explicit breakpoints, or a breakpoint type (e.g. `arms`, `1000000`).

Three utility functions are included with binning:

* `add_cns_loc`: For a CNS df, will add the center (`mid`) of a segment within the chromosome, within the linear genome (`cum_mid`) and the length of the segment (`length`).
* `group_bins`: Will group all samples in a a CNS df, by default using the `mean` function. Segments are grouped on `cum_mid`.
* `sum_cns`: Finds columns that match the CN column pattern and sums them up into the column `total_cn`.


### Analyze

The analyze module calculates statistics for the CNS data.

* `coverage`: Calculates the proportion of genome with assigned (not NaN) CN values.
* `ploidy`: Calculates the proportion of genome with aneuploid CN values (different from 2 or 1 for male sex chroms).
* `signatures`: Calculates the signatures related statistics - currently it only calculates breakpoints per sample/chromosome.

### Display

Display functions are in three categories:

* `fig`: A Whole figure with labels.
* `plot`: Only plots an individual axis.
* `labels`: Either labels or backgrounds for individual axes.

There are two main plot types:

* joint plots: `fig_line, fig_dots, fig_bars` - these display joint CN data from aggregates. Line plot is the primary plot, compared to a normal line plot, segments are connected only if they overlap. Dots are more practical for small regions, e.g. genes. Bars are useful for coarse data, e.g. arms.

![Lines Plot](./docs/lines.png)
![Dots Plot](./docs/dots.png)
![Bars Plot](./docs/bars.png)

* individual plots: `fig_CN_tracks` - this plot inserts all bins into a heatmap per sample, being more practical for dense data or equalized representation. Unlike the joint plots, the position on the genome is not considered and there are no gaps, so the sizes of chromosomes and the position therein are only for orientation.

![CN Tracks](./docs/tracks.png)

For the figures, the first parameter is always the `CNS_df`, or a list thereof in the case of joint plots. Following optional parameters:
* column: a string describing the column to plot or a list thereof. If none is specified, all columns matching the CN column pattern are used.
* chrom: a string describing the chromosome to plot. If none is specified, all chromosomes are used.
* label: string describing the CNS_df or a list thereof. If none is specified, no label is output.
* width: width of the figure in inches, is calculated automatically if not specified.
* max_cn: The CN tracks values are cropped to avoid outliers collapsing the color scale. If not specified, `10` is used.


### Utils

Utils contain the specification for the hg19, hg38 assemblies, including the gaps and cytobands. In addition, functions for files and data are provided:

#### Files

* `load_cns/save_cns`: Load/Save CNS data from a TSV file, with optional header and sample_id. By default, this moves from 1-based to 0-based coordinates.
* `load_regions/save_regions`: Load/Save regions from a TSV file, reading only the `chrom, start, end` columns. By default, this moves from 1-based to 0-based coordinates.
* `load_samples/save_samples`: Load samples from a TSV file. The first column is used as "sample_id" index and should match the CNS sample names.
* `get_cn_columns`: Get the CN columns from a CNS DataFrame (start or end with `CN` or `cn`).


#### Selection

Functions to select samples set from CNS df (head, tail, random), to filter chromosomes (autosomes, sex chromosomes...) and samples/CNS by type.

#### Conversions

Converts between CNS df, breakpoints, and segments. 

### Data Utils

* Functions to load the datasets (PCAWG, TCGA, TRACERx) and gene sets (Ensembl, COSMIC).
* Default filtering to remove samples from the datasets (low coverage, diploid, blacklisted, ...)
* Loading binned data / processed samples.

## REFERENCE

### Data

The data are obtained and processed using the code in `./notebooks/data/data_obtain.ipynb`. The sources are as follows:

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

### Please cite

...

### The MIT License

Copyright © 2023 Dr. Adam Streck, adam.streck@gmail.com

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the “Software”), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.