![Example Data](./docs/Logo.png)

CNSistent is a Python tool for processing and analyzing copy number data. It is designed to work with data from a variety of sources, including TCGA, PCAWG, TRACERx, and COSMIC. The tool is designed to be easy to use, and to provide a comprehensive set of analyses and visualizations.

Full documentation including examples is available [here](./docs/build/html/index.html).

# Example

The below data are directly taken from the published [TRACERx Zenodo archive](https://zenodo.org/records/7649257).

## 1. Load CNS Data and Display Heatmap
Load CNS data from a CSV file and visualize the first 5 rows using a heatmap.

```python
import cns
raw_df = cns.load_cns("./data/20220803_TxPri_mphase_by_sample_df.reduced.csv", cn_columns=["nMajor", "nMinor"], sep=",", print_info=True)
cns.fig_heatmap(cns.cns_head(raw_df, 5), max_cn=6)
```

![Raw Data Heatmap](./docs/files/intro_1.png)

## 2. Impute Missing Segments
Fill in missing segments in the data, impute using the extension method, and display a heatmap for the first 5 rows.

```python
imp_df = cns.main_fill_imp(raw_df, print_info=True)
cns.fig_heatmap(cns.cns_head(imp_df, 5), max_cn=6)
```
![Imputed Data Heatmap](./docs/files/intro_2.png)

## 3. Create 3 mb Segments and convert to a feature array
Aggregate the imputed CNS data into 3 MB segments and convert it into a feature array.

```python
seg_df = cns.main_seg_agg(imp_df, split_size=3_000_000, print_info=True)
features, rows, columns = cns.bins_to_features(seg_df)
print("(alleles, samples, bins):", features.shape)
```
`(alleles, samples, bins): (2, 403, 960)`

## 4. Group Segments by Cancer Type
Group the CNS data by cancer type, calculate the total CN, and visualize mean linear profiles.

```python
sample_df = cns.load_samples("./data/20221109_TRACERx421_all_patient_df.tsv")
type_groups = {c: cns.select_cns_by_type(seg_df, sample_df, c, "histology_multi_full") for c in ["LUAD", "LUSC"]}
groups_df = cns.stack_groups([cns.group_samples(v, group_name=k) for k, v in type_groups.items()])
cns.fig_lines(cns.add_total_cn(groups_df), cn_columns="total_cn")
```
![Raw Data Heatmap](./docs/files/intro_3.png)

The example code is also in `example_API.py`.

# Usage

CNSistent is available both as a [stand-alone](#cli) tool and as a [Python library](#library).

Requires:

* Git LFS
* Python 3.8+
* Pip 21.3+

> Git Clone the will also download all source datasets, in total around ~1GB.

## Quickstart (from Bash)

1. Create the conda environment: `conda env create -f ./cnsistent.yml`
2. Activate the environment: `conda activate cns`
3. Install the package from location: `pip install -e .`
4. Process data: `bash ./scripts/data_process.sh`

NOTES:

* By default, 16 threads are used, if that causes problems (crashes), reduce the number of threads in the `data_process.sh` and `data_aggregate.sh` scripts.
* The `example_API.py` is split into cells that can be run individually in an IDE.
* You can also install the package with `pip install .`, however there is a set of utility functions for loading data in `cns.data_utils.py` that will not be accesible then.


# API (Python library)

## Quickstart

The repository contains a set of utility functions for loading the available data in `cns.data_utils`. 

* `main_load`: The entry point function. It loads and filters samples from the datasets (PCAWG, TCGA, TRACERx) together with the imputed CNS data. If bins value is provided, instead of CNS data, binned data is loaded for the perspective bin size.
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

### File accecss:
* `out_path`: Returns the default output path where output of processing/binning files are stored.
* `data_path`: Returns the default data path where the input files (bundled together with the repository) are stored.
* `load_cns_out`: Loads CNS file from the output folder.
* `load_samples_out`: Loads samples file from the output folder.
* `load_bins`: Loads bins from the output folder.

## Library

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

Breakpoints are stored in a dictionary mapping a set of breakpoints to a chromosome. Boundaries are normally excluded, e.g.
```
{
    'chr1' : [125000000],
    ...
}
```
would be the arm breakpoints for hg19 chromosome 1.

#### Segments

Segments are analogous to breakpoints, instead of position defining a range. The range is start-inclusive, end-exclusive. E.g.

```
{
    'chr1' : [(0, 125000000), (125000000, 249250621)],
    ...
}
```
would be the arm breakpoints for hg19 chromosome 1.

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

Functions for working with segments. Segments are dictionaries of tuples `{chr: [(start, end), ...], ...}`, where the start is inclusive, and the end is exclusive.

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

#### Aggregation

Aggregation will produce segments of a certain size, aggregating the copy number values of the segment chunks into a single segment.

There are the following aggregate functions: `mean`, `min`, `max`, and `none`. The `none` function will just split exiting bins, without additional aggregation. This is useful if you want to introduce additional breakpoints into the data.

Aggregation can be done either using explicit segments, explicit breakpoints, or a breakpoint type (e.g. `arms`, `1000000`).


### Analyze

The analyze module calculates statistics for the CNS data.

* `coverage`: Calculates the proportion of genome with assigned (not NaN) CN values.
* `ploidy`: Calculates the proportion of genome with aneuploid CN values (different from 2 or 1 for male sex chroms).
* `breakage`: Calculates the signatures related statistics - currently it only calculates breakpoints per sample/chromosome.

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

#### Filtering (cutoff)

To help with defining cutoff regions for samples, it is possible to use a utility function that finds a knee/elbow point for a dataset. 

The primary function is `find_bends`, which will convert a value (e.g. coverage) into a cumulative distribution between `min_val` and `max_val` with `steps` steps. The function will then find the knee/elbow point in the distribution.

The function will return the point where the slope of the cumulative distribution has the highest convex (knee) or concave (elbow) curvature. To avoid finding local minimum, `dist` can be set to consider the angle between the mean of the nearest `dist` points. If `allow_pad` is set to true, the endpoinst are also considered as potential knees/elbows, with the slope at the beginning and the end being 0.

![CN Tracks](./docs/example_knee.png)

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

The data are obtained and processed using the code in `./notebooks/data/data_obtain.ipynb`. Sources:

TCGA data obtained from ASCATv3 at: https://github.com/VanLoo-lab/ascat/tree/master/ReleasedData    
Cite: https://www.pnas.org/doi/full/10.1073/pnas.1009843107   
The results published here are in part based upon data generated by the TCGA Research Network: https://www.cancer.gov/tcga.  

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

TBD, pre-print expected summmer 2024.

### The MIT License

Copyright © 2023 Dr. Adam Streck, adam.streck@gmail.com

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the “Software”), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.