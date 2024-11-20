.. _API:

Library usage (API)
===================

#### Assemblies

Assembly a class that provides information about the species. CNSistent currently supports `hg19` and `hg38` assemblies, if you want to use it with a different assembly, you need to create a new object.

Note that sex chromosomes are always expected to be named `chrX` and `chrY`.

## API

The repository contains a set of utility functions for loading the available data in `cns.data_utils`. 

* `main_load`: The entry point function. It loads and filters samples from the datasets (PCAWG, TCGA, TRACERx) together with the imputed CNS data. If bins value is provided, instead of CNS data, binned data is loaded for the perspective bin size.
	

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



.. toctree::
   :maxdepth: 2
   :caption: Contents:

   tutorial
   plotting