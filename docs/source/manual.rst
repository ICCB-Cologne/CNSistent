Manual
======

Data Types
----------

There are five main data formats used within CNSistent.

``cns_df``
``````````
 A pandas DataFrame with the following columns: ``sample_id, chrom, start, end, CN, ...``. 
 The ``CN`` columns are the copy number values for each segment. The ``chrom`` column is expected to be in the format ``chr1``, ``chr2``, ..., ``chrX``, ``chrY``, ``chrM``. The ``start`` and ``end`` columns are 0-based coordinates.

.. code-block:: python

    sample_id  chrom     start       end  CN1  CN2
    0        s1  chr19         0  13000000    1    1
    1        s1  chr19  13000000  59128983    3    1
    2        s2  chr19         0  26500000    2    0
    3        s2  chr19  26500000  59128983    0    0


``samples_df``
``````````````
 A pandas DataFrame with the following columns: ``sample_id, sex``. The ``sex`` column is expected to be ``xy`` or ``xx``.

 .. code-block:: python

                sex
    sample_id	
    s1	        xx
    s2      	xy


``segments``
````````````
A dictionary of lists with keys being chromosomes and each segment being 0-indexed triples of start and end coordinates and a string id.

.. code-block:: python

    {'chr19': [(0, 13000000, 'chr19_0'), (13000000, 26500000, 'chr19_1'), (26500000, 59128983, 'chr19_2')]}

``breakpoints``
```````````````

A dictionary of lists with keys being chromosomes and each breakpoint being a 0-indexed position in the chromosome.

.. code-block:: python

    {'chr19': [0, 13000000, 26500000, 59128983]}

``Assembly``
````````````
Assembly is a class that provides information about the species. 
CNSistent currently supports ``hg19`` and ``hg38`` assemblies. If you want to use it with a different assembly, you need to create a new object.

Note that sex chromosomes are always expected to be named ``chrX`` and ``chrY``.

.. code-block:: python

    Assembly(name, chr_lens, chr_x, chr_y, gaps, cytobands)

* ``chr_lens`` is a dictionary with chromosome names as keys and lengths as values. The ``chr_x`` and ``chr_y`` are
* ``chr_x``, ``chr_y`` are the string ids for sex chromosomes. "chrX" and "chrY" are used by default.     
* ``gaps``, ``cytobands`` are segment dictionaries for gaps and cytobands, respectively. 
These can be null unless you use ``regions_select("bands")`` or ``regions_select("gaps")``.


Process
-------

Pipelines
~~~~~~~~~


Segment
~~~~~~~

Functions for working with segments. Segments are dictionaries of tuples ``{chr: [(start, end), ...], ...}``, where the start is inclusive, and the end is exclusive.

Note that you can pass longer tuples, but the result will discard the 4th and further elements.

Notable functions:

* ``merge_segments``: Will merge overlapping segments, merging is possible if ``end==start`` for two consecutive segments on the same chromosome. Note that if the segments are not sorted, you need to set ``sort=True`` to sort them first.
* ``split_segments``: Will split into equidistant chunks based on specified size (useful for binning).
* ``segment_difference``: Will remove regions from a list of segments found in another list of segments.
* ``filter_min_size``: Will remove segments strictly smaller than the specified size.

Imputation
~~~~~~~~~~

Functions for adding missing segments and values in the CNS data. The process is to first add missing regions with NaN values and then impute the missing values.

There are separate functions to fill the telomeres, fill the gaps, and add missing chromosomes.

If guessing values in imputation is not desired, the ``fill_nans_with_zeros`` function can be used to simply fill with 0 instead.

Breakpoints
~~~~~~~~~~~

Creating of breakpoints (see Breakpoint data type above). The function ``make_breaks`` will create de-novo breakpoints of a certain size, whereas ``get_breaks`` will return the breakpoints for a given ``CNS_df``.

Aggregation
~~~~~~~~~~~

Aggregation will produce segments of a certain size, aggregating the copy number values of the segment chunks into a single segment.

There are the following aggregate functions: ``mean``, ``min``, ``max``, and ``none``. The ``none`` function will just split existing bins, without additional aggregation. This is useful if you want to introduce additional breakpoints into the data.

Aggregation can be done either using explicit segments, explicit breakpoints, or a breakpoint type (e.g. ``arms``, ``1000000``).

Analyze
-------

The analyze module calculates statistics for the CNS data.

* ``coverage``: Calculates the proportion of genome with assigned (not NaN) CN values.
* ``ploidy``: Calculates the proportion of genome with aneuploid CN values (different from 2 or 1 for male sex chromosomes).
* ``breakage``: Calculates the signatures related statistics - currently it only calculates breakpoints per sample/chromosome.

Display
-------

Display functions are in three categories:

* ``fig``: A Whole figure with labels.
* ``plot``: Only plots an individual axis.
* ``labels``: Either labels or backgrounds for individual axes.

There are two main plot types:

* joint plots: ``fig_line, fig_dots, fig_bars`` - these display joint CN data from aggregates. Line plot is the primary plot, compared to a normal line plot, segments are connected only if they overlap. Dots are more practical for small regions, e.g. genes. Bars are useful for coarse data, e.g. arms.

.. image:: files/lines.png
    :alt: Lines Plot

.. image:: files/dots.png
    :alt: Dots Plot

.. image:: files/bars.png
    :alt: Bars Plot

* individual plots: ``fig_CN_tracks`` - this plot inserts all bins into a heatmap per sample, being more practical for dense data or equalized representation. Unlike the joint plots, the position on the genome is not considered and there are no gaps, so the sizes of chromosomes and the position therein are only for orientation.

.. image:: files/tracks.png
    :alt: CN Tracks

For the figures, the first parameter is always the ``CNS_df``, or a list thereof in the case of joint plots. Following optional parameters:
* column: a string describing the column to plot or a list thereof. If none is specified, all columns matching the CN column pattern are used.
* chrom: a string describing the chromosome to plot. If none is specified, all chromosomes are used.
* label: string describing the CNS_df or a list thereof. If none is specified, no label is output.
* width: width of the figure in inches, is calculated automatically if not specified.
* max_cn: The CN tracks values are cropped to avoid outliers collapsing the color scale. If not specified, ``10`` is used.

Utils
-----

Utils contain the specification for the hg19, hg38 assemblies, including the gaps and cytobands. In addition, functions for files and data are provided:


Files
~~~~~

* ``load_cns/save_cns``: Load/Save CNS data from a TSV file, with optional header and sample_id. By default, this moves from 1-based to 0-based coordinates.
* ``load_regions/save_regions``: Load/Save regions from a TSV file, reading only the ``chrom, start, end`` columns. By default, this moves from 1-based to 0-based coordinates.
* ``load_samples/save_samples``: Load samples from a TSV file. The first column is used as "sample_id" index and should match the CNS sample names.
* ``get_cn_columns``: Get the CN columns from a CNS DataFrame (start or end with ``CN`` or ``cn``).

Selection
~~~~~~~~~

Functions to select samples set from CNS df (head, tail, random), to filter chromosomes (autosomes, sex chromosomes...) and samples/CNS by type.

Conversions
~~~~~~~~~~~

Converts between CNS df, breakpoints, and segments.

Data Utils
----------

* Functions to load the datasets (PCAWG, TCGA, TRACERx) and gene sets (Ensembl, COSMIC).
* Default filtering to remove samples from the datasets (low coverage, diploid, blacklisted, ...)
* Loading binned data / processed samples.