.. _CLI_usage:

Tool usage
==========

The command line interface uses the following pattern:

``cns [command] cns_file_path [options]``

The following commands are available (see details below):

* :ref:`fill_cmd`: Adds missing segments to the CNS data to match the assembly (fills gaps with NaNs).
* :ref:`impute_cmd`: Imputes missing values in the CNS data.
* :ref:`coverage_cmd`: Calculates coverage for filled (but not imputed) CNS data.
* :ref:`ploidy_cmd`: Calculates aneuploidy for CNS data (NaNs are ignored).
* :ref:`breakage_cmd`: Creates a clustering of breakpoints.
* :ref:`segment_cmd`: Create a consistent segmentation.
* :ref:`aggregate_cmd`: Creates bins for CNS data.

The ``cns_file_path`` must point to a CNS file as described in the following section.

.. _cli_data:

Data
----

The tool expects an unprocessed copy number dataset in the form of a ``TSV`` file with the following column scheme: ``sample_id, chrom, start, end, [*_cn]``.

* The ``sample_id`` is the identifier of the sample,
* ``chrom`` is the name of the chromosome,
* ``start`` and ``end`` are the start and end positions of the segment,
* ``[*_cn]`` is typically one or two copy number segments.

E.g.:

.. csv-table:: Raw CN data for two samples.

    sample_id, chrom, start, end, CN1, CN2
    s1, chr19, 1000000, 3000000, 1,
    s1, chr19, 3000000, 12000000, 1, 1
    s1, chr19, 12000000, 14000000, , 1
    s1, chr19, 14000000, 21000000, 3, 1
    s1, chr19, 21000000, 25000000, 3, 
    s1, chr19, 28000000, 58500000, 3,
    s2, chr19, 1000000, 24000000, 2,
    s2, chr19, 29000000, 58000000, 0,

.. _cns_raw_image:

.. figure:: ../cns_raw.png
    :width: 500px

    Raw copy number data for each sample and allele.

.. note::

    For conformation with the standard, the start and end positions are 1-based, and the end position is inclusive.
    However, for the sake of sanity of the author, these are converted to 0-based, and the end position is exclusive. See details on usage: TODO: link.

BED File Input
--------------

BED files do not have a sample identifier, so the ``sample_id`` column is not present. If you aim to process just a single sample, you can format it for input using the following command:

.. code-block:: bash

    awk 'BEGIN{FS="[ \t]+";OFS="\t"} {print "sample1", $1, $2, $3}' yourfile.txt | sed '1isample_id\tchrom\tstart\tend' > modified_file.tsv

Commands
--------

.. _fill_cmd:

``fill``
--------

Fills any gaps in *CNS* file with Nan values. The following steps are performed:

1. Added NaN segments to the telomeres.
2. Fill gaps in the data with NaN values.
3. Add missing chromosomes, if they are missing compared to the reference.
4. Merge neighbouring segments with the same copy numbers (or NaNs). Both minor and major must match.

.. _impute_cmd:

``impute``
----------

Replaces any NaNs in the *CNS* file with the values of the closest neighbouring region that is not NaN. The following steps are performed:

1. Assign telomeres the values of the closest neighbouring region is not NaN.
2. Split the gaps and to each side, assign the values of the closes neighbouring region that is not NaN, in the direction from the center towards the side (see example below).
3. If a whole chromosome is missing, or declared as NaN, its assigned to 0 for its whole length.
4. Merge neighbouring segments with the same copy numbers (or NaNs). Both minor and major CN values must match to be merged.

.. image:: ../cns_imputed.png
   :width: 640px

.. _coverage_cmd:

``coverage``
------------

.. note::

    For all sample statistics, the values are calculated for autosomes, sex chromosomes, and the total genome, with the values being suffixed with ``_aut``, ``_sex``, ``_tot``, respectively. If sex chromosomes are missing from data altogether, only ``_aut`` values are calculated.

Calculates the coverage of the *CNS* file. The coverage is calculated as the fraction of the genome that has a CN value assigned.

.. note::

    coverage should be run on a filled, but **not** imputed dataset.

The following statistics are calculated and stored in a *samples* file:

* ``sex``: ``xy`` for male, ``xx`` for female. If this information is not specified, ``xy`` is used if and only if ``chrY`` is present in the sample.
* ``chrom_count``: the number of autosomes that had any CN values assigned
* ``chrom_missing``: the list of chromosomes that have no CN values assigned
* ``bases_*``: total number of bases with CNS values assigned. ``*`` is for ``aut``, ``sex``, ``tot``, referring to autosomes, sex chromosomes, and the sum of both, respectively
* ``frac_*``: fraction of bases with CNS values assigned over the total number of bases

.. _ploidy_cmd:

``ploidy``
----------

Calculates the portions of the genome that are aneuploid, or for absent in case of male sex chromosomes.

.. note::

    ploidy should be run on an imputed dataset.

The following statistics are calculated and stored in a *samples* file:

* ``breaks_*``: the number of breakpoints, a healthy genome should have 0
* ``ane_+_*``: the number of bases that have a CN different from normal (so 2, or 1 for male sex chromosomes), ``+`` expands into CN columns names
* ``ane_+_frac_*``: the fraction of aneuploid bases over the total number of bases

.. _breakage_cmd:

``breakage``
------------

TODO

.. _segment_cmd:

``segment``
-----------

Conducts the binning - for selected segments, the aggregate CN value for selected samples are calculated.

.. note::

    binning should be run on an imputed dataset.

Binning can be done on the whole genome, or on selected segments. Additionally, segments can be removed from the dataset before binning. The following steps are performed:

1. If ``--select`` is provided, only the selected segments are used for binning. The segments are selected based on the ``chrom``, ``start``, and ``end`` columns. The segments can be selected by chromosome, chromosome arm, or chromosome band. If select is not provided, the whole genome is used, based on the assembly.
    1.1 If ``--filter`` is provided, segments that are strictly smaller than the value are removed.
2. If ``--remove`` is provided, these segments are subtracted from the selection. The segments are removed based on the ``chrom``, ``start``, and ``end`` columns. The segments can be removed by chromosome, chromosome arm, or chromosome band.
    2.1 If ``--filter`` is provided, segments that are strictly smaller than the value are removed both before and after the subtraction process, i.e. a if a remove segment is smaller than the filter value, it is not used in subtraction. If the subtraction results in a segments smaller than the filter, it is likewise not used for binning.
3. If ``--bins`` is provided, the data is binned into segments of the given size. The segments are created by aggregating the CN values of the selected segments.
    3.1 The aggregation can be done using one of the following ``--aggregate`` functions: ``mean``, ``min``, ``max``. For ``mean``, the aggregation is weighted by the segment length.
4. The binned data is stored in a TSV file with the following columns: ``sample_id, chrom, start, end, major_cn, minor_cn``. The names of the columns can be different, but the order must be the same.
    4.1 If ``--segfile`` is provided, only the segments are created, with the columns ``chrom, start, end``.

.. image:: ../cns_segmented.png
   :width: 640px

.. _aggregate_cmd:

``aggregate``
-------------

TODO

.. code-block:: bash

    cns aggregate cns.tsv --segments segments.tsv

.. image:: ../cns_aggregated.png
   :width: 640px