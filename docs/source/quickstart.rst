.. _quickstart:

Quickstart
==========

Processing the bundled data
---------------------------

The `cnsistent` repository contains data from the PCAWG, TCGA, and TRACERx studies, which can be preprocessed and segmented using bundled scripts under the following steps:

1. Start a ``bash``-compatible shell.
2. Create the conda environment: ``conda env create -f ./cnsistent.yml``
3. Activate the environment: ``conda activate cns``
4. Install the package from location: ``pip install -e .``
5. Process data: ``bash ./scripts/data_process.sh``
6. Optional: Aggregate processed data: ``bash ./scripts/data_aggregate.sh``

**NOTES:**

* The aggregation step will take several hours and produce ~40GB of data.
* By default, 16 threads are used, if that causes problems (crashes), reduce the number of threads in the ``data_process.sh`` and ``data_aggregate.sh`` scripts.
* You can also install the package with ``pip install``, however there is a set of utility functions for loading data in ``cns.data_utils.py``` that will not be accessible then.

Accessing the bundled data
--------------------------

The above produced data can be accessed using the ``cns.data_utils`` module. For example, to load the imputed data:

.. code-block:: python

    from cns.data_utils import main_lod
    samples_df, cns_df = main_load("imp")

The ``samples_df`` and ``cns_df`` are Pandas dataframes. 
The former contains information about each samples as well as its statistics (e.g. ``ane_hom_all`` for homozygous aneuploidy across all chromosomes).
The latter contains the copy number segments for each sample in the form of ``sample_id``, ``chrom``, ``start``, ``end``, ``major_cn``, ``minor_cn``, ``name`` where ``name`` identifies each segment. 
For example to load CNs for the COSMIC genes, data you can use the same function:

.. code-block:: python

    samples_df, cns_df = main_load("COSMIC")
    cns_df.head()

would produce

.. code-block:: python

      sample_id chrom start    end     major_cn  minor_cn name
    0 SP101724  chr1  2160133  2241558 2         2        SKI
    1 SP101724  chr1  2487077  2496821 2         2        TNFRSF14
    2 SP101724  chr1  2985731  3355185 2         2        PRDM16
    3 SP101724  chr1  6241328  6269449 2         2        RPL22
    4 SP101724  chr1  6845383  7829766 2         2        CAMTA1

Basic tool usage
----------------

CNSistent reads SCNA profiles as ``.tsv`` files. Have an example file ``data.tsv``:

.. code-block:: tsv

    sample_id    chrom   start   end     total_cn
    sample1      chr1    100     200     1       
    ...

.. note::
    Column naming is fully describe in the :ref:`Column Naming` section.

To preprocess the segments:

.. code-block:: bash

    cns fill data.tsv --out filled.tsv
    cns impute filled.tsv --out imputed.csv

To create statistics:

.. code-block:: bash

    cns coverage data.tsv --out samples.tsv
    cns ploidy imputed.tsv --samples samples.tsv --out samples.tsv
    cns signatures imputed.tsv --samples samples.tsv --out samples.tsv

To calculate the mean ploidy per chromosome arm:

.. code-block:: bash

    cns segment arms --out arms.bed
    cns aggregate imputed.tsv --segments --out a_bins.tsv

To conduct breakpoint clustering with 1 mb distance:

.. code-block:: bash

    cns cluster imputed.tsv --merge 1000000 --out clust.bed
    cns aggregate imputed.tsv  --segments clust.bed --out c_bins.tsv

To conduct segmentation using 5 mb bins:

.. code-block:: bash

    cns segment whole --step 5000000 --out clust.bed
    cns aggregate data.tsv  --segments clust.bed --out c_bins.tsv


Basic library usage
-------------------
TODO