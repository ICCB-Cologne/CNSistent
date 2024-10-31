.. _cli_quickstart:

CLI Quickstart
==============

.. code-block:: bash

   1. Create the conda environment: `conda env create -f ./cnsistent.yml`
   2. Activate the environment: `conda activate cns`
   3. Install the package from location: `pip install -e .`
   4. Process data: `bash ./scripts/data_process.sh`
   5. Aggregate processed data: `bash ./scripts/data_aggregate.sh`
   6. Run example analysis: `python ./example_API.py`

NOTES:

* By default, 16 threads are used, if that causes problems (crashes), reduce the number of threads in the `data_process.sh` and `data_aggregate.sh` scripts.
* The `example_API.py` is split into cells that can be run individually in an IDE.
* You can also install the package with `pip install .`, however there is a set of utility functions for loading data in `cns.data_utils.py` that will not be accessible then.