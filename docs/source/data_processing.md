.. _data_processing:

# Data Processing

1. Start a `bash`-compatible shell.	
2. Create the conda environment: `conda env create -f ./cnsistent.yml`
3. Activate the environment: `conda activate cns`
4. Install the package from location: `pip install -e .`
5. Process data: `bash ./scripts/data_process.sh`
6. Optional: Aggregate processed data: `bash ./scripts/data_aggregate.sh`

**NOTES:**

* By default, 16 threads are used, if that causes problems (crashes), reduce the number of threads in the `data_process.sh` and `data_aggregate.sh` scripts.
* The `example_API.py` is split into cells that can be run individually in an IDE.
* You can also install the package with `pip install .`, however there is a set of utility functions for loading data in `cns.data_utils.py` that will not be accesible then.
