![CNSistent Logo](https://cnsistent.readthedocs.io/en/latest/_images/Logo.png)


# <div align="center"> [**~ READ THE DOCS HERE ~**](https://cnsistent.readthedocs.io/en/latest/)</span>


CNSistent is a Python tool for processing and analyzing copy number data. It is designed to work with data from a variety of sources. The tool is designed to be easy to use, and to provide a comprehensive set of analyses and visualizations.

> This repository also contains data from the PCAWG, TCGA, and TRACERx datasets, as well as gene sets from COSMIC and Ensembl, therefore git clone will download around 1GB of data. If you only want to use the library, download via PyPI [TODO].


## Repository Data Quickstart

This repository contains raw data from PCAWG, TCGA, TRACERx, as well as genomic locations. The data needs to be processed first before it can be used.

### Requirements

* Git LFS
* Python 3.8+
* Pip 21.3+
* (Optional) Conda for environment creation

### Processing

1. Clone the repository: `git clone https://bitbucket.org/schwarzlab/cnsistent`
2. Install dependencies (`pip install -r requirements.txt`) or create a Conda environment (`conda env create -f cnsistent.yml`).
3. Install the package from location: `pip install -e .`
The `-e` will make sure that the data files can be accessed under the `cns` package.
4. [Optional] Process data: `bash ./scripts/data_process.sh` - will create imputed data and sample statistics.
5. [Optional] Aggregate data: `bash ./scripts/data_aggregate.sh` - will aggregate the imputed data using 15 different segmentation strategies.

### Usage

To load the data use:

``` python
from cns.data_utils import main_load
samples_df, cns_df = main_load("imp")
```

This will load the imputed and filtered data for all datasets. 

The ``samples_df`` and ``cns_df`` are Pandas dataframes. 
The former contains information about each samples as well as its statistics (e.g. ``ane_both_ all`` for homozygous aneuploidy across all chromosomes).
The latter contains the copy number segments for each sample in the form of ``sample_id``, ``chrom``, ``start``, ``end``, ``major_cn``, ``minor_cn``, ``name`` where ``name`` identifies each segment. 
For example to load CNs for the COSMIC genes, data you can use the same function:

``` python

    samples_df, cns_df = main_load("COSMIC")
    cns_df.head()
```

would produce

``` python

      sample_id chrom start    end     major_cn  minor_cn name
    0 SP101724  chr1  2160133  2241558 2         2        SKI
    1 SP101724  chr1  2487077  2496821 2         2        TNFRSF14
    2 SP101724  chr1  2985731  3355185 2         2        PRDM16
    3 SP101724  chr1  6241328  6269449 2         2        RPL22
    4 SP101724  chr1  6845383  7829766 2         2        CAMTA1
```

Alternativelly you can call:

* `main_load` to only load samples,
* `main_load("raw")` to load the raw data,
* `main_load("imp")` to load the imputed data,
* `main_load(agg_type)` to load the aggregated bins, if the aggregation has been done, which can be one of: `["1MB", "2MB", "3MB", "5MB", "10MB", "250KB", "500KB", "whole", "arms", "bands", "COSMIC", "ENSEMBL"]`.

### Notes

* By default, 16 threads are used, if that causes problems (crashes), reduce the number of threads in the `data_process.sh` and `data_aggregate.sh` scripts.
* The `example_API.py` is split into cells that can be run individually in an IDE.
* You can also install the package with `pip install .`, however there is a set of utility functions for loading data in `cns.data_utils.py` that will not be accesible then.
* Conda is optional, you can also install required packages manually using PIP based on the list in [cnsistent.yml](./cnsistent.yml).
* Additionally, 5 of the PCAWG medulloblastoma samples have been labeled as female in the source, however they contained CN calls for chromosome Y and we have therefore re-labelled them as male. 

## Repository Structure
  
**`.`**

*  `cnsistent.yml`: Conda environment file for the CNSistent package, references `requirements.txt`.
*  `requirements.txt`: Packages required to run the CNSistent package.
*  `example_API.py`: Example code for using the CNSistent package. 
*  `example_CLI.sh`: Example code for using the CNSistent package from the command line.
*  `pyproject.toml`: Configuration for packaging tools.
  
**`cns/`**

Contains the main code for the CNSistent package.

**`data/`**

Contains the raw data from PCAWG, TCGA, TRACERx, as well as genomic locations, also a notebook used to obtain them or merge source files.

**`docs/`**

Contains the documentation for the CNSistent package. The documentation is built using Sphinx, with the source in the `./docs/source` folder. The documentation can be built using the `make html` command in the `./docs` folder, provided the requirements in `./docs/requirements.txt` are met.

**`notebooks/`**

Contains notebooks used for data processing and analysis:

*  `analyze_break_clusters.ipynb`: A notebook used to analyze the breakpoint clustering, based on the distance between merged breakpoints. 
*  `analyze_CN_clipping.ipynb`: Evaluation of result of clipping the CN segment values, in particular the effects on distribution and proportion that is clipped of. 
*  `analyze_coverage.ipynb`: Calculates the proportion of the genome that is covered by segments and locations where it applies. 
*  `analyze_features.ipynb`	: Calculates and plots features across datasets. 
*  `analyze_lung.ipynb` : Plots the lung cancer data across datasets and cancer types, in particular for chromosome 3 and genes that have been established as important by IG method. 
*  `analyze_SOX2_overlay.ipynb`: Plots the SOX2 gene overlay on the lung cancer data. 
*  `analyze_types.ipynb`: Plots the distribution of cancer types and overall CN across datasets. 
*  `data_obtain.ipynb`: A notebook that has been used to obtain the raw data and potentially merge files where needed. 
*  `docs_illustrations.ipynb`: A notebook used to create illustrations for the documentation. 
*  `docs_knee_detection.ipynb`: A demo of the kneepoint detection algorithm. 
*  `docs_runtime.ipynb`: Calculates the runtime of the data processing across 1-32 threads (log scale). 

**`scripts/`**

*  `data_process.sh`: Fills and imputes the raw data. Also calculates the data stats, in particular coverage and aneuploidy. 
*  `data_aggregate.sh`: Creates various segmentations, and aggregates the preprocessed data based on these segmentations. 
Depends on `data_cluster.py` for breakpoint clustering.
*  `data_time.sh`: Run time tests for the data processing across `1-32` threads (log scale). 

**`tests/`**

*  `in` and `out`: Contains the input and output data for the tests. Output is generated using `example_CLI.sh`.  
*  `test_cli.sh` : Executes the tests and outputs to `./tests/temp`.  
*  `test_*`: unittest based tests of the public API. 

## Reference

### Data

The contents of the data folder were obtained by processing the following sources, accessed in December 2023.

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

TBD, pre-print expected 2024.

### The MIT License

Copyright © 2023 Dr. Adam Streck, adam.streck@gmail.com

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the “Software”), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.