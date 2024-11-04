"""
Utility functions for loading and managing copy number segment (CNS) data.

This module provides functions to load, filter, and manage CNS data from different sources
(PCAWG, TRACERx, TCGA). It includes utilities for handling samples, loading binned data,
and managing file paths.
"""
import pandas as pd
from os.path import join as pjoin, abspath, dirname
from cns.utils.canonization import get_cn_cols
from cns.utils.logging import log_info
from cns.utils.selection import select_CNS_samples
from cns.utils.files import load_cns, load_samples, load_segments
from cns.utils.anomaly import find_bends, z_score_filter
import matplotlib.pyplot as plt


def get_root_path():
    """Get the root path of the CNSistent package.

    Returns
    -------
    str
        Absolute path to the package root directory.
    """
    return abspath(pjoin(dirname(__file__), ".."))


def load_cns_out(filename, raw=False):
    """Load CNS data from the output directory.

    Parameters
    ----------
    filename : str
        Name of the file to load from the output directory.
    raw : bool, optional
        If True, return the raw DataFrame without renaming columns.

    Returns
    -------
    pd.DataFrame
        DataFrame containing CNS data.
    """
    cns_df = load_cns(pjoin(out_path, filename))
    if raw:
        return cns_df
    return _rename_cns_columns(cns_df)


def load_samples_out(filename):
    """Load sample data from the output directory.

    Parameters
    ----------
    filename : str
        Name of the sample file to load.

    Returns
    -------
    pd.DataFrame
        DataFrame containing sample data.
    """
    return load_samples(pjoin(out_path, filename))


def load_all_samples(filter=True, retype=True, print_info=False):
    """Load samples from all available datasets.

    Parameters
    ----------
    filter : bool, optional
        Whether to filter samples based on coverage and aneuploidy.
    retype : bool, optional
        Whether to standardize cancer type labels.
    print_info : bool, optional
        Whether to print progress information.

    Returns
    -------
    dict
        Dictionary containing DataFrames for each dataset.
    """
    samples = {
        "PCAWG": load_samples_out("PCAWG_samples.tsv"),
        "TRACERx": load_samples_out("TRACERx_samples.tsv"),
        "TCGA_hg19": load_samples_out("TCGA_hg19_samples.tsv")
    }
    # ... rest of the function ...


def get_cns_for_type(cns, samples, type):
    """Extract CNS data for a specific cancer type.

    Parameters
    ----------
    cns : pd.DataFrame
        CNS data.
    samples : pd.DataFrame
        Sample information.
    type : str
        Cancer type to extract.

    Returns
    -------
    pd.DataFrame
        CNS data filtered for the specified cancer type.
    """
    query = f"type == '{type}'"
    ids = samples.query(query).index
    select_cns = cns.set_index("sample_id").loc[ids].reset_index()
    return select_cns


def load_merged_samples(filter=True, retype=True, print_info=False):
    """Load and merge samples from all datasets.

    Parameters
    ----------
    filter : bool, optional
        Whether to filter samples based on coverage and aneuploidy.
    retype : bool, optional
        Whether to standardize cancer type labels.
    print_info : bool, optional
        Whether to print progress information.

    Returns
    -------
    pd.DataFrame
        Combined DataFrame containing all samples.
    """
    samples = load_all_samples(filter, retype, print_info=print_info)
    for k, v in samples.items():
        v["source"] = k
    all_samples = pd.concat(samples.values())        
    log_info(print_info, f"Total samples: {len(all_samples)}")
    return all_samples


def load_merged_bins(select_samples, segment_size):
    """Load and merge binned CNS data.

    Parameters
    ----------
    select_samples : pd.DataFrame
        Samples to include.
    segment_size : int
        Size of the bins.

    Returns
    -------
    pd.DataFrame
        Combined binned CNS data.
    """
    cns = {
        "PCAWG": load_cns_out(f"PCAWG_bin_{segment_size}.tsv"),
        "TRACERx": load_cns_out(f"TRACERx_bin_{segment_size}.tsv"),
        "TCGA_hg19": load_cns_out(f"TCGA_hg19_bin_{segment_size}.tsv")
    }
    all_cns = pd.concat(cns.values())
    if select_samples is not None:
        all_cns = select_CNS_samples(all_cns, select_samples)
    return all_cns.reset_index(drop=True)


def main_load_data(segment_type=None):
    """Main function to load both samples and CNS data.

    Parameters
    ----------
    segment_type : int, optional
        If provided, load binned data with this segment size.

    Returns
    -------
    tuple
        (samples DataFrame, CNS DataFrame)
    """
    samples = load_merged_samples()
    if segment_type == None:
        cns = load_merged_cns(samples)
    else:
        cns = load_merged_bins(samples, segment_type)
    return samples, cns


def load_COSMIC():
    """Load COSMIC consensus genes.

    Returns
    -------
    pd.DataFrame
        DataFrame containing COSMIC gene data.
    """
    return load_segments(pjoin(data_path, "COSMIC_consensus_genes.bed"))


def load_ENSEMBL():
    """Load ENSEMBL coding genes.

    Returns
    -------
    pd.DataFrame
        DataFrame containing ENSEMBL gene data.
    """
    return load_segments(pjoin(data_path, "ENSEMBL_coding_genes.bed"))


def save_cns_fig(fig_name, fig=None):
    """Save a figure to the images directory.

    Parameters
    ----------
    fig_name : str
        Name of the figure file (without extension).
    fig : matplotlib.figure.Figure, optional
        Figure to save. If None, uses current figure.
    """
    if fig == None:
        fig = plt.gcf()
    fig.savefig(f"{img_path}/{fig_name}.png", bbox_inches="tight", transparent=True, dpi=300)
    fig.savefig(f"{img_path}/{fig_name}.pdf", bbox_inches="tight", transparent=True)


def save_doc_fig(fig_name, fig=None):
    """Save a figure to the documentation directory.

    Parameters
    ----------
    fig_name : str
        Name of the figure file (without extension).
    fig : matplotlib.figure.Figure, optional
        Figure to save. If None, uses current figure.
    """
    if fig == None:
        fig = plt.gcf()
    fig.savefig(f"{docs_path}/{fig_name}.png", bbox_inches="tight", transparent=True, dpi=300)
    fig.savefig(f"{docs_path}/{fig_name}.pdf", bbox_inches="tight", transparent=True)