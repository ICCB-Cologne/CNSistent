"""
Utility functions for loading and managing copy number segment (CNS) data.

This module provides functions to load, filter, and manage CNS data from different sources
(PCAWG, TRACERx, TCGA). It includes utilities for handling samples, loading binned data,
and managing file paths.
"""

import pandas as pd
from os.path import join as pjoin, abspath, dirname
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


img_path = pjoin(get_root_path(), "img")
out_path = pjoin(get_root_path(), "out")
data_path = pjoin(get_root_path(), "data")
docs_path = pjoin(get_root_path(), "docs/img")


def load_cns_out(filename, print_info=False):
    """Load CNS data from the output directory.

    Parameters
    ----------
    filename : str
        Name of the file to load from the output directory.
    print_info : bool, optional
        If True, print information about the loading process.

    Returns
    -------
    pd.DataFrame
        DataFrame containing CNS data.
    """
    log_info(print_info, f"Loading CNS data from {filename}")
    return load_cns(pjoin(out_path, filename))


def load_samples_out(filename, use_filter=True, print_info=False):
    """Load sample data from the output directory.

    Parameters
    ----------
    filename : str
        Name of the sample file to load.
    filter : bool, optional
        Whether to filter samples based on coverage and aneuploidy.
    print_info : bool, optional
        Whether to print progress information.

    Returns
    -------
    pd.DataFrame
        DataFrame containing sample data.
    """
    log_info(print_info, f"Loading samples from {filename}")
    samples_df = load_samples(pjoin(out_path, filename))

    if use_filter:
        # calculate bend for aneuploidy
        ane_bends = find_bends(samples_df["ane_het_aut"])
        ane_min_frac = ane_bends[0][ane_bends[2]]

        # calculate the z-score for coverage
        cover_filtered = z_score_filter(samples_df["cover_het_aut"])
        cover_min_frac = cover_filtered.min()

        samples_df = _filter_samples(samples_df, ane_min_frac, cover_min_frac, print_info)


    return samples_df


def load_segs_out(filename, print_info=False):
    """Load segment data from the output directory.

    Parameters
    ----------
    filename : str
        Name of the segment file to load.
    print_info : bool, optional
        Whether to print progress information.
    
    Returns
    -------
    pd.DataFrame
        DataFrame containing segment data.
    """
    log_info(print_info, f"Loading segments from {filename}")
    return load_segments(pjoin(out_path, filename))



def _filter_samples(samples_df, ane_min_frac=0.001, cover_min_frac=0.95, print_info=False):
    """Filter samples based on aneuploidy and coverage criteria.

    Parameters
    ----------
    samples : pd.DataFrame
        DataFrame containing sample information.
    ane_min_frac : float, optional
        Minimum fraction for aneuploidy filtering.
    cover_min_frac : float, optional
        Minimum fraction for coverage filtering.
    print_info : bool, optional
        Whether to print progress information.

    Returns
    -------
    pd.DataFrame
        Filtered DataFrame containing sample information.
    """
    log_info(print_info, f"Total samples: {len(samples_df)}")
    
    cn_neutral = samples_df.query(f"ane_het_aut < {ane_min_frac}").index
    log_info(print_info, f"{len(cn_neutral)} samples are CN neutral (below {ane_min_frac:.5f})")

    # Find samples with low coverage (below 95% in autosomes)
    low_coverage = samples_df.query(f"cover_het_aut < {cover_min_frac}").index
    log_info(print_info, f"{len(low_coverage)} samples have low coverage (below {cover_min_frac:.5f})")

    if "TCGA_type" in samples_df.columns:
        samples_df["type"] = samples_df["TCGA_type"]    
        samples_df.drop(columns=["TCGA_type"], inplace=True)
    samples_df["type"] = samples_df["type"].replace({"LUADx2": "LUAD"}).replace({"LUADx3": "LUAD"})
    untyped = samples_df[samples_df["type"].fillna('').apply(lambda x: any(not c.isupper() for c in x))].index
    log_info(print_info, f"{len(untyped)} samples do not have exact type")

    filtered_df = samples_df.query("(index not in @untyped) & (index not in @cn_neutral) & (index not in @low_coverage)")

    log_info(print_info, f"Filtered samples: {len(filtered_df)}")

    return filtered_df.copy()


def main_load(segment_type=None, use_filter=True, print_info=False):
    """
    Load samples and CNS data for PCAWG, TRACERx, and TCGA.

    The type of segments depends on the segment type. (If none is specified, only samples are loaded.)

    Possible segments for non-aggregated data are:

    * imp: Imputed segments
    * preprocess: Preprocessed segments
    * fill: Filled segments

    Possible segments for aggregated data are:

    * whole: Whole genome segments
    * arm: Arm-level segments
    * bands: Cytoband-level segments
    * COSMIC: COSMIC consensus genes
    * ENSEMBL: ENSEMBL coding genes
    * 10MB, 5MB, 3MB, 2MB, 1MB, 500KB, 250KB: Binned segments of this size
    * merge_10000000, merge_5000000, merge_250000: Clustered breakpoints of this size.

    Parameters
    ----------
    segment_type : string, optional
        If provided, load binned data with this segment type. Default is None.
    use_filter : bool, optional
        Whether to filter samples based on coverage and aneuploidy. Default is True.
    retype : bool, optional
        Whether to standardize cancer type labels. Default is True.
    print_info : bool, optional
        Whether to print progress information. Default is False.

    Returns
    -------
    samples_df : pd.DataFrame
        Combined DataFrame containing all samples.
    cns_df : pd.DataFrame
        Combined DataFrame containing all CNS data. If segment_type is None, this is None.
    """
    datasets = ["PCAWG", "TRACERx", "TCGA_hg19"]

    samples_list = []
    for dataset in datasets:
        samples = load_samples_out(f"{dataset}_samples.tsv", use_filter, print_info)
        samples["source"] = dataset
        samples_list.append(samples)
    samples_df = pd.concat(samples_list)
    log_info(print_info, f"Loaded samples: {len(samples_df)}")

    if segment_type is None:
        return samples_df, None
        
    file_type = "cns" if segment_type in ["imp", "preprocess", "fill"] else "bin"
    cns_dict = [ load_cns_out(f"{dataset}_{file_type}_{segment_type}.tsv", print_info) for dataset in datasets ]
    cns_df = pd.concat(cns_dict)
    log_info(print_info, f"Total CNS segments: {len(cns_df)}")
    cns_df = select_CNS_samples(cns_df, samples_df).reset_index(drop=True)
    log_info(print_info, f"Total CNS segments (after filtering): {len(cns_df)}")
    
    return samples_df, cns_df



def select_cns_by_type(cns, samples, type):
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


def load_COSMIC():
    return load_segments(pjoin(data_path, "COSMIC_consensus_genes.bed"))


def load_ENSEMBL():
    return load_segments(pjoin(data_path, "ENSEMBL_coding_genes.bed"))