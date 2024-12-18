"""
This module provides functions for processing Copy Number Segment (CNS) data. It includes functions to fill gaps, impute missing values, aggregate data, calculate coverage, and more.

Functions
---------
- main_fill: Fills gaps in the CNS data and adds missing chromosomes.
- main_impute: Imputes missing values in the CNS data.
- main_fill_imp: Fills gaps, adds missing chromosomes, and imputes missing values in the CNS data.
- main_breakage: Identifies breakpoints in CNS data.
- main_coverage: Calculates coverage statistics for CNS data.
- main_ploidy: Calculates ploidy statistics for CNS data.
- main_segment: Segments CNS data based on specified parameters.
- main_aggregate: Aggregates CNS data over specified genomic segments.
- main_seg_agg: Segments CNS data and aggregates the results.
"""

from .analyze import *
from .process import *
from .utils import *


def main_fill(cns_df, samples_df=None, cn_columns=None, assembly=hg19, add_missing_chromosomes=True, print_info=False):
    """
    Fills gaps in the CNS data and adds missing chromosomes.

    Parameters
    ----------
    cns_df : pandas.DataFrame
        DataFrame containing CNS (Copy Number Segment) data.
    samples_df : pandas.DataFrame, optional
        DataFrame containing sample information. If None, samples are created from `cns_df`.
    cn_columns : list of str, optional
        List of column names for copy number data. If None, columns are inferred from `cns_df`.
    assembly : object, optional
        Genome assembly to use. Default is `hg19`.
    add_missing_chromosomes : bool, optional
        If True, adds missing chromosomes to the data. Default is True.
    print_info : bool, optional
        If True, prints informational messages during processing. Default is False.

    Returns
    -------
    pandas.DataFrame
        DataFrame with filled gaps and added missing chromosomes.

    Notes
    -----
    This function performs the following steps:
    1. Adds tails to the CNS data to cover chromosome ends.
    2. Fills gaps between segments.
    3. Optionally adds missing chromosomes.
    4. Removes outlier segments.
    5. Merges neighboring segments with the same copy number.

    Examples
    --------
    >>> filled_cns = main_fill(cns_df)
    """
    if not isinstance(cns_df, pd.DataFrame):       
        raise ValueError(f"cns_df must be a DataFrame, got {type(cns_df)}") 
    if samples_df is None:
        log_info(print_info, "No samples provided, creating samples from CNS data.")
        samples_df = samples_df_from_cns_df(cns_df)
    elif not isinstance(samples_df, pd.DataFrame):
        raise ValueError(f"samples_df must be a DataFrame, got {type(samples_df)}")
    cn_columns = get_cn_cols(cns_df, cn_columns)
    cns_tailed_df = add_tails(cns_df, assembly, print_info=print_info)
    cns_filled_df = fill_gaps(cns_tailed_df, print_info=print_info)
    if add_missing_chromosomes:
        cns_filled_df = add_missing(cns_filled_df, samples_df, assembly, print_info=print_info)
    cns_cleared_df = remove_outliers(cns_filled_df, assembly, print_info=print_info)
    res_df = merge_cns_df(cns_cleared_df, cn_columns, print_info=print_info)
    return res_df


def main_impute(cns_df, samples_df=None, method="extend", cn_columns=None, print_info=False):
    """
    Imputes missing values in the CNS data.

    Parameters
    ----------
    cns_df : pandas.DataFrame
        DataFrame containing CNS data.
    samples_df : pandas.DataFrame, optional
        DataFrame containing sample information. Required if `method` is "diploid".
    method : str, optional
        Imputation method to use. Options are "extend" or "diploid". Default is "extend".
    cn_columns : list of str, optional
        List of column names for copy number data. If None, columns are inferred from `cns_df`.
    print_info : bool, optional
        If True, prints informational messages during processing. Default is False.

    Returns
    -------
    pandas.DataFrame
        DataFrame with imputed copy number values.

    Notes
    -----
    This function performs the following steps:
    1. Imputes missing values based on the specified method.
    2. Fills any remaining NaNs with zeros.
    3. Merges neighboring segments with the same copy number.

    Examples
    --------
    >>> imputed_cns = main_impute(cns_df, method="diploid")
    """
    if not isinstance(cns_df, pd.DataFrame):       
        raise ValueError(f"cns_df must be a DataFrame, got {type(cns_df)}") 
    cn_columns = get_cn_cols(cns_df, cn_columns)
    if samples_df is None:
        if method == "diploid":
            log_info(print_info, "Diploid imputation method requires samples_df, but none provided, creating samples from CNS data.")
            samples_df = samples_df_from_cns_df(cns_df)    
    elif not isinstance(samples_df, pd.DataFrame):
        raise ValueError(f"samples_df must be a DataFrame, got {type(samples_df)}")
    imputed_df = cns_impute(cns_df, samples_df, method, cn_columns=cn_columns, print_info=print_info)
    filled_df = fill_nans_with_zeros(imputed_df, cn_columns=cn_columns, print_info=print_info)
    res_df = merge_cns_df(filled_df, cn_columns=cn_columns, print_info=print_info)
    return res_df


def main_fill_imp(
    cns_df,
    samples_df=None,
    cn_columns=None,
    assembly=hg19,
    add_missing_chromosomes=True,
    method="extend",
    print_info=False,
):
    """
    Fills gaps in the CNS data, adds missing chromosomes, and imputes missing values.

    Parameters
    ----------
    cns_df : pandas.DataFrame
        DataFrame containing CNS (Copy Number Segment) data.
    samples_df : pandas.DataFrame, optional
        DataFrame containing sample information. If None, samples are created from `cns_df`.
    cn_columns : list of str, optional
        List of column names for copy number data. If None, columns are inferred from `cns_df`.
    assembly : object, optional
        Genome assembly to use. Default is `hg19`.
    add_missing_chromosomes : bool, optional
        If True, adds missing chromosomes to the data. Default is True.
    method : str, optional
        Imputation method to use. Options are "extend" or "diploid". Default is "extend".
    print_info : bool, optional
        If True, prints informational messages during processing. Default is False.

    Returns
    -------
    pandas.DataFrame
        DataFrame with filled gaps, added missing chromosomes, and imputed values.
    """
    res_df = main_fill(cns_df, samples_df, cn_columns, assembly, add_missing_chromosomes, print_info)
    res_df = main_impute(res_df, samples_df, method, cn_columns, print_info)
    return res_df


# any: if True, based is considered as covered if any CN column has values assigned
def main_coverage(cns_df, samples_df=None, cn_columns=None, segs=None, assembly=hg19, print_info=False):
    """
    Calculates coverage statistics for CNS data.

    Parameters
    ----------
    cns_df : pandas.DataFrame
        DataFrame containing CNS data.
    samples_df : pandas.DataFrame, optional
        DataFrame containing sample information. If None, samples are created from `cns_df`.
    cn_columns : list of str, optional
        List of column names for copy number data. If None, columns are inferred from `cns_df`.
    segs : segments dictionary, optional
        Dictionary of segments used for selective masking. Default is None.
    assembly : Assembly object, optional
        Genome assembly to use. Default is `hg19`.
    print_info : bool, optional
        If True, prints informational messages during processing. Default is False.

    Returns
    -------
    pandas.DataFrame
        DataFrame containing coverage statistics for each sample.

    Notes
    -----
    The function calculates coverage metrics such as the fraction of the genome covered by CNS data.

    Examples
    --------
    >>> coverage_stats = main_coverage(cns_df)
    """
    if not isinstance(cns_df, pd.DataFrame):       
        raise ValueError(f"cns_df must be a DataFrame, got {type(cns_df)}") 
    if samples_df is None:
        log_info(print_info, "No samples provided, creating samples from CNS data.")
        samples_df = samples_df_from_cns_df(cns_df)    
    elif not isinstance(samples_df, pd.DataFrame):
        raise ValueError(f"samples_df must be a DataFrame, got {type(samples_df)}")
    
    res_df = samples_df.copy()
    cn_columns = get_cn_cols(cns_df, cn_columns)

    if segs is not None:
        cns_df = aggregate_by_segments(cns_df, segs, "none", cn_columns, print_info)
    norm_sizes = get_norm_sizes(segs, assembly)

    # Select the rows where copy-numbers are not Not a Number (NaN == NaN) is false
    any_nan_df = cn_not_nan(cns_df, cn_columns, True)
    if len(cn_columns) == 2:
        both_nan_df = cn_not_nan(cns_df, cn_columns, False)

    res_df = get_missing_chroms(any_nan_df, res_df, segs, assembly)

    res_df = get_covered_bases(any_nan_df, res_df, True)
    res_df = normalize_feature(res_df, "cover_any", norm_sizes)

    if len(cn_columns) == 2:
        res_df = get_covered_bases(both_nan_df, res_df, False)
        res_df = normalize_feature(res_df, "cover_both", norm_sizes)
    return res_df


def main_breakage(cns_df, samples_df=None, cn_columns=None, segs=None, assembly=hg19, print_info=False):
    """
    Identifies breakpoints in CNS data.

    Parameters
    ----------
    cns_df : pandas.DataFrame
        DataFrame containing CNS data.
    threshold : float, optional
        Threshold for detecting breakpoints. Default is 0.5.
    cn_columns : list of str, optional
        List of column names for copy number data. If None, columns are inferred from `cns_df`.        
    segs : segments dictionary, optional
        Dictionary of segments used for selective masking. Default is None.    
    assembly : Assembly object, optional
        Genome assembly to use. Default is `hg19`.
    print_info : bool, optional
        If True, prints informational messages during processing. Default is False.

    Returns
    -------
    pandas.DataFrame
        DataFrame containing breakpoint information.

    Notes
    -----
    This function detects breakpoints in the CNS data based on changes in copy number values.

    Examples
    --------
    >>> breakpoints = main_breakage(cns_df, threshold=1.0)
    """
    if not isinstance(cns_df, pd.DataFrame):       
        raise ValueError(f"cns_df must be a DataFrame, got {type(cns_df)}") 
    if samples_df is None:
        log_info(print_info, "No samples provided, creating samples from CNS data.")
        samples_df = samples_df_from_cns_df(cns_df)    
    elif not isinstance(samples_df, pd.DataFrame):
        raise ValueError(f"samples_df must be a DataFrame, got {type(samples_df)}")
    
    res_df = samples_df.copy()
    cn_columns = get_cn_cols(cns_df, cn_columns)
    if segs is not None:
        cns_df = aggregate_by_segments(cns_df, segs, "none", cn_columns, print_info)

    # check if non of the cn_columns are NaN
    if cns_df[cn_columns].isna().any().any():
        raise RuntimeError("Cannot calculate breakage with NaN values in CN columns, impute first.")

    total_added = False
    if len(cn_columns) == 2:
        cns_df["total_cn"] = cns_df[cn_columns].sum(axis=1)
        cn_columns.append("total_cn")
        total_added = True

    for cn_col in cn_columns:
        cns_subset_df = cns_df[["sample_id", "chrom", "start", "end", cn_col]]
        segs_df = merge_cns_df(cns_subset_df, cn_col, False)
        res_df = calc_breaks_per_sample(segs_df, res_df, cn_col, assembly)
        res_df = calc_step_per_sample(segs_df, res_df, cn_col, assembly)

    if total_added:
        cns_df.drop(columns=["total_cn"], inplace=True)

    return res_df


def main_ploidy(cns_df, samples_df=None, cn_columns=None, segs=None, assembly=hg19, print_info=False):
    """
    Calculates ploidy statistics for CNS data.

    Parameters
    ----------
    cns_df : pandas.DataFrame
        DataFrame containing CNS data.
    samples_df : pandas.DataFrame, optional
        DataFrame containing sample information. If None, samples are created from `cns_df`.
    cn_columns : list of str, optional
        List of column names for copy number data. If None, columns are inferred from `cns_df`.   
    segs : segments dictionary, optional
        Dictionary of segments used for selective masking. Default is None.    
    assembly : Assembly object, optional
        Genome assembly to use. Default is `hg19`.
    print_info : bool, optional
        If True, prints informational messages during processing. Default is False.

    Returns
    -------
    pandas.DataFrame
        DataFrame containing ploidy statistics for each sample.

    Notes
    -----
    This function calculates ploidy metrics such as the fraction of the genome that is aneuploid.

    Examples
    --------
    >>> ploidy_stats = main_ploidy(cns_df)
    """
    if not isinstance(cns_df, pd.DataFrame):       
        raise ValueError(f"cns_df must be a DataFrame, got {type(cns_df)}") 
    cn_columns = get_cn_cols(cns_df, cn_columns)

    if samples_df is None:
        log_info(print_info, "No samples provided, creating samples from CNS data.")
        samples_df = samples_df_from_cns_df(cns_df)    
    elif not isinstance(samples_df, pd.DataFrame):
        raise ValueError(f"samples_df must be a DataFrame, got {type(samples_df)}")
    res_df = samples_df.copy()
    
    if segs is not None:
        log_info(print_info, "Aggregating CN data by provided segments.")
        cns_df = aggregate_by_segments(cns_df, segs, "none", cn_columns, print_info)

    norm_sizes = get_norm_sizes(segs, assembly)

    if cns_df[cn_columns].isna().any().any():
        log_warn("NaNs are not considered in ploidy calculations, it is recommended to impute first.")
        cns_df = cns_df[cns_df[cn_columns].notna().all(axis=1)]

    log_info(print_info, "Calculating LOH for each sample.")
    res_df = calc_loh_bases(res_df, cns_df, cn_columns, "both", assembly)
    res_df = normalize_feature(res_df, "loh_both", norm_sizes)
    res_df = calc_loh_bases(res_df, cns_df, cn_columns, "any", assembly)
    res_df = normalize_feature(res_df, "loh_any", norm_sizes)
    log_info(print_info, "Calculating aneuploidy for each sample.")
    res_df = calc_ane_bases(res_df, cns_df, cn_columns, "any", assembly)
    res_df = normalize_feature(res_df, "ane_any", norm_sizes)
    if len(cn_columns) == 2:
        res_df = calc_ane_bases(res_df, cns_df, cn_columns, "both", assembly)
        res_df = normalize_feature(res_df, "ane_both", norm_sizes)
        log_info(print_info, "Calculating imbalance for each sample.")
        for col_i in range(2):
            res_df = calc_imb_bases(cns_df, res_df, cn_columns, col_index=col_i, assembly=assembly)
            res_df = normalize_feature(res_df, f"imb_{cn_columns[col_i]}", norm_sizes)

    log_info(print_info, "Calculating ploidy for each sample.")
    for cn_col in cn_columns:
        res_df[f"ploidy_{cn_col}"] = calc_ploidy_per_column(cns_df, cn_col)
    return res_df


def main_segment(
    segment_source="whole",
    remove_segs=None,
    split_size=-1,
    cluster_dist=-1,
    filter_size=-1,
    assembly=hg19,
    print_info=False,
):
    """
    Creates a segmentation based on specific segments.

    Parameters
    ----------
    segment_source : a cns_df, or a segments dictionary, or one of of ["whole", "arms", "bands", "centromeres"]
        What to create the segmentation based on. If a CNS DataFrame is provided, unique segments are inferred from it. If a built-in type is provided, segments are created based on the type.
    remove_segs : segments dictionary, optional
        DataFrame containing segments to remove from the selection.
    split_size : int, optional
        Size in base pairs to split segments. Default is -1 (no splitting).
    cluster_dist : int, optional
        Distance in base pairs to merge nearby segments. Default is -1 (no merging).
    filter_size : int, optional
        Minimum size in base pairs to filter segments. Default is -1 (no filtering).
    assembly : Assembly object, optional
        Genome assembly to use. Default is `hg19`.
    print_info : bool, optional
        If True, prints informational messages during processing. Default is False.

    Returns
    -------
    list of tuples
        List of segments after processing.

    Examples
    --------
    >>> segmented_cns = main_segment(cns_df)  # consistent segmentation
    """
    # if input data is a DataFrame, convert it to unique segments
    if isinstance(segment_source, str) and segment_source in ["whole", "arms", "bands"]:            
        log_info(print_info, f"Creating {segment_source} segments...")
        segment_source = regions_select(segment_source, assembly)
    elif isinstance(segment_source, pd.DataFrame):        
        input_segs = cns_df_to_segments(segment_source)
        input_breaks = segments_to_breaks(input_segs)
        segment_source = breaks_to_segments(input_breaks)
    elif not isinstance(segment_source, dict):
        raise ValueError(f"input_data must be a CNS DataFrame, a segments dictionary or one of of ['whole', 'arms', 'bands'], got {type(segment_source)}")
    res = process_segments(segment_source, remove_segs, filter_size)
    if cluster_dist > 0:
        res = cluster_segments(res, cluster_dist, True, print_info)
    if split_size > 0:
        res = split_segments(res, split_size)
    return res


def main_aggregate(cns_df, segs, how="mean", cn_columns=None, print_info=False):
    """
    Aggregates CNS data over specified genomic segments.

    Parameters
    ----------
    cns_df : pandas.DataFrame
        DataFrame containing CNS data.
    segs : pandas.DataFrame
        DataFrame containing segments over which to aggregate CNS data.
    how : str, optional
        Aggregation method. Options are "mean", "min", "max", or "none". Default is "mean".
    cn_columns : list of str, optional
        List of column names for copy number data. If None, columns are inferred from `cns_df`.
    print_info : bool, optional
        If True, prints informational messages during processing. Default is False.

    Returns
    -------
    pandas.DataFrame
        DataFrame with aggregated copy number values.

    Notes
    -----
    If `how` is not "none" and there are NaNs in `cns_df`, a warning is issued because NaNs are not considered in aggregation.

    Examples
    --------
    >>> aggregated_cns = main_aggregate(cns_df, segs, how="max")
    """
    if not isinstance(cns_df, pd.DataFrame):       
        raise ValueError(f"cns_df must be a DataFrame, got {type(cns_df)}") 
    cn_columns = get_cn_cols(cns_df, cn_columns)
    if how not in ["", "none"] and cns_df[cn_columns].isna().any().any():
        log_warn("NaNs found, it is recommended to impute first.")
    return aggregate_by_segments(cns_df, segs, how, cn_columns, print_info)


def main_seg_agg(
    cns_df,
    segment_source="whole",
    cn_columns=None,
    remove_segs=None,
    how="mean",
    split_size=-1,
    cluster_dist=-1,
    filter_size=-1,
    assembly=hg19,
    print_info=False,
):
    """
    Segments CNS data and aggregates the results.

    Parameters
    ----------
    cns_df : pandas.DataFrame
        DataFrame containing CNS (Copy Number Segment) data.
    segment_source : a cns_df (possibly the same as cns_df), a segments dictionary or one of of ["whole", "arms", "bands", "centromeres"]
        What to create the segmentation based on. If a CNS DataFrame is provided, unique segments are inferred from it. If a built-in type is provided, segments are created based on the type.
    cn_columns : list of str, optional
        List of column names for copy number data. If None, columns are inferred from `cns_df`.
    remove_segs : segments dictionary, optional
        DataFrame containing segments to remove from the selection.
    how : str, optional
        Aggregation method to use. Default is "mean".
    split_size : int, optional
        Size in base pairs to split segments. Default is -1 (no splitting).
    cluster_dist : int, optional
        Distance in base pairs to merge nearby segments. Default is -1 (no merging).
    filter_size : int, optional
        Minimum size in base pairs to filter segments. Default is -1 (no filtering).
    assembly : Assembly object, optional
        Genome assembly to use. Default is `hg19`.
    print_info : bool, optional
        If True, prints informational messages during processing. Default is False.

    Returns
    -------
    pandas.DataFrame
        DataFrame with aggregated CNS data.
    """
    segs = main_segment(segment_source, remove_segs, split_size, cluster_dist, filter_size, assembly, print_info)
    res_df = main_aggregate(cns_df, segs, how, cn_columns, print_info)
    return res_df


