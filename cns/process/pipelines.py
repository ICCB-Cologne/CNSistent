import pandas as pd
import numpy as np


from cns.analyze.aneuploidy import calc_ane_bases, calc_imb_bases, calc_loh_bases
from cns.analyze.coverage import normalize_feature, get_covered_bases, get_missing_chroms, get_not_nan
from cns.analyze.breakpoints import calc_breaks_per_sample, calc_seg_size_per_sample, calc_step_per_sample, prepare_segments
from cns.process.binning import bin_by_segments
from cns.process.breakpoints import get_breaks
from cns.process.cluster import calc_clusters, cluster_within_segments
from cns.process.imputation import add_missing, add_tails, cns_impute, fill_gaps, fill_nans_with_zeros, merge_neighbours, remove_outliers
from cns.process.segments import get_genome_segments, regions_remove, regions_select, split_segments
from cns.utils.canonization import find_cn_cols_if_none, is_hap_spec, rename_cn_cols
from cns.utils.files import samples_df_from_cns_df
from cns.utils.logging import log_info
from cns.utils.assemblies import hg19


def main_fill(cns_df, samples_df=None, cn_columns=None, assembly=hg19, add_missing_chromosomes=True, print_info=False):
    if samples_df is None:
        log_info(print_info, "No samples provided, creating samples from CNS data.")
        samples_df = samples_df_from_cns_df(cns_df)
    cn_columns = find_cn_cols_if_none(cns_df, cn_columns)
    cns_tailed_df = add_tails(cns_df, assembly.chr_lens, print_info=print_info)
    cns_filled_df = fill_gaps(cns_tailed_df, print_info=print_info)
    if add_missing_chromosomes:
        cns_filled_df = add_missing(cns_filled_df, samples_df, assembly.chr_lens, print_info=print_info)
    cns_cleared_df = remove_outliers(cns_filled_df, assembly.chr_lens, print_info=print_info)
    res = merge_neighbours(cns_cleared_df, cn_columns, print_info=print_info)
    return res


def main_impute(cns_df, samples_df=None, method='extend', cn_columns=None,  print_info=False):
    cn_columns = find_cn_cols_if_none(cns_df, cn_columns)
    if method == 'diploid' and samples_df is None:
        log_info(print_info, "Diploid imputation requires samples, but none provided, creating samples from CNS data.")
        samples_df = samples_df_from_cns_df(cns_df)
    imputed_df = cns_impute(cns_df, samples_df, method, cn_columns=cn_columns, print_info=print_info)
    filled_df = fill_nans_with_zeros(imputed_df, cn_columns=cn_columns, print_info=print_info)
    res = merge_neighbours(filled_df, cn_columns=cn_columns, print_info=print_info)
    return res


def main_bin(cns_df, segs, fun_type='mean', cn_columns=None, print_info=False):
    return bin_by_segments(cns_df, segs, fun_type, cn_columns, print_info)


# any: if True, based is considered as covered if any CN column has values assigned
def main_coverage(cns_df, samples_df, cn_columns=None, assembly=hg19, print_info=False):
    res = samples_df.copy()
    cn_columns = find_cn_cols_if_none(cns_df, cn_columns)    
    if "length" not in cns_df.columns:
        cns_df["length"] = cns_df["end"] - cns_df["start"]

    res = get_missing_chroms(cns_df, res, assembly)
    # Select the rows where copy-numbers are not Not a Number (NaN == NaN) is false
    
    hom_nan_df = get_not_nan(cns_df, cn_columns, False)
    res = get_covered_bases(hom_nan_df, res, False)
    res = normalize_feature(res, "cover_hom", assembly)

    if len(cn_columns) == 2:
        hom_nan_df = get_not_nan(cns_df, cn_columns, True)
        res = get_covered_bases(hom_nan_df, res, True)
        res = normalize_feature(res, "cover_het", assembly)
    return res


def main_signatures(cns_df, samples_df, cn_columns=None, assembly=hg19, print_info=False):
    res = samples_df.copy()
    cn_columns = find_cn_cols_if_none(cns_df, cn_columns)    
    if "length" not in cns_df.columns:
        cns_df["length"] = cns_df["end"] - cns_df["start"]

    if len(cn_columns) == 2:
        cns_df["total_cn"] = cns_df[cn_columns].sum(axis=1)
        cn_columns.append("total_cn")
    
    for cn_col in cn_columns:
        segs_df = prepare_segments(cns_df, cn_col)
        res = calc_breaks_per_sample(segs_df, res, cn_col, assembly)
        res = calc_step_per_sample(segs_df, res, cn_col, assembly)
        res = calc_seg_size_per_sample(segs_df, res, cn_col, assembly)
    
    return res

# NaNs are not considered for ploidy calculation
def main_ploidy(cns_df, samples_df, cn_columns=None, assembly=hg19, print_info=False):
    cn_columns = find_cn_cols_if_none(cns_df, cn_columns)
    if "length" not in cns_df.columns:
        cns_df["length"] = cns_df["end"] - cns_df["start"]
    
    samples_df = calc_ane_bases(cns_df, samples_df, cn_columns, assembly)
    samples_df = normalize_feature(samples_df, "ane_hom", assembly)
    samples_df = calc_loh_bases(cns_df, samples_df, cn_columns, assembly)
    samples_df = normalize_feature(samples_df, "loh_hom", assembly)
    if len(cn_columns) == 2:
        samples_df = normalize_feature(samples_df, "ane_het", assembly)
        samples_df = normalize_feature(samples_df, "loh_het", assembly)
        samples_df = calc_imb_bases(cns_df, samples_df, cn_columns, col_index=0, assembly=assembly)
        samples_df = normalize_feature(samples_df, f"imb_{cn_columns[0]}", assembly)
        if is_hap_spec(cn_columns):
            samples_df = calc_imb_bases(cns_df, samples_df, cn_columns, col_index=1, assembly=assembly)
            samples_df = normalize_feature(samples_df, f"imb_{cn_columns[1]}", assembly)
    return samples_df


def main_segment(cns_df, select, remove, split_size=0, merge_dist=0, filter_size=0, assembly=hg19, print_info=False):
    select = regions_select(select, assembly)
    remove = regions_remove(remove, assembly)
    segs = get_genome_segments(select, remove, filter_size)
    if merge_dist > 0:
        breaks = get_breaks(cns_df, keep_ends=False, assembly=assembly)
        return cluster_within_segments(breaks, segs, merge_dist, assembly, print_info)
    if split_size > 0:
        segs = split_segments(segs, split_size)
    segs