import pandas as pd
from cns.analyze.aneuploidy import *
from cns.analyze.coverage import *
from cns.analyze.breakpoints import *
from cns.process.aggregation import *
from cns.process.breakpoints import *
from cns.process.cluster import *
from cns.process.imputation import *
from cns.process.normalize import get_norm_sizes, normalize_feature
from cns.process.segments import *
from cns.utils.canonization import *
from cns.utils.conversions import *
from cns.utils.files import *
from cns.utils.logging import *
from cns.utils.assemblies import *


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
    res_df = merge_neighbours(cns_cleared_df, cn_columns, print_info=print_info)
    return res_df


def main_impute(cns_df, samples_df=None, method="extend", cn_columns=None, print_info=False):
    cn_columns = find_cn_cols_if_none(cns_df, cn_columns)
    if method == "diploid" and samples_df is None:
        log_info(print_info, "Diploid imputation requires samples, but none provided, creating samples from CNS data.")
        samples_df = samples_df_from_cns_df(cns_df)
    imputed_df = cns_impute(cns_df, samples_df, method, cn_columns=cn_columns, print_info=print_info)
    filled_df = fill_nans_with_zeros(imputed_df, cn_columns=cn_columns, print_info=print_info)
    res_df = merge_neighbours(filled_df, cn_columns=cn_columns, print_info=print_info)
    return res_df


def main_aggregate(cns_df, segs, how="mean", cn_columns=None, print_info=False):
    cn_columns = find_cn_cols_if_none(cns_df, cn_columns)
    if how != "" and how != "none" and cns_df[cn_columns].isna().any().any():
        log_warn("NaNs are not considered in aggregation calculations, it is recommended to impute first.")
    return aggregate_by_segments(cns_df, segs, how, cn_columns, print_info)


# any: if True, based is considered as covered if any CN column has values assigned
def main_coverage(cns_df, samples_df=None, cn_columns=None, segs=None, assembly=hg19, print_info=False):
    if samples_df is None:
        log_info(print_info, "No samples provided, creating samples from CNS data.")
        samples_df = samples_df_from_cns_df(cns_df)
    res_df = samples_df.copy()
    cn_columns = find_cn_cols_if_none(cns_df, cn_columns)

    if segs is not None:
        cns_df = aggregate_by_segments(cns_df, segs, "none", cn_columns, print_info)
    norm_sizes = get_norm_sizes(segs, assembly)

    # Select the rows where copy-numbers are not Not a Number (NaN == NaN) is false
    het_nan_df = get_not_nan(cns_df, cn_columns, True)
    if len(cn_columns) == 2:        
        hom_nan_df = get_not_nan(cns_df, cn_columns, False)    
    
    res_df = get_missing_chroms(het_nan_df, res_df, segs, assembly)
    
    res_df = get_covered_bases(het_nan_df, res_df, True)
    res_df = normalize_feature(res_df, "cover_het", norm_sizes)

    if len(cn_columns) == 2:
        res_df = get_covered_bases(hom_nan_df, res_df, False)
        res_df = normalize_feature(res_df, "cover_hom", norm_sizes)
    return res_df


def main_signatures(cns_df, samples_df=None, cn_columns=None, segs=None, assembly=hg19, print_info=False):
    if samples_df is None:
        log_info(print_info, "No samples provided, creating samples from CNS data.")
        samples_df = samples_df_from_cns_df(cns_df)
    res_df = samples_df.copy()
    cn_columns = find_cn_cols_if_none(cns_df, cn_columns)
    if segs is not None:
        cns_df = aggregate_by_segments(cns_df, segs, "none", cn_columns, print_info)

    # check if non of the cn_columns are NaN
    if cns_df[cn_columns].isna().any().any():
        raise RuntimeError("Cannot calculate signatures with NaN values in CN columns, impute first.")

    if len(cn_columns) == 2:
        cns_df["total_cn"] = cns_df[cn_columns].sum(axis=1)
        cn_columns.append("total_cn")

    for cn_col in cn_columns:
        segs_df = prepare_segments(cns_df, cn_col)
        res_df = calc_breaks_per_sample(segs_df, res_df, cn_col, assembly)
        res_df = calc_step_per_sample(segs_df, res_df, cn_col, assembly)
    return res_df


# NaNs are no
def main_ploidy(cns_df, samples_df=None, cn_columns=None, segs=None, assembly=hg19, print_info=False):
    if samples_df is None:
        log_info(print_info, "No samples provided, creating samples from CNS data.")
        samples_df = samples_df_from_cns_df(cns_df)
    res_df = samples_df.copy()
    cn_columns = find_cn_cols_if_none(cns_df, cn_columns)
    if segs is not None:
        cns_df = aggregate_by_segments(cns_df, segs, "none", cn_columns, print_info)    
    norm_sizes = get_norm_sizes(segs, assembly)

    if cns_df[cn_columns].isna().any().any():
        log_warn("NaNs are not considered in ploidy calculations, it is recommended to impute first.")
        cns_df = cns_df[cns_df[cn_columns].notna().all(axis=1)]

    res_df = calc_ane_bases(res_df, cns_df, cn_columns, "het", assembly)
    res_df = normalize_feature(res_df, "ane_het", norm_sizes)
    res_df = calc_loh_bases(res_df, cns_df, cn_columns, "hom", assembly)
    res_df = normalize_feature(res_df, "loh_hom", norm_sizes)
    res_df = calc_loh_bases(res_df, cns_df, cn_columns, "het", assembly)    
    res_df = normalize_feature(res_df, "loh_het", norm_sizes)
    if len(cn_columns) == 2:
        res_df = calc_ane_bases(res_df, cns_df, cn_columns, "hom", assembly)
        res_df = normalize_feature(res_df, "ane_hom", norm_sizes)
        for col_i in range(2):
            res_df = calc_imb_bases(cns_df, res_df, cn_columns, col_index=col_i, assembly=assembly)
            res_df = normalize_feature(res_df, f"imb_{cn_columns[col_i]}", norm_sizes)
    return res_df


def main_segment(
    cns_df, select_segs=None, remove_segs=None, split_size=0, merge_dist=0, filter_size=0, assembly=hg19, print_info=False
):
    if select_segs == None:
        select_segs = genome_to_segments(assembly)
    res = get_genome_segments(select_segs, remove_segs, filter_size)
    if merge_dist > 0:
        breaks = get_breaks_from_cns(cns_df, keep_ends=False, assembly=assembly)
        res = cluster_within_segments(breaks, res, merge_dist, print_info)
    if split_size > 0:
        res = split_segments(res, split_size)
    return res
