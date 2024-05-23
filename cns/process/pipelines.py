import pandas as pd
import numpy as np

from cns.analyze.aneuploidy import calc_ane_per_chrom, calc_ane_per_sample, norm_aut_aneuploidy, norm_sex_aneuploidy
from cns.analyze.coverage import get_base_frac, get_covered_bases, get_missing_chroms
from cns.analyze.signatures import add_breaks_per_sample
from cns.process.binning import add_cns_loc, bin_by_segments
from cns.process.breakpoints import calc_arm_breaks, calc_cytoband_breaks
from cns.process.cluster import created_merged_segs, get_breaks
from cns.process.imputation import add_missing, add_tails, create_imputed_entries, fill_gaps, fill_nans_with_zeros, merge_neighbours
from cns.process.segments import filter_min_size, segment_difference, split_segments
from cns.utils.conversions import genome_to_segments, breaks_to_segments, tuples_to_segments
from cns.utils.files import get_ane_cols_if_none, load_regions, samples_df_from_cns_df, find_cn_cols_if_none
from cns.utils.assemblies import hg19


def main_fill(cns_df, cn_columns=None, samples_df=None, assembly=hg19, add_missing_chromosomes=True, print_info=False):
    samples_df = samples_df_from_cns_df(cns_df)
    cn_columns = find_cn_cols_if_none(cns_df, cn_columns)

    cns_tailed_df = add_tails(cns_df, assembly.chr_lens, print_info=print_info)
    cns_filled_df = fill_gaps(cns_tailed_df, print_info=print_info)
    if add_missing_chromosomes:
        cns_filled_df = add_missing(cns_filled_df, samples_df, assembly.chr_lens, print_info=print_info)
    res = merge_neighbours(cns_filled_df, cn_columns, print_info=print_info)
    return res


def main_impute(cns_df, cn_columns=None, print_info=False):
    cn_columns = find_cn_cols_if_none(cns_df, cn_columns)

    imputed_df = create_imputed_entries(cns_df, cn_columns, print_info=print_info)
    filled_df = fill_nans_with_zeros(imputed_df, cn_columns, print_info=print_info)
    res = merge_neighbours(filled_df, cn_columns, print_info=print_info)
    assert len(res[res.isnull().any(axis=1)]) == 0, "NaNs still present in final_df."
    return res


def main_bin(cns_df, segs, fun_type='mean', print_info=False):
    return bin_by_segments(cns_df, segs, fun_type, print_info)


def main_coverage(cns_df, samples_df, cn_columns=None, assembly=hg19, any_nan=True, print_info=False):
    cn_columns = find_cn_cols_if_none(cns_df, cn_columns)
    # Select the rows where copy-numbers are not Not a Number (NaN == NaN) is false
    nan_vals = cns_df[cn_columns].isna()
    nan_filter = ~nan_vals.any(axis=1) if any_nan else ~nan_vals.all(axis=1)
    cns_vals = cns_df.loc[nan_filter].copy()
    coverage = get_missing_chroms(cns_vals, samples_df, assembly)
    coverage = get_covered_bases(cns_vals, coverage)
    coverage = get_base_frac(coverage, assembly)
    return coverage


def main_ploidy(cns_df, samples, cn_columns=None, assembly=hg19, print_info=False):
    cn_columns = find_cn_cols_if_none(cns_df, cn_columns)
    if len(cn_columns) > 2:
        raise ValueError("To calculate aneuploidy, only one (sum) or two (major, minor) CN columns are allowed.")
    elif len(cn_columns) == 2:
        cn_columns_df = cns_df[cn_columns]
        major_col = cn_columns_df.values.max(axis=1)  
        minor_col = cn_columns_df.values.min(axis=1)
        cns_df = cns_df.drop(columns=cn_columns)
        cns_df["major_cn"] = major_col
        cns_df["minor_cn"] = minor_col
        cns_df["total_cn"] = cns_df["major_cn"] + cns_df["minor_cn"]
    elif len(cn_columns) == 1:
        cns_df.rename(columns={cn_columns[0]: "total_cn"}, inplace=True)
    else:
        raise ValueError("No CN columns found.")

    samples = add_breaks_per_sample(cns_df, samples, assembly)
    cns_df = add_cns_loc(cns_df, assembly)
    pre_chr = calc_ane_per_chrom(cns_df, samples, cn_columns)
    ane_cols = get_ane_cols_if_none(cns_df, cn_columns)
    autosomes_sum, sex_chrom_sum = calc_ane_per_sample(pre_chr, ane_cols, assembly)
    autosomes_sum = norm_aut_aneuploidy(autosomes_sum, ane_cols, assembly)
    sex_chrom_sum = norm_sex_aneuploidy(samples, sex_chrom_sum, ane_cols, assembly)
    merged_df = autosomes_sum.merge(sex_chrom_sum, left_index=True, right_index=True, suffixes=('_aut', '_sex'))
    res = samples.merge(merged_df, left_index=True, right_index=True)
    return res


def main_cluster(cns_df, dist, assembly=hg19, print_progress=False):    
    dict_start = get_breaks(cns_df)

    if print_progress:
        orig_count = sum(len(values) for values in dict_start.values())
        print(f"Reducing {orig_count} breakpoints, merge distance {dist} ... ")   

    res = created_merged_segs(dict_start, dist, assembly)

    if print_progress:
        new_count = len(res) - len(dict_start) # number of segments is breakpoints + 1 per chrom
        reduction = np.round(new_count / orig_count, 2)
        print(f"Merged breakpoints: {new_count}, reduced by {1 - reduction:.2%}.")

    df = pd.DataFrame(res, columns=["chrom", "start", "end"]).set_index("chrom")
    return df


def regions_select(select, assembly=hg19):    
    if select == "arms":
        breaks = calc_arm_breaks(assembly)
        return breaks_to_segments(breaks)
    elif select == "bands":
        breaks = calc_cytoband_breaks(assembly)
        return breaks_to_segments(breaks)
    elif select =="":
        return genome_to_segments(assembly)
    else:
        return load_regions(select)
    

def regions_remove(remove, assembly=hg19):
    if remove == "gaps":
        breaks = assembly.gaps
        return tuples_to_segments(breaks)
    elif remove == "":
        return None
    else:
        return load_regions(remove)


def get_genome_segments(select, bin_size=0, remove=None, filter_size=0):
    res = select
    if filter_size > 0:
        res = filter_min_size(res, filter_size)
    if remove != None:
        if filter_size > 0:
            remove = filter_min_size(remove, filter_size)
        res = segment_difference(res, remove)
        if filter_size > 0:
            res = filter_min_size(res, filter_size)
    if bin_size > 0:
        res = split_segments(res, bin_size)
    return res

