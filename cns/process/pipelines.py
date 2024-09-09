import pandas as pd
import numpy as np

from cns.analyze.aneuploidy import get_ane_bases
from cns.analyze.coverage import normalize_feature, get_covered_bases, get_missing_chroms, get_non_nan
from cns.analyze.signatures import add_breaks_per_sample
from cns.process.binning import bin_by_segments
from cns.process.breakpoints import calc_arm_breaks, calc_cytoband_breaks, get_breaks
from cns.process.cluster import created_merged_segs
from cns.process.imputation import add_missing, add_tails, cns_impute, fill_gaps, fill_nans_with_zeros, merge_neighbours
from cns.process.segments import filter_min_size, segment_difference, split_segments
from cns.utils.conversions import genome_to_segments, breaks_to_segments, tuples_to_segments
from cns.utils.files import load_regions, samples_df_from_cns_df, find_cn_cols_if_none, rename_cn_cols
from cns.utils.logging import log_info
from cns.utils.assemblies import hg19


def main_fill(cns_df, samples_df=None, cn_columns=None, assembly=hg19, add_missing_chromosomes=True, print_info=False):
    if samples_df is None:
        log_info(print_info, print("No samples provided, creating samples from CNS data."))
        samples_df = samples_df_from_cns_df(cns_df)
    cn_columns = find_cn_cols_if_none(cns_df, cn_columns)
    cns_tailed_df = add_tails(cns_df, assembly.chr_lens, print_info=print_info)
    cns_filled_df = fill_gaps(cns_tailed_df, print_info=print_info)
    if add_missing_chromosomes:
        cns_filled_df = add_missing(cns_filled_df, samples_df, assembly.chr_lens, print_info=print_info)
    res = merge_neighbours(cns_filled_df, cn_columns, print_info=print_info)
    return res


def main_impute(cns_df, samples_df=None, cn_columns=None, method='extend', print_info=False):
    cn_columns = find_cn_cols_if_none(cns_df, cn_columns)
    if method == 'diploid' and samples_df is None:
        log_info(print_info, print("Diploid imputation requires samples, but none provided, creating samples from CNS data."))
        samples_df = samples_df_from_cns_df(cns_df)
    imputed_df = cns_impute(cns_df, samples_df, method, cn_columns=cn_columns, print_info=print_info)
    filled_df = fill_nans_with_zeros(imputed_df, cn_columns=cn_columns, print_info=print_info)
    res = merge_neighbours(filled_df, cn_columns=cn_columns, print_info=print_info)
    assert len(res[res.isnull().any(axis=1)]) == 0, "NaNs still present in final_df."
    return res


def main_bin(cns_df, segs, fun_type='mean', print_info=False):
    return bin_by_segments(cns_df, segs, fun_type, print_info)


# any: if True, based is considered as covered if any CN column has values assigned
def main_coverage(cns_df, samples_df, cn_columns=None, any=True, assembly=hg19, print_info=False):
    cn_columns = find_cn_cols_if_none(cns_df, cn_columns)
    # Select the rows where copy-numbers are not Not a Number (NaN == NaN) is false
    non_nan_df = get_non_nan(cns_df, cn_columns, any)
    samples_df = get_missing_chroms(non_nan_df, samples_df, assembly)
    samples_df = get_covered_bases(non_nan_df, samples_df)
    samples_df = normalize_feature(samples_df, "cover", assembly)
    return samples_df


def main_signatures(cns_df, samples_df, assembly=hg19, print_info=False):
    res = add_breaks_per_sample(cns_df, samples_df, assembly)
    return res


def main_ploidy(cns_df, samples_df, cn_columns=None, assembly=hg19, print_info=False):
    cns_df, cn_columns = rename_cn_cols(cns_df, cn_columns, assembly)
    samples_df = get_ane_bases(cns_df, samples_df, cn_columns, assembly)
    samples_df = normalize_feature(samples_df, "ane_hom", assembly)
    if len(cn_columns) == 2:
        samples_df = normalize_feature(samples_df, "ane_het", assembly)
    return samples_df


def main_cluster(cns_df, dist, assembly=hg19, print_info=False):    
    dict_start = get_breaks(cns_df, keep_ends=False, assembly=assembly)

    if print_info:
        orig_count = sum(len(values) for values in dict_start.values())
        print(f"Reducing {orig_count} breakpoints, merge distance {dist} ... ")   

    res = created_merged_segs(dict_start, dist, assembly, extend=True)

    if print_info:
        new_count = len(res) - len(dict_start) # number of segments is breakpoints + 1 per chrom
        reduction = np.round(new_count / orig_count, 2)
        print(f"Merged breakpoints: {new_count}, reduced by {1 - reduction:.2%}.")

    df = pd.DataFrame(res, columns=["chrom", "start", "end"])
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

