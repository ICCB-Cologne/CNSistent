import pandas as pd
import numpy as np

from cns.analyze.aneuploidy import calc_ane_per_chrom, calc_ane_per_sample, norm_aut_aneuploidy, norm_sex_aneuploidy
from cns.analyze.coverage import get_base_frac, get_covered_bases, get_missing_chroms
from cns.analyze.signatures import add_breaks_per_sample
from cns.process.binning import add_derived
from cns.process.breakpoints import calc_arm_breaks, calc_cytoband_breaks
from cns.process.cluster import created_merged_segs, get_breaks
from cns.process.imputation import add_missing, add_tails, create_imputed_entries, fill_gaps, fill_nans_with_zeros, merge_neighbours
from cns.process.segments import breaks_to_segments, genome_to_segments, tuples_to_segments
from cns.utils.assemblies import get_assembly
from cns.utils.files import load_regions, samples_df_from_cna_df


def main_fill(cna_df, samples_df=None, assembly='hg38', cn_columns=('major_cn', 'minor_cn'), add_missing_chromosomes=True, print_info=False):
    if isinstance(assembly, str):
        assembly = get_assembly(assembly)
    if samples_df is None:
        samples_df = samples_df_from_cna_df(cna_df)

    # Fill missing values with NaNs
    cna_df = cna_df.reset_index(drop=True)
    cna_tailed_df = add_tails(cna_df,assembly.chr_lens, cn_columns, print_info=print_info)
    cna_filled_df = fill_gaps(cna_tailed_df, print_info=print_info)
    if add_missing_chromosomes:
        cna_filled_df = add_missing(cna_filled_df, samples_df, assembly.chr_lens, print_info=print_info)
    res = merge_neighbours(cna_filled_df, cn_columns=cn_columns, print_info=print_info)
    return res


def main_impute(cna_df, cn_columns=('major_cn', 'minor_cn'), print_info=False):
    # Impute
    cna_df = cna_df.reset_index(drop=True)
    imputed_df = create_imputed_entries(cna_df, cn_columns, print_info=print_info)
    filled_df = fill_nans_with_zeros(imputed_df, cn_columns=cn_columns, print_info=print_info)
    res = merge_neighbours(filled_df, cn_columns=cn_columns, print_info=print_info)
    assert len(res[res.isnull().any(axis=1)]) == 0, "NaNs still present in final_df."
    return res


# TODO: should have any and all option (any is nan, or all are nan)
def main_coverage(cna_df, samples_df, assembly, print_info=False):
    # Select the rows where copy-numbers are not Not a Number (NaN == NaN) is false
    cna_vals = cna_df.loc[~cna_df.isna().any(axis=1)].copy()
    coverage = get_missing_chroms(cna_vals, samples_df, assembly)
    coverage = get_covered_bases(cna_vals, coverage)
    coverage = get_base_frac(coverage, assembly)
    return coverage


def main_ploidy(cns, samples, assembly, print_progress=False):
    samples = add_breaks_per_sample(cns, samples, assembly)
    cns = add_derived(cns, assembly)
    pre_chr = calc_ane_per_chrom(cns, samples)
    autosomes_sum, sex_chrom_sum = calc_ane_per_sample(pre_chr, assembly)
    autosomes_sum = norm_aut_aneuploidy(autosomes_sum, assembly)
    sex_chrom_sum = norm_sex_aneuploidy(samples, sex_chrom_sum, assembly)
    merged_df = autosomes_sum.merge(sex_chrom_sum, left_index=True, right_index=True, suffixes=('_aut', '_sex'))
    res = samples.merge(merged_df, left_index=True, right_index=True)
    return res


def main_cluster(cns, dist, assembly, print_progress=False):
    dict_start = get_breaks(cns)

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


def regions_select(select, assembly):    
    if select == "arms":
        breaks = calc_arm_breaks(assembly)
        return breaks_to_segments(breaks)
    elif select == "bands":
        breaks = calc_cytoband_breaks(assembly)
        return breaks_to_segments(breaks)
    elif select =="":
        return genome_to_segments(assembly.chr_lens)
    else:
        return load_regions(select)
    

def regions_remove(remove, assembly):
    if remove == "gaps":
        breaks = assembly.gaps
        return tuples_to_segments(breaks)
    elif remove == "":
        return None
    else:
        return load_regions(remove)

