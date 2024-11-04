import numpy as np
import pandas as pd
from cns.utils.assemblies import hg19
from cns.utils.conversions import calc_lengths
from cns.utils.selection import get_chr_sets


def _check_total(row, col_names, sex, feature, allele_spec, assembly):
    val = row[col_names[0]]
    chrom = row["chrom"]
    if feature == "ane":
        if sex == "xy":
            if chrom == assembly.chr_x or chrom == assembly.chr_y:
                return val != 1
            else:
                return val != 2
        else:
            if chrom == assembly.chr_y:
                return val != 0
            else:    
                return val != 2
    else:
        if sex == "xy":
            if chrom == assembly.chr_x or chrom == assembly.chr_y:
                return val == 0
            else:
                return val < 2 if allele_spec == "het" else val == 0
        else:
            if chrom == assembly.chr_y:
                return False
            else:    
                return val < 2 if allele_spec == "het" else val == 0
            

def _check_alleles(row, col_names, sex, feature, allele_spec, assembly):
    vals = sorted(row[col_names], reverse=True)
    if sex == "xy":
        if row["chrom"] == assembly.chr_x or row["chrom"] == assembly.chr_y:
            expected = [1, 0]
        else:
            expected = [1, 1]
    else:
        if row["chrom"] == assembly.chr_y:
            expected = [0, 0]
        else:
            expected = [1, 1]
    
    if feature == "ane":
        checks = [vals[0] != expected[0], vals[1] != expected[1]]
    else:  # check_type == "loh"
        checks = [vals[0] < expected[0], vals[1] < expected[1]]
    
    return any(checks) if allele_spec == "het" else all(checks)


def _get_feature_per_seg(cns_df, samples_df, cn_columns, feature, allele_spec, assembly=hg19):
    if len(cn_columns) == 1:
        return cns_df.apply(lambda row: _check_total(row, cn_columns, samples_df.loc[row["sample_id"]]["sex"], feature, allele_spec, assembly), axis=1)       
    else:
        return cns_df.apply(lambda row: _check_alleles(row, cn_columns, samples_df.loc[row["sample_id"]]["sex"], feature, allele_spec, assembly), axis=1)


def _calc_bases_per_chr_group(res, masked_cns_df, label, groups):
    for suffix, names in groups.items():
        subset = masked_cns_df.query("chrom in @names")
        length = calc_lengths(subset)
        res[f"{label}_{suffix}"] = length.groupby(subset["sample_id"]).sum()
        res[f"{label}_{suffix}"] = res[f"{label}_{suffix}"].fillna(0).astype(np.int64)
    return res


def _count_bases_with_feature(res_df, cns_df, cn_columns, feature, allele_spec, assembly):
    label = feature + "_" + allele_spec
    mask = _get_feature_per_seg(cns_df, res_df, cn_columns, feature, allele_spec, assembly)
    chr_sets = get_chr_sets(cns_df, assembly)
    return _calc_bases_per_chr_group(res_df, cns_df[mask], label, chr_sets)


def calc_loh_bases(samples_df, cns_df, cn_columns, allele_spec, assembly=hg19):
    res_df = samples_df.copy()
    return _count_bases_with_feature(res_df, cns_df, cn_columns, "loh", allele_spec, assembly)


def calc_ane_bases(samples_df, cns_df, cn_columns, allele_spec, assembly=hg19):
    res_df = samples_df.copy()
    return _count_bases_with_feature(res_df, cns_df, cn_columns, "ane", allele_spec, assembly)


def calc_imb_bases(cns_df, samples_df, cn_columns, col_index=0, assembly=hg19):
    res = samples_df.copy()
    if len(cn_columns) != 2:
        raise ValueError("There must be two CN columns to calculate imbalance score")
    cn_col1 = cn_columns[col_index]
    cn_col2 = cn_columns[1 - col_index]
    mask = cns_df[cn_col1] > cns_df[cn_col2]
    label = "imb_" + cn_col1
    chr_sets = get_chr_sets(cns_df, assembly)
    res = _calc_bases_per_chr_group(res, cns_df[mask], label, chr_sets)
    return res


def calc_ploidy_per_column(cns_df, cn_column):
    grouped  = cns_df.groupby('sample_id')
    # get the ploidy for each sample
    res = {}
    for sample_id, group_df in grouped:
        legnths = group_df["end"] - group_df["start"]
        ploidy = (group_df[cn_column] * legnths).sum() / legnths.sum()
        res[sample_id] = ploidy
    return pd.Series(dict(res))
