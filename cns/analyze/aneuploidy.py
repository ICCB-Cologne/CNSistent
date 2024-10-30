import numpy as np
from cns.process.normalize import get_chr_sets
from cns.utils.assemblies import hg19
from cns.utils.conversions import calc_lenghts


def get_expected_ploidy(column, chrom, sex, assembly=hg19):
    if chrom == assembly.chr_x:
        if sex == "xy":
            if column == "minor_cn" or column == "hap2":
                return 0
            else:
                return 1
        else:
            if column == "total_cn":
                return 2
            else:
                return 1
    elif chrom == assembly.chr_y:
        if sex == "xy":
            if column == "major_cn" or column == "hap1":
                return 1
            else:
                return 0
        else:
            return 0
    else:
        if column == "total_cn":
            return 2
        else:
            return 1


def get_ane_for_col(col, cns_row, sample_row, assembly=hg19):
    return get_expected_ploidy(col, cns_row["chrom"], sample_row["sex"], assembly) != cns_row[col] and not np.isnan(cns_row[col])


def get_het_loh_for_col(col, cns_row, sample_row, assembly=hg19):
    return get_expected_ploidy(col, cns_row["chrom"], sample_row["sex"], assembly) > cns_row[col]


def get_hom_loh_for_col(col, cns_row, sample_row, assembly=hg19):
    return get_expected_ploidy(col, cns_row["chrom"], sample_row["sex"], assembly) != 0 and cns_row[col] == 0
    

def is_total_seg_loh(row, col_name, sex, assembly):
    val = row[col_name]
    chrom = row["chrom"]
    if sex == "xy":
        if chrom == assembly.chr_x or chrom == assembly.chr_y:
            return val < 1
        else:
            return val < 2
    else:
        if chrom == assembly.chr_y:
            return val < 0
        else:    
            return val < 2
        

def is_total_seg_ane(row, col_names, sex, assembly):
    val = row[col_names[0]]
    chrom = row["chrom"]
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
        

def are_hap_seg_loh(row, col_names, sex, is_het, assembly):
    vals = sorted(row[col_names], reverse=True)
    if sex == "xy":
        if row["chrom"] == assembly.chr_x or row["chrom"] == assembly.chr_y:
            expected = [0, 1]
        else:
            expected = [1, 1]	
    else:
        if row["chrom"] == assembly.chr_y:
            expected = [0, 0]
        else:
            expected = [1, 1]
    loh = [vals[0] < expected[0], vals[1] < expected[1]]
    return any(loh) if is_het else all(loh)
        

def are_hap_seg_ane(row, col_names, sex, is_het, assembly):
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
    ane = [vals[0] != expected[0], vals[1] != expected[1]]
    return any(ane) if is_het else all(ane)


def get_ane_per_seg(cns_df, samples_df, cn_columns, is_het, assembly=hg19):
    if len(cn_columns) == 1:
        return cns_df.apply(lambda row: is_total_seg_ane(row, cn_columns, samples_df.loc[row["sample_id"]]["sex"], assembly), axis=1)       
    else:
        return cns_df.apply(lambda row: are_hap_seg_ane(row, cn_columns, samples_df.loc[row["sample_id"]]["sex"], is_het, assembly), axis=1)


def get_loh_per_seg(cns_df, samples_df, cn_columns, is_het, assembly=hg19):
    if len(cn_columns) == 1:
        return cns_df.apply(lambda row: is_total_seg_loh(row, cn_columns, samples_df.loc[row["sample_id"]]["sex"], assembly), axis=1)       
    else:
        return cns_df.apply(lambda row: are_hap_seg_loh(row, cn_columns, samples_df.loc[row["sample_id"]]["sex"], is_het, assembly), axis=1)


def _get_check_fun(feature):
    if feature == "ane":
        return get_ane_per_seg
    elif feature == "loh":
        return get_loh_per_seg
    else:
        raise ValueError("feature must be one of ['ane', 'loh']")


def _calc_bases_per_chr_group(res, masked_cns_df, label, groups):
    for suffix, names in groups.items():
        subset = masked_cns_df.query("chrom in @names")
        length = calc_lenghts(subset)
        res[f"{label}_{suffix}"] = length.groupby(subset["sample_id"]).sum()
        res[f"{label}_{suffix}"] = res[f"{label}_{suffix}"].fillna(0).astype(np.int64)
    return res


def _count_bases_with_feature(res_df, cns_df, cn_columns, allele_spec, feature, assembly):
    label = feature + "_" + allele_spec
    function = _get_check_fun(feature)
    mask = function(cns_df, res_df, cn_columns, allele_spec == "het", assembly)
    chr_sets = get_chr_sets(cns_df, assembly)
    return _calc_bases_per_chr_group(res_df, cns_df[mask], label, chr_sets)


def calc_loh_bases(samples_df, cns_df, cn_columns, allele_spec, assembly=hg19):
    res_df = samples_df.copy()
    return _count_bases_with_feature(res_df, cns_df, cn_columns, allele_spec, "loh", assembly)


def calc_ane_bases(samples_df, cns_df, cn_columns, allele_spec, assembly=hg19):
    res_df = samples_df.copy()
    return _count_bases_with_feature(res_df, cns_df, cn_columns, allele_spec, "ane", assembly)


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
