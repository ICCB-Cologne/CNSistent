import numpy as np
from cns.process.normalize import get_chr_sets
from cns.utils.assemblies import hg19
from cns.utils.conversions import calc_lenghts


def get_expected_ploidy(column, chrom, is_xy, assembly=hg19):
    if chrom == assembly.chr_x:
        if is_xy:
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
        if is_xy:
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


def get_ane_for_col(col, row, samples_df, assembly=hg19):
    return get_expected_ploidy(col, row["chrom"], samples_df.loc[row["sample_id"]]["sex"] == "xy", assembly) != row[
        col
    ] and not np.isnan(row[col])


def is_seg_ane(cns_df, samples_df, cn_columns, het, assembly=hg19):
    is_ane = [cns_df.apply(lambda x: get_ane_for_col(col, x, samples_df, assembly), axis=1) for col in cn_columns]
    return np.any(is_ane, axis=0) if het else np.all(is_ane, axis=0)


def get_het_loh_for_col(col, row, samples_df, assembly=hg19):
    return get_expected_ploidy(col, row["chrom"], samples_df.loc[row["sample_id"]]["sex"] == "xy", assembly) > row[col]


def get_hom_loh_for_col(col, row, samples_df, assembly=hg19):
    return get_expected_ploidy(col, row["chrom"], samples_df.loc[row["sample_id"]]["sex"] == "xy", assembly) != 0 and row[col] == 0


def is_seg_loh(cns_df, samples_df, cn_columns, het, assembly=hg19):
    if len(cn_columns) == 2:
        is_ane = [cns_df.apply(lambda x: get_het_loh_for_col(col, x, samples_df, assembly), axis=1) for col in cn_columns]
        return np.any(is_ane, axis=0) if het else np.all(is_ane, axis=0)
    else:
        fun = get_het_loh_for_col if het else get_hom_loh_for_col
        is_loh = cns_df.apply(lambda x: fun(cn_columns[0], x, samples_df, assembly), axis=1)
        return is_loh
    

def _get_check_fun(feature):
    if feature == "ane":
        return is_seg_ane
    elif feature == "loh":
        return is_seg_loh
    else:
        raise ValueError("feature must be one of ['ane', 'loh']")


def _calc_bases_per_chr_group(res, masked_cns_df, label, groups, assembly=hg19):
    for suffix, names in groups.items():
        subset = masked_cns_df.query("chrom in @names")
        length = calc_lenghts(subset)
        res[f"{label}_{suffix}"] = length.groupby(subset["sample_id"]).sum()
        res[f"{label}_{suffix}"] = res[f"{label}_{suffix}"].fillna(0).astype(np.int64)
    return res


def _calc_bases_per_column(res, cns_df, cn_columns, het, feature, assembly=hg19):
    label = feature + "_" + ("het" if het else "hom")
    function = _get_check_fun(feature)
    mask = function(cns_df, res, cn_columns, het, assembly)
    chr_sets = get_chr_sets(cns_df, assembly)
    return _calc_bases_per_chr_group(res, cns_df[mask], label, chr_sets, assembly)


def calc_ane_bases(cns_df, samples_df, cn_columns, assembly=hg19):
    res = samples_df.copy()
    if len(cn_columns) == 2:
        res = _calc_bases_per_column(res, cns_df, cn_columns, False, "ane", assembly)
    res = _calc_bases_per_column(res, cns_df, cn_columns, True, "ane", assembly)
    return res


def calc_loh_bases(cns_df, samples_df, cn_columns, assembly=hg19):
    res = samples_df.copy()
    res = _calc_bases_per_column(res, cns_df, cn_columns, False, "loh", assembly)
    res = _calc_bases_per_column(res, cns_df, cn_columns, True, "loh", assembly)
    return res


def calc_imb_bases(cns_df, samples_df, cn_columns, col_index=0, assembly=hg19):
    res = samples_df.copy()
    if len(cn_columns) != 2:
        raise ValueError("There must be two CN columns to calculate imbalance score")
    cn_col1 = cn_columns[col_index]
    cn_col2 = cn_columns[1 - col_index]
    mask = cns_df[cn_col1] > cns_df[cn_col2]
    label = "imb_" + cn_col1
    chr_sets = get_chr_sets(cns_df, assembly)
    res = _calc_bases_per_chr_group(res, cns_df[mask], label, chr_sets, assembly)
    return res
