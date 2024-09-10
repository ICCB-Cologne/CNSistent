import numpy as np
from cns.utils.assemblies import hg19
from cns.utils.selection import only_aut, only_sex


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
    return get_expected_ploidy(col, row["chrom"], samples_df.loc[row["sample_id"]]["sex"] == "xy", assembly) != row[col]


def get_ane_for_cols(cns_df, samples_df, cn_columns, assembly=hg19):
    is_ane = [cns_df.apply(lambda x: get_ane_for_col(col, x, samples_df, assembly), axis=1) for col in cn_columns]
    return is_ane


def get_ane_for_samples(cns_df, samples_df, cn_columns, het, assembly=hg19):
    is_ane = get_ane_for_cols(cns_df, samples_df, cn_columns, assembly)
    cns_df["ane"] = np.any(is_ane, axis=0) if het else np.all(is_ane, axis=0)
    res = cns_df[cns_df["ane"]].groupby("sample_id")["length"].sum()
    return res


def get_ane_for_chroms(res, aut_df, sex_df, samples_df, cn_columns, het, assembly=hg19):
    label = "het" if het else "hom"
    # Group the differences by sample_id and compute the sum for each group
    res[f"ane_{label}_aut"] = get_ane_for_samples(aut_df, samples_df, cn_columns, het, assembly).reindex(res.index).fillna(0).astype(np.int64)
    res[f"ane_{label}_sex"] = get_ane_for_samples(sex_df, samples_df, cn_columns, het, assembly).reindex(res.index).fillna(0).astype(np.int64)
    res[f"ane_{label}_tot"] = res[f"ane_{label}_aut"] + res[f"ane_{label}_sex"]
    return res

# Note: missing values are NOT considered to be aneuploid, to consider missing segments, first impute
def get_ane_bases(cns_df, samples_df, cn_columns, assembly=hg19):
    res = samples_df.copy()
    cns_df["length"] = cns_df["end"] - cns_df["start"]

    aut_df = only_aut(cns_df).copy()
    sex_df = only_sex(cns_df).copy()

    if len(cn_columns) == 2:
        res = get_ane_for_chroms(res, aut_df, sex_df, samples_df, cn_columns, True, assembly)
    res = get_ane_for_chroms(res, aut_df, sex_df, samples_df, cn_columns, False, assembly)    
    return res
