import numpy as np
from cns.utils.assemblies import hg19
from cns.utils.selection import only_aut, only_sex


def get_expected_ploidy(column, chrom, is_xy, assembly=hg19):
    if chrom == assembly.chr_x:
        if is_xy:
            if column == "minor_cn" or column == "hapY_cn":
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
            if column == "minor_cn" or column == "hapY_cn":
                return 0
            else:
                return 1
        else:
            return 0
    else:
        if column == "total_cn":
            return 2
        else:
            return 1
        

def get_ane_for_cols(cns_df, samples_df, cn_columns, assembly=hg19):
    is_ane = [cns_df.apply(lambda x: 
                           get_expected_ploidy( col, x["chrom"], samples_df.loc[x["sample_id"]]["sex"] == "xy", assembly) != x[col], 
                           axis=1) for col in cn_columns]
    return is_ane


def get_ane_for_samples(cns_df, samples_df, cn_columns, any=True, assembly=hg19):
    is_ane = get_ane_for_cols(cns_df, samples_df, cn_columns, assembly)
    cns_df["ane"] = np.any(is_ane, axis=0) if any else np.all(is_ane, axis=0)
    res = cns_df[cns_df["ane"]].groupby("sample_id")["length"].sum()
    return res


# Note: missing values are NOT considered to be aneuploid, to consider missing segments, first impute
def get_ane_bases(cns_df, samples_df, cn_columns, any=True, assembly=hg19):
    res = samples_df.copy()
    cns_df["length"] = cns_df["end"] - cns_df["start"]

    aut_df = only_aut(cns_df).copy()
    sex_df = only_sex(cns_df).copy()

    # Group the differences by sample_id and compute the sum for each group
    res["ane_aut"] = get_ane_for_samples(aut_df, samples_df, cn_columns, any)
    res["ane_aut"] = res["ane_aut"].fillna(0).astype(np.int64)
    res["ane_sex"] = get_ane_for_samples(sex_df, samples_df, cn_columns, any)
    res["ane_sex"] = res["ane_sex"].fillna(0).astype(np.int64)
    res["ane_tot"] = res["ane_aut"] + res["ane_sex"]
    return res.fillna(0)
