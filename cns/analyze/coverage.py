import numpy as np
import pandas as pd
from cns.utils.conversions import calc_len
from cns.utils.selection import only_aut, only_sex
from cns.utils.assemblies import hg19


def get_not_nan(cns_df, cn_columns, het):
    nan_vals = cns_df[cn_columns].isna()
    nan_filter = ~nan_vals.all(axis=1) if het else ~nan_vals.any(axis=1)
    non_nan_df = cns_df.loc[nan_filter]
    return non_nan_df


def get_covered_bases(nan_bases_df, samples_df, het):
    res = samples_df.copy()
    label = "het" if het else "hom"
    aut_df = only_aut(nan_bases_df)
    aut_df_len = calc_len(aut_df)
    sex_df = only_sex(nan_bases_df)
    sex_df_len = calc_len(sex_df)
    # Group the differences by sample_id and compute the sum for each group
    res[f"cover_{label}_aut"] = (
        aut_df_len.groupby(aut_df["sample_id"]).sum().reindex(res.index).fillna(0).astype(np.int64)
    )
    res[f"cover_{label}_sex"] = (
        sex_df_len.groupby(sex_df["sample_id"]).sum().reindex(res.index).fillna(0).astype(np.int64)
    )
    res[f"cover_{label}_all"] = res[f"cover_{label}_aut"] + res[f"cover_{label}_sex"]
    return res


def get_missing_chroms(cns_df, samples_df, segs=None, assembly=hg19):
    res = samples_df.copy()
    # create a serise where the value is sex_xy if expected_chrs == 'xy' lese it is sex_xx
    xy_names = assembly.aut_names + ["chrX", "chrY"]
    xx_names = assembly.aut_names + ["chrX"]
    if segs is not None:
        seg_chrs = segs.keys()
        xy_names = [x for x in xy_names if x in seg_chrs]
        xx_names = [x for x in xx_names if x in seg_chrs]

    expected_chrs = res["sex"].map({"xy": xy_names, "xx": xx_names, "NA": xx_names})
    tot_chrs = cns_df.groupby("sample_id")["chrom"].unique()

    merged = pd.DataFrame([expected_chrs, tot_chrs]).T
    diff = merged.apply(lambda x: np.setdiff1d(x.iloc[0], x.iloc[1]), axis=1)

    res["chrom_count"] = tot_chrs.apply(lambda x: len(x)).reindex(res.index).fillna(0).astype(np.int64)
    res["chrom_missing"] = diff
    return res
