import numpy as np
import pandas as pd
from cns.utils.selection import only_aut, only_sex
from cns.utils.assemblies import hg19

def normalize_feature(samples, feature, assembly=hg19):
    res = samples.copy()
    x_length = assembly.chr_lens[assembly.chr_x]
    y_length = assembly.chr_lens[assembly.chr_y]
    res[f"{feature}_aut"] = res[f"{feature}_aut"] / assembly.aut_len
    res[f"{feature}_sex"] = np.where(
        res["sex"] == "xy",
        res[f"{feature}_sex"] / (x_length + y_length),
        res[f"{feature}_sex"] / (x_length),
    )
    res[f"{feature}_tot"] = np.where(
        res["sex"] == "xy",
        res[f"{feature}_tot"] / (assembly.aut_len + x_length + y_length),
        res[f"{feature}_tot"] / (assembly.aut_len + x_length),
    )
    return res


def get_covered_bases(cns_df, samples):
    res = samples.copy()
    aut_df = only_aut(cns_df)
    sex_df = only_sex(cns_df)
    # Compute the differences between end and start
    aut_lens = aut_df["end"] - aut_df["start"]
    sex_lens = sex_df["end"] - sex_df["start"]

    # Group the differences by sample_id and compute the sum for each group
    res["cover_aut"] = aut_lens.groupby(aut_df["sample_id"]).sum()
    res["cover_aut"] = res["cover_aut"].fillna(0).astype(np.int64)
    res["cover_sex"] = sex_lens.groupby(sex_df["sample_id"]).sum()
    res["cover_sex"] = res["cover_sex"].fillna(0).astype(np.int64)
    res["cover_tot"] = res["cover_aut"] + res["cover_sex"]
    return res


def get_missing_chroms(cns, samples, assembly=hg19):
    res = samples.copy()
    # create a serise where the value is sex_xy if exhpected_chrs == 'xy' lese it is sex_xx
    xy_names = assembly.aut_names + ["chrX", "chrY"]
    xx_names = assembly.aut_names + ["chrX"]
    expected_chrs = res["sex"].map({"xy": xy_names, "xx": xx_names, "NA": xx_names})
    tot_chrs = cns.groupby("sample_id")["chrom"].unique()

    merged = pd.DataFrame([expected_chrs, tot_chrs]).T
    diff = merged.apply(lambda x: np.setdiff1d(x.iloc[0], x.iloc[1]), axis=1)

    res["chrom_count"] = tot_chrs.apply(lambda x: len(x))
    res["chrom_missing"] = diff
    return res
