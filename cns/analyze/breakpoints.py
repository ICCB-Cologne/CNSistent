from cns.analyze.aneuploidy import get_ane_for_col
from cns.process.imputation import merge_neighbours
from cns.utils.assemblies import hg19
import numpy as np
import pandas as pd


def prepare_segments(cns_df, cn_col):
    cns_subset_df = cns_df[["sample_id", "chrom", "start", "end",  cn_col]]
    merged_cns_df = merge_neighbours(cns_subset_df, cn_col, False)
    merged_cns_df["length"] = merged_cns_df["end"] - merged_cns_df["start"]
    return merged_cns_df


# count segments per chromosome and subtract 1
def calc_breaks_per_chr(cns):
    breaks = cns.reset_index().groupby(["sample_id", "chrom"]).size()
    breaks = breaks.reset_index(name="breaks")
    breaks["breaks"] = breaks["breaks"] - 1
    return breaks


def calc_breaks_per_sample(cns_df, samples_df, cn_col, assembly=hg19):
    res = samples_df.copy()
    breaks_per_chr = calc_breaks_per_chr(cns_df)
    chrom_types = {"aut": assembly.aut_names, "sex": assembly.sex_names}

    for suffix, names in chrom_types.items():
        res[f"breaks_{cn_col}_{suffix}"] = (
            breaks_per_chr.query("chrom in @names")
            .groupby("sample_id")["breaks"]
            .sum()
            .reindex(res.index)
            .fillna(0)
            .astype(np.int64)
        )

    res[f"breaks_{cn_col}_all"] = res[f"breaks_{cn_col}_aut"] + res[f"breaks_{cn_col}_sex"]
    return res


def calc_step_per_chr(cns_df, cn_col):
    groups = cns_df.groupby(["sample_id", "chrom"])
    res = []
    for group in groups:
        vals = group[1][cn_col].values
        val = 0 if len(vals) < 2 else np.abs(np.diff(vals)).sum()
        res.append((group[0][0], group[0][1], val, len(vals) - 1))
    return pd.DataFrame(res, columns=["sample_id", "chrom", "step", "count"])


def calc_step_per_sample(cns_df, samples_df, cn_col, assembly=hg19):
    res = samples_df.copy()
    chrom_types = {"aut": assembly.aut_names, "sex": assembly.sex_names, "tot": assembly.chr_names}
    step_per_chr = calc_step_per_chr(cns_df, cn_col)

    for suffix, names in chrom_types.items():
        res[f"step_{cn_col}_{suffix}"] = (
            step_per_chr.query("chrom in @names").groupby("sample_id")
            .apply(lambda x: x["step"].sum() / x["count"].sum() if x["count"].sum() != 0 else 0)
            .reset_index(name="cnstep").set_index("sample_id")
        )
    return res.fillna(0)


def calc_seg_size_per_sample(cns_df, samples_df, cn_col, assembly=hg19):
    res = samples_df.copy()    
    chrom_types = {"aut": assembly.aut_names, "sex": assembly.sex_names, "tot": assembly.chr_names}
    mask = cns_df.apply(lambda x: get_ane_for_col(cn_col, x, samples_df, assembly), axis=1)
    for suffix, names in chrom_types.items():
        subset = cns_df[mask].query("chrom in @names")
        res[f"segsize_{cn_col}_{suffix}"] = subset.groupby("sample_id")["length"].mean()
    return res.fillna(0)