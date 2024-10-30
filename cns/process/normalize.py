import numpy as np
from cns.utils import hg19


def normalize_feature(samples, feature, norm_sizes):
    res_df = samples.copy()
    if norm_sizes["aut"] == 0:
        raise ValueError("Total autosome size is 0 - check the input segments.")
    res_df[f"{feature}_aut"] = res_df[f"{feature}_aut"] / norm_sizes["aut"]
    if f"{feature}_sex" in res_df.columns:
        res_df[f"{feature}_sex"] = np.where(
            res_df["sex"] == "xy",
            res_df[f"{feature}_sex"] / norm_sizes["sexXY"] if norm_sizes["sexXY"] != 0 else 0,
            res_df[f"{feature}_sex"] / norm_sizes["sexXX"] if norm_sizes["sexXX"] != 0 else 0,
        )
        res_df[f"{feature}_all"] = np.where(
            res_df["sex"] == "xy",
            res_df[f"{feature}_all"] / norm_sizes["allXY"],
            res_df[f"{feature}_all"] / norm_sizes["allXX"],
        )
    return res_df


def calc_chr_sizes(chr_segs_df):
    return {chr: sum([seg[1] - seg[0] for seg in chr_segs]) for chr, chr_segs in chr_segs_df.items()}


def calc_group_sizes(segs_df, assembly=hg19):
    chr_sizes = calc_chr_sizes(segs_df)
    groups = {"sexXX": [assembly.chr_x], "aut": assembly.aut_names, "sexXY": [assembly.chr_y, assembly.chr_x]}
    return {
        key: sum([chr_sizes[chrom] if chrom in chr_sizes else 0 for chrom in group]) for key, group in groups.items()
    }


def get_norm_sizes(segs, assembly=hg19):
    sizes = (
        {
            "aut": assembly.aut_len,
            "sexXX": assembly.chr_lens[assembly.chr_x],
            "sexXY": assembly.chr_lens[assembly.chr_y] + assembly.chr_lens[assembly.chr_x],
        }
        if segs is None
        else calc_group_sizes(segs, assembly)
    )
    sizes["allXX"] = sizes["aut"] + sizes["sexXX"]
    sizes["allXY"] = sizes["aut"] + sizes["sexXY"]
    return sizes


    
    