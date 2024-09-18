
from cns.process.segments import calc_group_sizes, hg19
import numpy as np
from cns.utils.assemblies import hg19

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


def get_chr_sets(cns_df, assembly=hg19):
    chroms = cns_df["chrom"].unique().tolist()
    aut_selection = [chrom for chrom in chroms if chrom in assembly.aut_names]
    if len(aut_selection) == 0:
        print(cns_df)
        raise ValueError("No autosomes found in the input segments.")    
    res_dict = { "aut": aut_selection}
    sex_selection = [chrom for chrom in chroms if chrom in assembly.sex_names]
    if len(sex_selection) != 0:
        res_dict["sex"] = sex_selection
        res_dict["all"] = aut_selection + sex_selection
    return res_dict
    
    