import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from cns.display.label import plot_cytobands, plot_x_ticks
from cns.utils.selection import only_aut, only_sex
from cns.utils.assemblies import hg19


def plot_nans(nans_df, dataset_name):
    fig, ax = plt.subplots(1, figsize=(16, 4), dpi = 200)
    nans_df = nans_df.copy()
    plot_x_ticks(ax)
    ax.hist(nans_df["cum_mid"], weights=nans_df["len"], bins=310)
    # get current max y
    max_y = ax.get_ylim()[1]
    plot_cytobands(ax, 0, max_y)
    ax.hist(nans_df["cum_mid"], weights=nans_df["len"], bins=310, color="blue")
    ax.set_title(f"Distribution of NaNs by chromosome in {dataset_name} dataset")
    ax.set_xlabel("Chromosome and position")
    ax.set_ylabel("Length of NaNs")


def plot_coverage(samples_df, dataset, lower_lim=0.9):
    # create 2 rows of 2 plots each
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    included = sum(samples_df["tot_frac"] >= lower_lim)
    total = len(samples_df)

    top_range = np.linspace(lower_lim, 1.0, 51)

    # plots for "tot_frac"
    axes[0].hist(samples_df["tot_frac"], bins=top_range)
    axes[0].set_title(
        "Histogram of total coverage fraction ({} of {})".format(included, total)
    )
    axes[0].set_ylabel("Number of samples")
    axes[0].set_xlabel("Total coverage fraction")

    # plots for "aut_frac"
    axes[1].hist(samples_df["aut_frac"], bins=top_range)
    axes[1].set_title(
        "Histogram of autosome coverage fraction ({} of {})".format(included, total)
    )
    axes[1].set_xlabel("Autosome coverage fraction")

    plt.suptitle(f"Coverage {dataset} dataset ({len(samples_df)})", fontsize=16)


def plot_segments_per_chr(info_df, dataset):
    # create two plots next to each other, get axes
    plt.figure(figsize=(12, 5))
    plt.subplot(1, 2, 1)
    plt.hist(info_df["extra_segs"], bins=50)
    plt.xlabel("Extra segments")
    plt.ylabel("Count")
    plt.title("Extra segments (segments-chromosomes)")
    plt.subplot(1, 2, 2)
    plt.hist(info_df["anueploidies"], bins=50)
    plt.xlabel("Anueploidies")
    plt.ylabel("Count")
    plt.title("Anueploidies (on whole chromosomes)")

    # set title to "CNS of the PCAWG dataset"
    plt.suptitle(f"CINs of the {dataset} dataset ({len(info_df)})", fontsize=16)


def get_base_frac(samples, assembly=hg19):
    res = samples.copy()
    x_length = assembly.chr_lens["chrX"]
    y_length = assembly.chr_lens["chrY"]
    res["cover_frac_aut"] = res["cover_bases_aut"] / assembly.aut_len
    res["cover_frac_sex"] = np.where(
        res["sex"] == "xy",
        res["cover_bases_sex"] / (x_length + y_length),
        res["cover_bases_sex"] / (x_length),
    )
    res["cover_frac_tot"] = np.where(
        res["sex"] == "xy",
        res["cover_bases_tot"] / (assembly.aut_len + x_length + y_length),
        res["cover_bases_tot"] / (assembly.aut_len + x_length),
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
    res["cover_bases_aut"] = aut_lens.groupby(aut_df["sample_id"]).sum()
    res["cover_bases_aut"] = res["cover_bases_aut"].fillna(0).astype(np.int64)
    res["cover_bases_sex"] = sex_lens.groupby(sex_df["sample_id"]).sum()
    res["cover_bases_sex"] = res["cover_bases_sex"].fillna(0).astype(np.int64)
    res["cover_bases_tot"] = res["cover_bases_aut"] + res["cover_bases_sex"]
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
