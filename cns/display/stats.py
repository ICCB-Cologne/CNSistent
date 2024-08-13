import numpy as np
from cns.display.label import plot_cytobands, plot_x_ticks
import matplotlib.pyplot as plt


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


def plot_bends(cutoffs, counts, knee_point, elbow_point):
    # set plot aspect ratio 1:1 
    plt.figure(figsize=(6, 6))
    plt.plot(cutoffs, counts)
    plt.axvline(x=cutoffs[knee_point], color='r', linestyle='--')
    plt.axvline(x=cutoffs[elbow_point], color='g', linestyle='--')
    plt.xlabel("Cutoff")
    plt.ylabel("Count")
    plt.title("Knee and Elbow points")
    plt.legend(["Counts", f"Knee: {cutoffs[knee_point]}", f"Elbow: {cutoffs[elbow_point]}"])
    plt.show()