import matplotlib.pyplot as plt
import numpy as np

from cns.display.label import plot_chr_bg, plot_x_ticks
from cns.utils.conversions import column_to_label
from cns.utils.assemblies import hg19


# Will plot lines between consecutive segments
def line_plot(ax, grouped, column, color="red", label=None, alpha=1, line_width=1, sel_chrom=None):
    chroms = grouped["chrom"].unique() if sel_chrom is None else [sel_chrom]
    for chr in chroms:
        df = grouped.query(f"chrom == '{chr}'").copy()
        df["is_consecutive"] = df["start"] - df["end"].shift(1) != 0
        # plot consecutive segments
        for _, group_df in df.groupby(df["is_consecutive"].cumsum()):
            ax.plot(
                group_df["cum_mid" if sel_chrom is None else "mid"],
                group_df[column],
                c=color,
                linewidth=line_width,
                label=label,
                alpha=alpha,
            )
            label = None  # only use label for the first chromosome
    return ax


# TODO: Scale by the length of the segment
def scatter_plot(ax, grouped, column, color="red", label=None, alpha=1, dot_size=1):
    ax.scatter(
        grouped["cum_mid"],
        grouped[column],
        s=dot_size,
        label=label,
        color=color,
        alpha=alpha,
    )


def fig_genome(
    grouped,
    assembly=hg19,
    label=None,
    column="total_cn",
    min_cn=0,
    max_cn=10,
    width=16,
    dpi=100,
    colored=False,
    size=1,
    alpha=1,
    color="red"
):
    fig, ax = plt.subplots(1, figsize=(width, width / 4), dpi=dpi)
    line_width = size * width * 100 / len(grouped)
    plot_chr_bg(ax, assembly, min_cn, max_cn, colored)
    plot_x_ticks(ax, assembly)
    line_plot(ax, grouped, column, color, label, alpha, line_width)
    if label != None:
        ax.legend(loc="upper right")
    y_label = column_to_label(column)
    ax.set_ylabel(y_label, fontdict={"fontsize": width})
    ax.set_xlabel(
        "position on a chromosome", fontdict={"fontsize": width}, labelpad=width / 2
    )
    ax.set_title(f"{y_label} per bin on linear genome")
    return fig, ax


def fig_genome_groups(
    grouped_per_group,
    label_per_group,
    assembly=hg19,
    column="total_cn",
    min_cn=0,
    max_cn=10,
    width=16,
    dpi=100,
    colored=False,
    size=1,
    alpha=1
):
    fig, ax = plt.subplots(1, figsize=(width, width / 4), dpi=dpi)
    group_count = len(grouped_per_group)
    assert group_count == len(label_per_group)
    colors = plt.cm.nipy_spectral(np.linspace(0.05, 0.95, group_count))
    line_width = size * width * 100 / len(grouped_per_group[0])
    plot_chr_bg(ax, assembly, min_cn, max_cn, colored)
    plot_x_ticks(ax, assembly)
    alpha *= (1 / group_count) ** (1 / 3)
    for i, group in enumerate(grouped_per_group):
        color = colors[i]
        label = label_per_group[i]
        line_plot(ax, group, column, color, label, alpha, line_width)
    y_label = column_to_label(column)
    ax.set_ylabel(y_label, fontdict={"fontsize": width})
    ax.set_xlabel(
        "position on a chromosome", fontdict={"fontsize": width}, labelpad=width / 2
    )
    ax.set_title(f"{y_label} per bin on linear genome")
    ax.legend(loc="upper right")
    return fig, ax


def fig_manhattan(
    grouped,
    label=None,
    assembly=hg19,
    column="total_cn",
    min_cn=0,
    max_cn=10,
    width=16,
    dpi=100,
    colored=False,
    size=1,
    alpha=0.3,
    color="red"
):
    fig, ax = plt.subplots(1, figsize=(width, width / 4), dpi=dpi)
    plot_chr_bg(ax, assembly, min_cn, max_cn, colored)
    plot_x_ticks(ax, assembly)
    dot_size = size * width * 100 / len(grouped)
    scatter_plot(ax, grouped, column, color, label, alpha, dot_size)
    if label != None:
        ax.legend(loc="upper right")
    y_label = column_to_label(column)
    ax.set_ylabel(y_label, fontdict={"fontsize": width})
    ax.set_xlabel(
        "position on a chromosome", fontdict={"fontsize": width}, labelpad=width / 2
    )
    ax.set_title(f"{y_label} per bin on linear genome")
    return fig, ax


def fig_manhattan_groups(
    grouped_per_group,
    label_per_group,
    assembly=hg19,
    column="total_cn",
    min_cn=0,
    max_cn=10,
    width=16,
    dpi=100,
    colored=False,
    size=1,
    alpha=0.3,
    color="red"
):
    fig, ax = plt.subplots(1, figsize=(width, width / 4), dpi=dpi)
    group_count = len(grouped_per_group)
    assert group_count == len(label_per_group)
    plot_chr_bg(ax, assembly, min_cn, max_cn, colored)
    plot_x_ticks(ax, assembly)
    colors = plt.cm.nipy_spectral(np.linspace(0.05, 0.95, group_count))
    alpha *= (1 / group_count) ** (1 / 3)
    dot_size = size * width * 100 / len(grouped_per_group[0])
    for i, grouped in enumerate(grouped_per_group):
        color = colors[i]
        label = label_per_group[i]
        scatter_plot(ax, grouped, column, color, label, alpha, dot_size)

    y_label = column_to_label(column)
    ax.set_ylabel(y_label, fontdict={"fontsize": width})
    ax.set_xlabel(
        "position on a chromosome", fontdict={"fontsize": width}, labelpad=width / 2
    )
    ax.set_title(f"{y_label} per bin on linear genome")
    ax.legend(loc="upper right")
    return fig, ax
