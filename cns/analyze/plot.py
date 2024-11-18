import numpy as np
import pandas as pd
from matplotlib import pyplot as plt, patches as mpatches
from collections.abc import Sequence
from itertools import accumulate

from cns.utils.cytobands import cytoband_color
from cns.utils.gaps import gap_color
from cns.utils import get_cn_cols, hg19


def _get_CN_color(cn, min_cn, max_cn):
    val_range = max_cn - min_cn
    if cn == 0:
        return (1, 0, 0)  # Red for zero
    else:
        ratio = 1 - min(1, (cn - min_cn) / val_range)
        return (ratio, ratio, 1)  # Blue gradient based on cn


def _get_CN_color_vector(min_cn, max_cn):
    # Vectorize with otypes for float output
    vectorized_func = np.vectorize(lambda cn: _get_CN_color(cn, min_cn, max_cn), otypes=[float, float, float])
    return lambda cn_array: np.stack(vectorized_func(cn_array), axis=-1)


def _get_start(chrom, assembly):
    return assembly.chr_starts[chrom]


def _get_start_vector(assembly):
    return np.vectorize(lambda chrom: _get_start(chrom, assembly), otypes=[np.uint32])


def plot_lines(ax, cns_df, cn_column, color="green", label=None, alpha=1, size=1, assembly=hg19):
    f_start_pos = _get_start_vector(assembly)
    is_consecutive = cns_df["start"] - cns_df["end"].shift(1) != 0
    # plot consecutive segments
    for _, group_df in cns_df.groupby(is_consecutive.cumsum()):
        length = group_df["end"] - group_df["start"]
        x = group_df["start"] + length / 2 + f_start_pos(group_df["chrom"])
        ax.plot(x, group_df[cn_column], c=color, linewidth=size, label=label, alpha=alpha)
        label = None  # only use label for the first segment
    return ax


def plot_dots(ax, cns_df, cn_column, color="green", label=None, alpha=1, size=1, assembly=hg19):
    length = cns_df["end"] - cns_df["start"]
    f_start_pos = _get_start_vector(assembly)
    x = cns_df["start"] + length / 2 + f_start_pos(cns_df["chrom"])
    ax.scatter(x, cns_df[cn_column], s=size, label=label, color=color, alpha=alpha)
    return ax


def plot_bars(ax, cns_df, cn_column, color="green", label=None, alpha=1, size=1, assembly=hg19):
    length = cns_df["end"] - cns_df["start"]
    f_start_pos = _get_start_vector(assembly)
    x = cns_df["start"] + length / 2 + f_start_pos(cns_df["chrom"])
    ax.bar(x, cns_df[cn_column], width=length, color=color, label=label, alpha=alpha)
    return ax


def plot_heatmap(ax, cns_df, cn_column, min_cn = 0, max_cn = 16, assembly=hg19):
    f_start_pos = _get_start_vector(assembly)
    # lowest value strictly greater than 0
    f_col = _get_CN_color_vector(min_cn, max_cn)

    ax.set_facecolor("gray")
    labels = []
    for i, ((sample_id), group_df) in enumerate(cns_df.groupby("sample_id")):
        y = i
        height = 1
        width = group_df["end"] - group_df["start"]
        left = group_df["start"] + f_start_pos(group_df["chrom"])
        color = f_col(group_df[cn_column])
        ax.barh(y, width, height, left=left, color=color)
        labels.append(sample_id)

    ax.set_yticks(np.arange(len(labels)))
    ax.set_yticklabels(labels)

    ax.autoscale(tight=True)
    return ax


def _get_columns(cns_df, cn_columns):
    if cn_columns == None:
        cn_columns = get_cn_cols(cns_df)
        if len(cn_columns) == 0:
            raise ValueError("If cn_columns is not specified, at least one column ending with '_cn' must exist in data")
    elif isinstance(cn_columns, str):
        if not cn_columns in cns_df.columns:
            raise ValueError("specified CN column must in cns_df.columns")
        cn_columns = [cn_columns]
    elif isinstance(cn_columns, Sequence):
        if len(cn_columns) <= 0:
            raise ValueError("cn_columns must be a string or a non-empty list of strings")
        elif not all(c in cns_df.columns for c in cn_columns):
            raise ValueError("all elements in cn_columns must be columns in data")
    else:
        raise ValueError("cn_columns must be a string or a list of strings")

    return cn_columns


def _get_colors(colors, line_count):
    if colors == None:
        if line_count == 1:
            colors = ["blue"]
        else:
            colors = plt.cm.nipy_spectral(np.linspace(0.05, 0.95, line_count))
    elif line_count == 1:
        colors = [colors]
    elif isinstance(colors, Sequence):
        if line_count != len(colors):
            raise ValueError("colors must be None or a list with the same length as the number of lines")
    else:
        raise ValueError("colors must be None or a list with the same length as the number of lines")
    return colors


def _fig_common(cns_df, f_plot, cn_columns=None, colors=None, size=1, assembly=hg19):
    cn_columns = _get_columns(cns_df, cn_columns)
    groups_df = cns_df.groupby("sample_id")
    line_count = len(groups_df)
    colors = _get_colors(colors, line_count)
    alpha = (1 / line_count) ** (1 / 3) if f_plot == plot_lines else 1 / line_count

    n_columns = len(cn_columns)
    x_min, x_max = x_limits(cns_df, assembly)
    width = max(4, (x_max - x_min) / 200_000_000)
    height = 4*n_columns
    fig, axes = plt.subplots(n_columns, 1, figsize=(width, height), sharex=True)

    for j, cn_column in enumerate(cn_columns):
        ax = axes[j] if n_columns > 1 else axes

        max_cn = cns_df[cn_column].max()
        plot_chr_bg(ax, assembly=assembly, y_min = 0, y_max=max_cn + 1, alpha=0.2)
        for i, (group_key, group_df) in enumerate(groups_df):
            color = colors[i]
            label = group_key
            f_plot(
                ax = ax,
                cns_df = group_df,
                cn_column=cn_column,
                color=color,
                label=label,
                size=size,
                alpha=alpha,
                assembly=assembly,
            )

        # if 1 < line_count <= 3:
        ax.legend(loc="upper right")
        # elif line_count > 3:
        #     ax.legend(loc="upper left", bbox_to_anchor=(1.0, 1.0))

        ax.set_ylabel(f"{cn_column}")
        ax.set_xlim(x_min, x_max)
        if j == n_columns - 1:
            plot_x_ticks(ax, assembly, x_min, x_max)   
            ax.set_xlabel("position on the linear genome")
    
    plt.tight_layout()
    return fig, axes


def fig_lines(cns_df, cn_columns=None, colors=None, size=1, assembly=hg19):
    return _fig_common(cns_df, plot_lines, cn_columns, colors, size, assembly)


def fig_dots(cns_df, cn_columns=None, colors=None, size=1, assembly=hg19):
    return _fig_common(cns_df, plot_dots, cn_columns, colors, size, assembly)


def fig_bars(cns_df, cn_columns=None, colors=None, assembly=hg19):
    return _fig_common(cns_df, plot_bars, cn_columns, colors, 1, assembly)
 

def _make_layout(width, height, n_columns):
    vertical = width > height
    if vertical:
        height = height * n_columns
        sharex, sharey = True, False
        n_rows, n_cols = n_columns, 1
    else:
        width = width * n_columns
        sharex, sharey = False, True
        n_rows, n_cols = 1, n_columns
    width = max(3, width)
    height = max(3, height)
    return plt.subplots(n_rows, n_cols, figsize=(width, height), sharex=sharex, sharey=sharey)


def fig_heatmap(cns_df, cn_columns=None, min_cn = 0, max_cn = 10, vertical = True, assembly=hg19):
    cn_columns = _get_columns(cns_df, cn_columns)

    sample_count = len(cns_df["sample_id"].unique())
    n_columns = len(cn_columns)
    x_min, x_max = x_limits(cns_df, assembly)
    width = (x_max - x_min) / 200_000_000
    height = sample_count / 5
    fig, axes = _make_layout(width, height, n_columns)
    
    min_cn = cns_df[cn_columns][cns_df[cn_columns] > 0].min().min()
    max_cn = min(cns_df[cn_columns].max().max(), max_cn)

    for j, cn_column in enumerate(cn_columns):
        ax = axes[j] if n_columns > 1 else axes
        plot_heatmap(
            ax = ax,
            cns_df = cns_df,
            cn_column=cn_column,
            min_cn = min_cn,
            max_cn = max_cn,
            assembly=assembly,
        )
        if vertical:
            ax.set_ylabel(f"{cn_column}")
        elif j == 0:
            ax.set_ylabel("")
        if not vertical:
            ax.set_xlabel(f"{cn_column}")
        elif j == n_columns - 1:
            ax.set_xlabel(f"position on the linear genome")
        plot_x_lines(ax, assembly)
        plot_x_ticks(ax, assembly, x_min, x_max)   
        ax.set_xlim(x_min, x_max)
    
    ax.margins(x=0, y=0)

    # Add legend
    handles = []
    handles.append(mpatches.Patch(facecolor ='blue', label=f'{max_cn:.2f}', edgecolor='black'))
    handles.append(mpatches.Patch(facecolor ='white', label=f'{min_cn:.2f}', edgecolor='black'))
    handles.append(mpatches.Patch(facecolor ='red', label='0', edgecolor='black'))
    handles.append(mpatches.Patch(facecolor ='gray', label='NaN', edgecolor='black'))
    axes[0].legend(handles=handles, loc='upper left', bbox_to_anchor=(1, 1))
    
    return fig, axes


def x_limits(cns_df, assembly=hg19):
    offset = _get_start_vector(assembly)(cns_df["chrom"])
    min_x = (cns_df["start"] + offset).min()
    max_x = (cns_df["end"] + offset).max()
    return min_x, max_x


def y_limits(cns_df, column):
    min_y = cns_df[column].min()
    max_y = cns_df[column].max()
    return min_y, max_y


def _plot_rectangles(ax, items, y_min, y_max, assembly, color_func, alpha):
    height = y_max - y_min

    for item in items:
        chrom, start, end = item[:3]
        x_pos = start + assembly.chr_starts[chrom]
        color = color_func(item)
        width = end - start
        rect = mpatches.Rectangle((x_pos, y_min), width, height, color=color, alpha=alpha)
        ax.add_patch(rect)

    ax.set_ylim(y_min, y_max)
    ax.set_xlim(0, assembly.gen_len)


def plot_chr_bg(ax, y_min=0, y_max=2, assembly=hg19, alpha=0.2, colored=False):
    def color_func(item):
        chrom = item[0]
        is_even = item[1] % 2 == 0
        return assembly.chr_colors[chrom] if colored else ("darkgray" if is_even else "lightgray")

    items = [(chrom, i, i + length) for i, (chrom, length) in enumerate(assembly.chr_lens.items())]
    _plot_rectangles(ax, items, y_min, y_max, assembly, color_func, alpha=alpha)


def plot_cytobands(ax, y_min=0, y_max=2, assembly=hg19, alpha=0.2, color=None):
    f_color = lambda item: color if color is not None else cytoband_color[item[4]]
    _plot_rectangles(ax, assembly.cytobands, y_min, y_max, assembly, f_color, alpha)


def plot_gaps(ax, y_min=0, y_max=2, assembly=hg19, alpha=0.2, color=None):
    f_color = lambda item: color if color is not None else gap_color[item[3]]
    _plot_rectangles(ax, assembly.gaps, y_min, y_max, assembly, f_color, alpha)


def plot_x_ticks(ax, assembly=hg19, min_x=0, max_x=None):
    positions = list(assembly.chr_lens.items())
    if max_x is None:
        max_x = assembly.aut_len

    x_pos = 0
    major_tick_pos = []
    minor_tick_pos = []
    minor_tick_labels = []
    for chrom, length in positions:
        if min_x <= x_pos <= max_x:
            label_text = "\n" + chrom[3:]
            major_tick_pos.append(x_pos)
            minor_tick_pos.append(x_pos + length / 2)
            minor_tick_labels.append(label_text)
        # Update the x position for the next chromosome
        x_pos += length

    ax.set_xticks(major_tick_pos)
    ax.set_xticks(minor_tick_pos, minor=True)
    ax.set_xticklabels([" "] * len(major_tick_pos))
    ax.set_xticklabels(minor_tick_labels, minor=True)
    ax.set_xlim(min_x, max_x)
        # Hide the lines for the minor ticks
    ax.tick_params(axis='x', which='minor', length=0, pad=-8, labelsize=10)

    return major_tick_pos, minor_tick_pos


def add_cytoband_legend(ax):
    legend_elements = [mpatches.Patch(color=color, label=gie_stain) for gie_stain, color in cytoband_color.items()]
    ax.legend(handles=legend_elements, bbox_to_anchor=(1, 1), loc="upper left")


def add_gap_legend(ax):
    legend_elements = [mpatches.Patch(color=color, label=gap_type) for gap_type, color in gap_color.items()]
    ax.legend(handles=legend_elements, bbox_to_anchor=(1, 1), loc="upper left")


def no_y_ticks(ax):
    ax.set_yticks([])
    ax.set_yticklabels([])


def no_x_ticks(ax):
    ax.set_xticks([])
    ax.set_xticklabels([])


def plot_x_lines(ax, assembly=hg19, positions=None, width=1, alpha=0.5):
    positions = list(accumulate(assembly.chr_lens.values()))
    for pos in positions:
        ax.axvline(
            pos,
            color="black",
            linewidth=width,
            linestyle="--",
            alpha=alpha,
        )
