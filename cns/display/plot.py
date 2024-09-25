import numpy as np
import pandas as pd
from matplotlib import pyplot as plt, colors as mcolors
from collections.abc import Sequence

from cns.display.label import no_y_ticks, plot_chr_bg, plot_x_ticks, plot_x_lines, get_size_and_bounds
from cns.utils.canonization import find_cn_cols_if_none
from cns.utils.conversions import chrom_to_sortable
from cns.utils.assemblies import hg19


def plot_lines(ax, grouped, column, color="green", label=None, alpha=1, line_width=1, chrom=None, pos_col=None):
    chroms = grouped["chrom"].unique() if chrom is None else [chrom]
    pos_col = ("cum_mid" if chrom is None else "mid") if pos_col is None else pos_col
    for chr in chroms:
        df = grouped.query(f"chrom == '{chr}'").copy()
        df["is_consecutive"] = df["start"] - df["end"].shift(1) != 0
        # plot consecutive segments
        for _, group_df in df.groupby(df["is_consecutive"].cumsum()):
            x = group_df[pos_col]
            ax.plot(x, group_df[column], c=color, linewidth=line_width, label=label, alpha=alpha)
            label = None  # only use label for the first chromosome
    return ax


def plot_dots(ax, grouped, column, color="green", label=None, alpha=1, dot_size=1, chrom=None, pos_col=None):
    chroms = grouped["chrom"].unique() if chrom is None else [chrom]
    pos_col = ("cum_mid" if chrom is None else "mid") if pos_col is None else pos_col
    for chr in chroms:		
        df = grouped.query(f"chrom == '{chr}'")
        ax.scatter(df[pos_col], df[column], s=dot_size, label=label, color=color, alpha=alpha)
        label = None  # only use label for the first chromosome
    return ax


def plot_bars(ax, grouped, column, color="green", label=None, alpha=1, chrom=None, pos_col=None):
    chroms = grouped["chrom"].unique() if chrom is None else [chrom]
    pos_col = ("cum_mid" if chrom is None else "mid") if pos_col is None else pos_col
    for chrom in chroms:
        df = grouped.query(f"chrom == '{chrom}'")
        ax.bar(df[pos_col], df[column], width=df["length"], color=color, label=label, alpha=alpha)
        label = None  # only use label for the first chromosome
    return ax


def _check_fig_input(data, column, label, chrom, assembly, pos_col):
    if chrom != None and not (isinstance(chrom, str) or not chrom in assembly.keys()):
        raise ValueError("chrom must be None or a string")
    
    if isinstance(data, pd.DataFrame):
        if label == None:
            has_label = False
        elif isinstance(label, str):
            has_label = True
        else:
            raise ValueError("label must be None or a string")
        label = [label]
        data = [data]
    elif isinstance(data, Sequence) and len(data) > 0 and isinstance(data[0], pd.DataFrame):
        label = list(label) # numpy array would cause problems here
        if list(label) == None:
            has_label = False
        elif isinstance(label, Sequence) and len(label) == len(data) and all(isinstance(l, str) for l in label):
            has_label = True
        else:
            raise ValueError("label must be None or a list of strings with the same length as data")
    else:
        raise ValueError("data must be a pandas DataFrame or a list of pandas DataFrames")

    if column == None:
        column = find_cn_cols_if_none(data[0])
        if len(column) == 0:
            raise ValueError("If column is not specified, at least one column ending with '_cn' must exist in data")
    elif isinstance(column, str):
        if not column in data[0].columns:
            raise ValueError("column must be a column in data")
        column = [column]	
    elif isinstance(column, Sequence):
        if len(column) <= 0:
            raise ValueError("column must be a string or a non-empty list of strings")
        elif not all(c in data[0].columns for c in column):
            raise ValueError("all elements in column must be columns in data")
    else:
        raise ValueError("column must be a string or a list of strings")	
    
    for df in data:
        if pos_col not in df.columns:
            raise ValueError(f"all dataframes must have a column '{pos_col}'")
    
    line_count = len(data)*len(column)
    return data, label, column, line_count, has_label
    

def _get_min_max_cn(dfs, columns):
    min_cn = np.inf
    max_cn = -np.inf
    for df in dfs:
        for column in columns:
            min_cn = min(min_cn, df[column].min())
            max_cn = max(max_cn, df[column].max())
    return min_cn, max_cn


def _get_colors(colors, line_count):	
    if colors == None:
        if line_count == 1:
            colors = ["blue"]
        else:            
            colors = plt.cm.nipy_spectral(np.linspace(0.05, 0.95, line_count))
    elif line_count == 1:
        colors = [colors]
    elif line_count != len(colors):
        raise ValueError("colors must be None or a list with the same length as the number of lines")
    return colors


def _fig_main(data_df, plot_func, label=None, column=None, color=None, chrom=None, width=None, dpi=100, assembly=hg19, pos_col="cum_mid"):
    width = width if width != None else (18 if chrom == None else 4)
    height = width / 6 if chrom == None else width
    fig, ax = plt.subplots(1, figsize=(width, height), dpi=dpi)
    dfs, labels, columns, line_count, has_label = _check_fig_input(data_df, column, label, chrom, assembly, pos_col)
    colors = _get_colors(color, line_count)
    alpha = (1 / line_count) ** (1/3) if plot_func == plot_lines else 1 / line_count
    if chrom == None:
        min_cn, max_cn = _get_min_max_cn(dfs, columns) 
        plot_chr_bg(ax, assembly, min_cn * .95, max_cn * 1.05)
        plot_x_ticks(ax, assembly=assembly)
    for i in range(len(dfs)):
        for j in range(len(columns)):
            color = colors[i*len(columns) + j]
            if has_label:
                label = labels[i]
                if len(columns) > 1:
                    label += " - " + columns[j]
            else:
                label = None
            column = columns[j]
            plot_func(ax, dfs[i], column=column, color=color, label=label, chrom=chrom, alpha=alpha)
    if has_label: 
        if line_count > 3:
            ax.legend(loc='upper left', bbox_to_anchor=(1.0, 1.0))
        else:
            ax.legend(loc='upper right')
    ax.set_ylabel("mean CN per segment")
    ax.set_xlabel(f"position on {chrom if chrom != None else 'linear genome'}")
    return fig, ax


def fig_lines(data, label=None, column=None, color=None, chrom=None, width=None, dpi=100, assembly=hg19):
    return _fig_main(data, plot_lines, label, column, color, chrom, width, dpi, assembly)


def fig_dots(data, label=None, column=None, color=None, chrom=None, width=None, dpi=100, assembly=hg19):
    return _fig_main(data, plot_dots, label, column, color, chrom, width, dpi, assembly)


def fig_bars(data, label=None, column=None, color=None, chrom=None, width=None, dpi=100, assembly=hg19):
    return _fig_main(data, plot_bars, label, column, color, chrom, width, dpi, assembly)


# Define the colors for the custom heatmap
__heatmap_colors__ = ["red", "white", "blue"]


def _get_column_table(cn_bins, cn_column="total_cn", chrom=None):
    if chrom is not None:
        cn_bins = cn_bins.query(f"chrom == '{chrom}'")
    sample_heatmap_df = cn_bins[["sample_id", cn_column]].copy()
    sample_heatmap_df["position"] = (
        sample_heatmap_df.groupby("sample_id").cumcount() + 1
    )
    sample_pivot_df = sample_heatmap_df.pivot(
        index="sample_id", columns="position", values=cn_column
    ).fillna(0)
    return sample_pivot_df.values


# data to image
def _get_im_data(data, max_cn, pixels_per_row):
    data = np.clip(data, 0, max_cn) / max_cn
    data[data == 0] = -1
    return np.repeat(data, pixels_per_row, axis=0)


def _plot_heatmap_legend(fig, vpos, top_size, max_cn, font_size=None):
    Z_dummy = np.linspace(0, 1, 10).reshape(1, 10)
    width, height = fig.get_size_inches()
    ax_dummy = fig.add_axes(
        [0.5, vpos + top_size, 0.4, top_size]
    )  # [left, bottom, width, height]
    cmap = mcolors.LinearSegmentedColormap.from_list("custom_colormap", __heatmap_colors__[1:])
    img = ax_dummy.imshow(Z_dummy, cmap=cmap, aspect="auto")
    cbar = fig.colorbar(img, cax=ax_dummy, orientation="horizontal")

    # Add tick marks to the color bar
    tick_vals = list(range(max_cn + 1))
    cbar.set_ticks([val / max_cn for val in tick_vals])
    cbar.set_ticklabels(map(str, tick_vals))

    # Add red vertical line at the zero position of the color bar
    position_zero_in_colorbar = (0 - min(tick_vals)) / (max(tick_vals) - min(tick_vals))
    ax_dummy.axvline(position_zero_in_colorbar, color="red", linewidth=width / 2)

    # Make the zero tick label red, rescale
    for label in ax_dummy.get_xticklabels():
        if label.get_text() == "0":
            label.set_color("red")
        label.set_fontsize(font_size)


def _add_track_y_ticks(ax, sample_ids, pixels_per_row, font_size):
    tick_pos = np.arange(
        pixels_per_row / 2, (len(sample_ids)) * pixels_per_row + pixels_per_row / 2, pixels_per_row
    )
    ax.set_yticks(tick_pos, sample_ids, fontsize=font_size)
    ax.tick_params(axis="y", length=0)  # Hide y-tick marks
    return tick_pos


def plot_CN_heatmap(ax, cn_bins, cn_column="total_cn", chrom=None, ratio=0.02, max_cn=10):
    data = _get_column_table(cn_bins, cn_column, chrom)
    pixels_per_row = max(1, int(np.round(data.shape[1] * ratio)))  
    im_data = _get_im_data(data, max_cn, pixels_per_row)    
    cmap = mcolors.LinearSegmentedColormap.from_list("custom_colormap", __heatmap_colors__)
    ax.imshow(im_data, cmap=cmap, vmin=-1, vmax=1, interpolation="none", origin="lower", extent=[0, im_data.shape[1], 0, im_data.shape[0]])
    return im_data


def fig_CN_heatmap(
    data_df,
    label="CN Tracks",
    column=None,
    chrom=None,
    max_cn=10,
    width=None,
    dpi=100,
    assembly=hg19,	
    vertical_legend_correction=0.0025,
    print_info=False
):
    if chrom != None and not (isinstance(chrom, str) or not chrom in assembly.keys()):
        raise ValueError("chrom must be None or a string")
    data_df = data_df.query(f"chrom == '{chrom}'") if chrom is not None else data_df
    sample_ids = data_df["sample_id"].unique()
    pos_col = "cum_mid" 
    dfs, labels, columns, line_count, has_label = _check_fig_input(data_df, column, label, chrom, assembly, pos_col)
    max_cn = min(max_cn, int(np.ceil(_get_min_max_cn(dfs, columns)[1])))

    width = width if width != None else (18 if chrom == None else 6) 
    ratio= 0.02 if chrom is None else 0.2
    height = ((ratio * len(sample_ids) / len(columns) * width) + 1) 
    font_size = (width * 0.8 if chrom is None else width * 3) / len(columns)
    
    fig, axs = plt.subplots(1, len(columns), figsize=(width, height), dpi=dpi)

    im_datas = []
    for i, column in enumerate(columns):
        ax = axs[i] if len(columns) > 1 else axs

        im_data = plot_CN_heatmap(ax, data_df, column, chrom, ratio, max_cn)
        pixels_per_row = im_data.shape[0] // len(sample_ids)

        # Create y-ticks for each sample
        if i == 0:
            ax.set_ylabel("sample", fontsize=font_size*1.25)
            _add_track_y_ticks(ax, sample_ids, pixels_per_row, font_size)
        else:
            no_y_ticks(ax)

        pos_per_chr = data_df.groupby("chrom")[pos_col].nunique()
        chrom_sizes = pos_per_chr.sort_index(key=lambda x: x.map(chrom_to_sortable)).items()
        pos = plot_x_ticks(ax, positions=chrom_sizes, font_size=font_size)
        plot_x_lines(ax, positions=pos)
        ax.set_xlabel(column, fontsize=font_size*1.25)

        im_datas.append(im_data)

    plt.tight_layout()

    # Add a color bar for the heatmap
    ax = axs[1] if len(columns) > 1 else axs
    fig_size, _, bound_box = get_size_and_bounds(fig, ax, len(columns), print_info)
    top_size = 50 * ratio * fig_size[0] / (fig_size[1] * dpi) / len(columns)
    top_pos = bound_box[3] + vertical_legend_correction
    if print_info:
        print("image_pixes: ", [im_data.shape for im_data in im_datas])
        print("dpi: ", dpi)
        print("Top size: ", top_size)

    _plot_heatmap_legend(fig, top_pos, top_size, max_cn, font_size)

    # Add title to the top-left corner of the figure
    fig.text(
        bound_box[0],
        top_pos + top_size,
        label,
        verticalalignment="center",
        fontdict={"fontsize": font_size*2, "fontweight": "bold", "ha": "left"},
    )

    return fig, axs