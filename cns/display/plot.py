from collections.abc import Sequence

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt, colors as mcolors

from cns.display.label import plot_chr_bg, plot_x_ticks
from cns.utils.conversions import column_to_label, chrom_to_sortable
from cns.utils.assemblies import hg19


def plot_lines(ax, grouped, column, color="green", label=None, alpha=1, line_width=1, chrom=None):
	chroms = grouped["chrom"].unique() if chrom is None else [chrom]
	for chr in chroms:
		df = grouped.query(f"chrom == '{chr}'").copy()
		df["is_consecutive"] = df["start"] - df["end"].shift(1) != 0
		# plot consecutive segments
		for _, group_df in df.groupby(df["is_consecutive"].cumsum()):
			x = group_df["cum_mid" if chrom is None else "mid"]
			ax.plot(x, group_df[column], c=color, linewidth=line_width, label=label, alpha=alpha)
			label = None  # only use label for the first chromosome
	return ax


def plot_dots(ax, grouped, column, color="green", label=None, alpha=1, dot_size=1, chrom=None):
	chroms = grouped["chrom"].unique() if chrom is None else [chrom]
	for chr in chroms:		
		df = grouped.query(f"chrom == '{chr}'")
		ax.scatter(df["cum_mid"], df[column], s=dot_size, label=label, color=color, alpha=alpha)
		label = None  # only use label for the first chromosome
	return ax


def _check_fig_input(data, column, label, chrom, assembly):
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
		if label == None:
			has_label = False
		elif isinstance(label, Sequence) and len(label) == len(data) and all(isinstance(l, str) for l in label):
			has_label = True
		else:
			raise ValueError("label must be None or a list of strings with the same length as data")
	else:
		raise ValueError("data must be a pandas DataFrame or a list of pandas DataFrames")

	if column == None:
		# set column to all columns in data that end with "_cn"
		column = [c for c in data[0].columns if c.endswith("cn") or c.endswith("CN") or c.startswith("cn") or c.startswith("CN")]
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
		colors = plt.cm.nipy_spectral(np.linspace(0.05, 0.95, line_count))
	elif line_count == 1:
		colors = [colors]
	elif line_count != len(colors):
		raise ValueError("colors must be None or a list with the same length as the number of lines")
	return colors


def _fig_main(data, plot_func, label=None, column=None, color=None, chrom=None, width=18, dpi=100, assembly=hg19):
    height = width / 6 if chrom == None else width
    fig, ax = plt.subplots(1, figsize=(width, height), dpi=dpi)
    dfs, labels, columns, line_count, has_label = _check_fig_input(data, column, label, chrom, assembly)
    min_cn, max_cn = _get_min_max_cn(dfs, columns)    
    colors = _get_colors(color, line_count)
    alpha = (1 / line_count) ** (1/3) if plot_func == plot_lines else 1 / line_count
    plot_chr_bg(ax, assembly, min_cn * .95, max_cn * 1.05)
    plot_x_ticks(ax, assembly)
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
        if len(columns) > 1:
            ax.legend(loc='upper left', bbox_to_anchor=(1.0, 1.0))
        else:
            ax.legend(loc='upper right')
    ax.set_ylabel("mean CN per bin")
    ax.set_xlabel("position on a chromosome")
    return fig, ax


def fig_lines(data, label=None, column=None, color=None, chrom=None, width=18, dpi=100, assembly=hg19):
    return _fig_main(data, plot_lines, label, column, color, chrom, width, dpi, assembly)


def fig_dots(data, label=None, column=None, color=None, chrom=None, width=18, dpi=100, assembly=hg19):
    return _fig_main(data, plot_dots, label, column, color, chrom, width, dpi, assembly)


# Define the colors for the custom heatmap
__heatmap_colors__ = ["red", "white", "blue"]


def _get_table(cn_bins, cn_column="total_cn"):
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


def _plot_legend(fig, vpos, top_size, cn_column, max_cn):
    Z_dummy = np.linspace(0, 1, 10).reshape(1, 10)
    width, height = fig.get_size_inches()
    ax_dummy = fig.add_axes(
        [0.5, vpos + top_size * 2, 0.4, top_size]
    )  # [left, bottom, width, height]
    cmap = mcolors.LinearSegmentedColormap.from_list("custom_colormap", __heatmap_colors__[1:])
    img = ax_dummy.imshow(Z_dummy, cmap=cmap, aspect="auto")
    cbar = fig.colorbar(img, cax=ax_dummy, orientation="horizontal")

    # Add a label to the color bar
    ax_dummy.yaxis.set_label_position("left")
    font = {"fontsize": width, "fontweight": "bold"}
    ax_dummy.set_xlabel(
        column_to_label(cn_column), rotation=0, labelpad=-0.1, va="top", fontdict=font
    )

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
        label.set_fontsize(width * 0.8)


def _add_y_ticks(ax, sample_ids, pixels_per_row, width):
    tick_pos = np.arange(
        pixels_per_row / 2, (len(sample_ids)) * pixels_per_row + pixels_per_row / 2, pixels_per_row
    )
    ax.set_yticks(tick_pos, sample_ids, fontsize=width * 0.8)
    ax.tick_params(axis="y", length=0)  # Hide y-tick marks
    return tick_pos


def plot_CN_heatmap(ax, cn_bins, cn_column="total_cn", ratio=0.02, max_cn=10):
    data = _get_table(cn_bins, cn_column)
    pixels_per_row = max(1, int(np.round(data.shape[1] * ratio)))  
    im_data = _get_im_data(data, max_cn, pixels_per_row)    
    cmap = mcolors.LinearSegmentedColormap.from_list("custom_colormap", __heatmap_colors__)
    ax.imshow(im_data, cmap=cmap, vmin=-1, vmax=1, interpolation="none", origin="lower", extent=[0, im_data.shape[1], 0, im_data.shape[0]])
    return im_data


# for really long plots, the top legend is off by about the correction factor
# I could not figure out why, so I added a correction factor of .25%
# It has almost no effect on smaller plots
# 
# ratio = the height / width of a single sample
def fig_CN_heatmap(
    cn_bins,
    label="CN Tracks",
    cn_column="total_cn",
    width=10,
    max_cn=10,
    ratio=0.02,
    dpi=100,
    print_info=False,
    vertical_legend_correction = 0.0025
):
    sample_ids = cn_bins.index.unique()
    height = (ratio * len(sample_ids) * width) + 1  
    fig = plt.figure(figsize=(width, height), dpi=dpi)
    ax = fig.gca()

    im_data = plot_CN_heatmap(ax, cn_bins, cn_column, ratio, max_cn)
    pixels_per_row = im_data.shape[0] // len(sample_ids)

    # Create y-ticks for each sample
    _add_y_ticks(ax, sample_ids, pixels_per_row, width)
    ax.set_ylabel("sample", fontdict={"fontsize": width})

    # Create x-ticks for each chromosome
    pos_per_chr = cn_bins.groupby("chrom")["cum_mid"].nunique()
    chrom_sizes = pos_per_chr.sort_index(key=lambda x: x.map(chrom_to_sortable)).items()
    x_tick_pos =  plot_x_ticks(ax)

    # add vertical dashed line at each tick pos
    plot_x_lines(ax, x_tick_pos, width / len(x_tick_pos), alpha=0.5)
    ax.set_xlabel("chromosome", fontdict={"fontsize": width}, labelpad=5)

    # Add a color bar for the heatmap
    fig_size, _, bound_box = get_pos(fig, ax, print_info)
    top_size = 1.5 * fig_size[0] / (fig_size[1] * dpi)
    top_pos = bound_box[3] + vertical_legend_correction
    if print_info:
        print("image_pixes: ", im_data.shape)
        print("dpi: ", dpi)
        print("Top size: ", top_size)

    _plot_legend(fig, top_pos, top_size, cn_column, max_cn)

    # Add title to the top-left corner of the figure
    fig.text(
        0.10,
        top_pos + top_size,
        label,
        verticalalignment="center",
        horizontalalignment="left",
        fontdict={"fontsize": width * 1.5, "fontweight": "bold"},
    )
