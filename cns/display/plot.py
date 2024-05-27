from collections.abc import Sequence

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt, colors as mcolors

from cns.display.label import plot_chr_bg, plot_x_ticks
from cns.utils.files import get_cn_columns
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
		column = get_cn_columns(data[0])
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
		colors = plt.cm.nipy_spectral(np.linspace(0.05, 0.95, line_count))
	elif line_count == 1:
		colors = [colors]
	elif line_count != len(colors):
		raise ValueError("colors must be None or a list with the same length as the number of lines")
	return colors


def _fig_main(data, plot_func, label=None, column=None, color=None, chrom=None, width=18, dpi=100, assembly=hg19, pos_col="cum_mid"):
    height = width / 6 if chrom == None else width
    fig, ax = plt.subplots(1, figsize=(width, height), dpi=dpi)
    dfs, labels, columns, line_count, has_label = _check_fig_input(data, column, label, chrom, assembly, pos_col)
    min_cn, max_cn = _get_min_max_cn(dfs, columns)    
    colors = _get_colors(color, line_count)
    alpha = (1 / line_count) ** (1/3) if plot_func == plot_lines else 1 / line_count
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
    ax.set_ylabel("mean CN per bin")
    ax.set_xlabel("position on a chromosome")
    return fig, ax


def fig_lines(data, label=None, column=None, color=None, chrom=None, width=18, dpi=100, assembly=hg19):
    return _fig_main(data, plot_lines, label, column, color, chrom, width, dpi, assembly)


def fig_dots(data, label=None, column=None, color=None, chrom=None, width=18, dpi=100, assembly=hg19):
    return _fig_main(data, plot_dots, label, column, color, chrom, width, dpi, assembly)