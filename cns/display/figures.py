from collections.abc import Sequence

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

from cns.analyze.labels import plot_chr_bg, plot_x_ticks
from cns.utils.assemblies import hg19


def plot_line(ax, grouped, column, color="red", label=None, alpha=1, line_width=1, chrom=None):
	chroms = grouped["chrom"].unique() if chrom is None else [chrom]
	for chr in chroms:
		df = grouped.query(f"chrom == '{chr}'").copy()
		df["is_consecutive"] = df["start"] - df["end"].shift(1) != 0
		# plot consecutive segments
		for _, group_df in df.groupby(df["is_consecutive"].cumsum()):
			x = group_df["cum_mid" if chrom is None else "mid"]
			ax.plot(x, group_df[column], c=color, linewidth=line_width, label=label, alpha=alpha,)
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

	if isinstance(column, str):
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
	return min_cn - min_cn/10, max_cn + max_cn/10


def _get_colors(colors, line_count):	
	if colors == None:
		colors = plt.cm.nipy_spectral(np.linspace(0.05, 0.95, line_count))
	elif line_count == 1:
		colors = [colors]
	elif line_count != len(colors):
		raise ValueError("colors must be None or a list with the same length as the number of lines")
	return colors


def fig_line(
	data,
	assembly=hg19,
	label=None,
	column="total_cn",
	color=None,
	chrom=None,
	width=18,
	dpi=100,
):
	height = width / 6 if chrom == None else width
	fig, ax = plt.subplots(1, figsize=(width, height), dpi=dpi)
	dfs, labels, columns, line_count, has_label = _check_fig_input(data, column, label, chrom, assembly)
	min_cn, max_cn = _get_min_max_cn(dfs, columns)	
	colors = _get_colors(color, line_count)
	plot_chr_bg(ax, assembly, min_cn, max_cn)
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
			label = labels[i] if has_label else None
			plot_line(ax, dfs[i], columns[j], color, label, chrom)
	if has_label:
		ax.legend(loc="upper right")
	fontdict = {"fontsize": width}
	ax.set_ylabel("mean CN per bin", fontdict=fontdict)
	ax.set_xlabel("position on a chromosome", fontdict=fontdict, labelpad=width / 2)
	return fig, ax
