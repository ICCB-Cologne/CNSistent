import numpy as np
import pandas as pd
from matplotlib import pyplot as plt, colors as mcolors
from collections.abc import Sequence
from .label import no_y_ticks, plot_chr_bg, plot_x_ticks, plot_x_lines, get_size_and_bounds
from cns.utils.canonization import get_cn_cols
from cns.utils.conversions import chrom_to_sortable
from cns.utils.assemblies import hg19


def plot_lines(ax, cns_df, column, color="green", label=None, alpha=1, line_width=1, chrom=None, assembly=hg19):
    single = chrom is not None
    chroms = [chrom] if single else cns_df["chrom"].unique()
    for chrom in chroms:
        chr_start = 0 if single else assembly.chr_starts[chrom]
        chr_cns_df = cns_df.query(f"chrom == '{chrom}'")
        is_consecutive = chr_cns_df["start"] - chr_cns_df["end"].shift(1) != 0
        # plot consecutive segments
        for _, group_df in chr_cns_df.groupby(is_consecutive.cumsum()):
            length = group_df["end"] - group_df["start"]
            x = group_df["start"] + length / 2 + chr_start
            ax.plot(x, group_df[column], c=color, linewidth=line_width, label=label, alpha=alpha)
            label = None  # only use label for the first chromosome
    return ax


def plot_dots(ax, cns_df, column, color="green", label=None, alpha=1, dot_size=1, chrom=None, assembly=hg19):
    single = chrom is not None
    chroms = [chrom] if single else cns_df["chrom"].unique()
    for chrom in chroms:
        chr_start = 0 if single else assembly.chr_starts[chrom]
        group_df = cns_df.query(f"chrom == '{chrom}'")
        length = group_df["end"] - group_df["start"]
        x = group_df["start"] + length / 2 + chr_start
        ax.scatter(x, group_df[column], s=dot_size, label=label, color=color, alpha=alpha)
        label = None  # only use label for the first chromosome
    return ax


def plot_bars(ax, cns_df, column, color="green", label=None, alpha=1, chrom=None, assembly=hg19):
    single = chrom is not None
    chroms = [chrom] if single else cns_df["chrom"].unique()
    for chrom in chroms:
        chr_start = 0 if single else assembly.chr_starts[chrom]
        group_df = cns_df.query(f"chrom == '{chrom}'")
        length = group_df["end"] - group_df["start"]
        x = group_df["start"] + length / 2 + chr_start
        print(x, length)
        ax.bar(x, group_df[column], width=length, color=color, label=label, alpha=alpha)
        label = None  # only use label for the first chromosome
    return ax



def _check_fig_input(cns_df, cn_columns, label, chrom, assembly):
    if chrom != None and not (isinstance(chrom, str) or not chrom in assembly.keys()):
        raise ValueError("chrom must be None or a string")
    
    if cn_columns == None:
        cn_columns = get_cn_cols(cns_df[0])
        if len(cn_columns) == 0:
            raise ValueError("If cn_columns is not specified, at least one column ending with '_cn' must exist in data")
    elif isinstance(cn_columns, str):
        if not cn_columns in cns_df[0].columns:
            raise ValueError("specified CN column must in cns_df.columns")
        cn_columns = [cn_columns]	
    elif isinstance(cn_columns, Sequence):
        if len(cn_columns) <= 0:
            raise ValueError("cn_columns must be a string or a non-empty list of strings")
        elif not all(c in cns_df[0].columns for c in cn_columns):
            raise ValueError("all elements in cn_columns must be columns in data")
    else:
        raise ValueError("cn_columns must be a string or a list of strings")	
    
    line_count = len(cns_df)*len(cn_columns)
    return cns_df, label, cn_columns, line_count
    

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


def _fig_main(cns_df, plot_func, cn_columns=None, color=None, chrom=None, width=None, assembly=hg19):
    width = width if width != None else (18 if chrom == None else 4)
    height = width / 6 if chrom == None else width
    fig, ax = plt.subplots(1, figsize=(width, height))
    sample_df, columns, line_count = _check_fig_input(cns_df, cn_columns, chrom, assembly)
    colors = _get_colors(color, line_count)
    alpha = (1 / line_count) ** (1/3) if plot_func == plot_lines else 1 / line_count
    if chrom == None:
        min_cn, max_cn = _get_min_max_cn(sample_df, columns) 
        plot_chr_bg(ax, assembly, -0.05, max_cn * 1.05)
        plot_x_ticks(ax, assembly=assembly)
    for i in range(len(sample_df)):
        for j in range(len(columns)):
            color = colors[i*len(columns) + j]
            label = labels[i]
            if len(columns) > 1:
                label += " - " + columns[j]
            cn_columns = columns[j]
            plot_func(ax, sample_df[i], column=cn_columns, color=color, label=label, chrom=chrom, alpha=alpha)
    if has_label: 
        if line_count > 3:
            ax.legend(loc='upper left', bbox_to_anchor=(1.0, 1.0))
        else:
            ax.legend(loc='upper right')
    ax.set_ylabel("mean CN per segment")
    ax.set_xlabel(f"position on {chrom if chrom != None else 'linear genome'}")
    return fig, ax


def fig_lines(cns_df, cn_columns=None, color=None, chrom=None, width=None, assembly=hg19):
    return _fig_main(cns_df, plot_lines, cn_columns, color, chrom, width, assembly)


def fig_dots(cns_df, cn_columns=None, color=None, chrom=None, width=None, assembly=hg19):
    return _fig_main(cns_df, plot_dots, cn_columns, color, chrom, width, assembly)


def fig_bars(cns_df, cn_columns=None, color=None, chrom=None, width=None, assembly=hg19):
    return _fig_main(cns_df, plot_bars, cn_columns, color, chrom, width, assembly)


def fig_heatmap(data_df, cn_columns=None, chrom=None, width=None, dpi=100, assembly=hg19):
    width = width if width != None else (18 if chrom == None else 4)
    height = width / 6 if chrom == None else width
    fig, ax = plt.subplots(1, figsize=(width, height), dpi=dpi)
    dfs, labels, columns, line_count, has_label = _check_fig_input(data_df, cn_columns,  chrom, assembly)


    if chrom is not None:
        data_df = data_df[data_df['chrom'] == chrom]

    samples = data_df['sample_id'].unique()
    n_samples = len(samples)
    n_columns = len(cn_columns)

    fig_height = n_samples * 0.5  # Adjust as needed
    fig_width = n_columns * 5 if width is None else width
    fig, axes = plt.subplots(nrows=1, ncols=n_columns, figsize=(fig_width, fig_height), dpi=dpi)

    # Ensure axes is a 2D array
    if n_samples == 1 and n_columns == 1:
        axes = np.array([[axes]])
    elif n_samples == 1:
        axes = np.array([axes])
    elif n_columns == 1:
        axes = axes[:, np.newaxis]

    cmap = plt.cm.Blues

    for i, sample in enumerate(samples):
        sample_data = data_df[data_df['sample_id'] == sample].copy()
        if chrom is not None:
            genome_length = assembly.chr_lens[chrom]
        else:
            genome_length = sum(assembly.chr_lens.values())

        sample_data['length'] = sample_data['end'] - sample_data['start']
        sample_data.sort_values(['chrom', 'start'], inplace=True)

        for j, col in enumerate(cn_columns):
            ax = axes[i, j]
            segments = sample_data.copy()
            values = segments[col]
            vmin = segments[col].min()
            vmax = segments[col].max()
            norm = mcolors.Normalize(vmin=vmin, vmax=vmax)

            for idx, row in segments.iterrows():
                start = row['start']
                length = row['length']
                value = row[col]
                if pd.isna(value):
                    color = 'white'
                elif value == 0:
                    color = 'red'
                else:
                    color = cmap(norm(value))
                ax.barh(y=0, width=length, left=start, color=color, edgecolor=None, height=0.8)

            ax.set_xlim(0, genome_length)
            ax.set_ylim(-0.5, 0.5)
            if i == n_samples - 1:
                ax.set_xlabel(col)
            if j == 0:
                ax.set_ylabel(sample)
            else:
                ax.set_yticks([])
            ax.set_xticks([])
            ax.set_yticks([])

    plt.tight_layout()
    return fig, axes