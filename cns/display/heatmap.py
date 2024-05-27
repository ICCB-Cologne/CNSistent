import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np

from cns.display.label import plot_x_ticks, plot_x_lines, get_size_and_bounds
from cns.utils.conversions import chrom_to_sortable, column_to_label


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


def _plot_heatmap_legend(fig, vpos, top_size, cn_column, max_cn):
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



def _get_x_tick_pos(ax, chr_lens, fontsize):
    x_pos=0
    tick_pos = [x_pos]
    for chrom, length in chr_lens:
        label_text = "\n" + chrom[3:]
        ax.text(
            x_pos + length / 2,
            ax.get_ylim()[0],
            label_text,
            ha="center",
            va="center_baseline",
            fontsize=fontsize,
            linespacing=1.25,
        )
        # Update the x position for the next chromosome
        x_pos += length
        tick_pos.append(x_pos)




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
    sample_ids = cn_bins["sample_id"].unique()
    height = (ratio * len(sample_ids) * width) + 1  
    fig = plt.figure(figsize=(width, height), dpi=dpi)
    ax = fig.gca()

    im_data = plot_CN_heatmap(ax, cn_bins, cn_column, ratio, max_cn)
    pixels_per_row = im_data.shape[0] // len(sample_ids)

    # Create y-ticks for each sample
    _add_y_ticks(ax, sample_ids, pixels_per_row, width)
    ax.set_ylabel("sample")


    pos_per_chr = cn_bins.groupby("chrom")["cum_mid"].nunique()
    chrom_sizes = pos_per_chr.sort_index(key=lambda x: x.map(chrom_to_sortable)).items()
    pos = plot_x_ticks(ax, positions=chrom_sizes)
    plot_x_lines(ax, positions=pos)
    ax.set_xlabel("chromosome")

    # Add a color bar for the heatmap
    fig_size, _, bound_box = get_size_and_bounds(fig, ax, print_info)
    top_size = 1.5 * fig_size[0] / (fig_size[1] * dpi)
    top_pos = bound_box[3] + vertical_legend_correction
    if print_info:
        print("image_pixes: ", im_data.shape)
        print("dpi: ", dpi)
        print("Top size: ", top_size)

    _plot_heatmap_legend(fig, top_pos, top_size, cn_column, max_cn)

    # Add title to the top-left corner of the figure
    fig.text(
        0.10,
        top_pos + top_size,
        label,
        verticalalignment="center",
        horizontalalignment="left",
        fontdict={"fontsize": width * 1.5, "fontweight": "bold"},
    )
