
from matplotlib.patches import Rectangle
import matplotlib.patches as mpatches
from cns.utils.assemblies import hg19
from itertools import accumulate

# Get the tight bounding box of the specified ax, including tick labels, in display units
#
# Example print:
# Figure size in inches: [10.         10.17996442] (Width x Height)
# Ax size in inches: [8.41291667 0.63665135] (Width x Height)
# Ax coordinates within the figure:
#    Left: 0.07935416666666673
#    Bottom: 0.46873017812050666
#    Right: 0.9206458333333333
#    Top: 0.5312698218794933

def get_size_and_bounds(fig, ax, columns = 1, should_print=False):
    # Get the tight bounding box of the axes, including tick labels, in display units
    tight_bbox = ax.get_tightbbox(fig.canvas.get_renderer())
    # Transform the bounding box to figure fraction coordinates
    transformed_bbox = tight_bbox.transformed(fig.transFigure.inverted())
    # Calculate the size of the bounding box in inches
    fig_size = fig.get_size_inches()
    if should_print:
        print(f"Figure size in inches: {fig_size} (Width x Height)")
    axes_size_inches = fig_size * transformed_bbox.size
    if should_print:
        print(f"Ax size in inches: {axes_size_inches} (Width x Height)")
    # Calculate the position of the bounding box in figure fraction coordinates
    left_pos = (1 - axes_size_inches[0] / fig_size[0] * columns) / 2
    bottom_pos = (1 - axes_size_inches[1] / fig_size[1]) / 2
    right_position = 1 - left_pos
    top_position = 1 - bottom_pos
    bound_box = [left_pos, bottom_pos, right_position, top_position]
    if should_print:
        print(f"Ax coordinates within the figure:")
        for name, val in zip(["Left", "Bottom", "Right", "Top"], bound_box):
            print(f"\t{name}: {val}")
    return fig_size, axes_size_inches, bound_box


def cytoband_color(gie_stain):
    if gie_stain == "acen":
        return "red"
    elif gie_stain == "gneg":
        return "white"
    elif gie_stain == "gpos25":
        return "lightgray"
    elif gie_stain == "gpos50":
        return "gray"
    elif gie_stain == "gpos75":
        return "darkgray"
    elif gie_stain == "gpos100":
        return "black"
    elif gie_stain == "stalk":
        return "black"
    elif gie_stain == "gvar":
        return "yellow"
    else:
        return "white"
    

# 'centromere', 'clone', 'contig', 'heterochromatin', 'short_arm', 'telomere'
def gap_color(gap_type):
    if gap_type == "centromere":
        return "red"
    elif gap_type == "short_arm":
        return "yellow"
    elif gap_type == "telomere":
        return "blue"
    elif gap_type == "heterochromatin":
        return "gray"
    elif gap_type == "clone":
        return "purple"
    elif gap_type == "contig":
        return "green"
    else:
        return "white"


def plot_chr_bg(ax, assembly=hg19, y_min=0, y_max=1, colored=False):
    x_pos = 0
    is_even = True
    for chrom, length in assembly.chr_lens.items():        
        # Get the color for the current chromosome
        color = (
            assembly.chr_colors[chrom] if colored else ("darkgray" if is_even else "lightgray")
        )
        # Add a rectangle to the plot with the appropriate color and width
        rect = Rectangle((x_pos, y_min), length, y_max - y_min, color=color, alpha=0.2)
        ax.add_patch(rect)
        # Update the x position for the next chromosome
        x_pos += length
        is_even = not is_even

    ax.set_ylim(y_min, y_max)
    ax.set_xlim(0, assembly.gen_len)


def plot_cytobands(ax, assembly=hg19, y_min=0, y_max=1, alpha=.2, chrom=None):
    for band in assembly.cytobands:
        band_chr, start, end, name, gie_stain = band
        if (chrom is not None) and (band_chr != chrom):
            continue
        x_pos = (start + assembly.chr_starts[band_chr]) if chrom is None else start
        color = cytoband_color(gie_stain)
        length = end - start
        rect = Rectangle(
            (x_pos, y_min), length, y_max-y_min, color=color, alpha=alpha
        )
        ax.add_patch(rect)
   
    ax.set_ylim(y_min, y_max)     
    ax.set_xlim(0, assembly.gen_len if chrom is None else assembly.chr_lens[chrom])


def add_cytoband_legend(ax):
    legend_elements = [
        mpatches.Patch(color='red', label='Acen'),
        mpatches.Patch(color='white', label='Gneg'),
        mpatches.Patch(color='lightgray', label='Gpos25'),
        mpatches.Patch(color='gray', label='Gpos50'),
        mpatches.Patch(color='darkgray', label='Gpos75'),
        mpatches.Patch(color='black', label='Gpos100'),
        mpatches.Patch(color='yellow', label='Gvar'),
    ]
    # Add the legend to your ax
    ax.legend(handles=legend_elements, bbox_to_anchor=(1, 1), loc='upper left')


def plot_gaps(ax, assembly=hg19, y_min=0, y_max=1, alpha=.2, chrom=None, color=None):
    for gap in assembly.gaps:
        gap_chrom, start, end, gap_type, bridge = gap        
        if (chrom is not None) and (chrom != gap_chrom):
            continue
        x_pos = (start + assembly.chr_starts[gap_chrom]) if chrom is None else start
        color = gap_color(gap_type) if color is None else color
        length = end - start + 1
        rect = Rectangle(
            (x_pos, y_min), length, y_max-y_min, color=color, alpha=alpha
        )
        ax.add_patch(rect)
    
    ax.set_ylim(y_min, y_max)
    ax.set_xlim(0, assembly.gen_len if chrom is None else assembly.chr_lens[chrom])


def add_gap_legend(ax):
    legend_elements = [
        mpatches.Patch(color='red', label='Centromere'),
        mpatches.Patch(color='yellow', label='Short Arm'),
        mpatches.Patch(color='blue', label='Telomere'),
        mpatches.Patch(color='gray', label='Heterochromatin'),
        mpatches.Patch(color='purple', label='Clone'),
        mpatches.Patch(color='green', label='Contig'),
        mpatches.Patch(color='white', label='Other')
    ]
    # Add the legend to your ax
    ax.legend(handles=legend_elements, bbox_to_anchor=(1, 1), loc='upper left')


def plot_x_ticks(ax, assembly=hg19, positions=None, font_size=None):
    if positions is None:
        positions = list(assembly.chr_lens.items())
    x_pos=0
    tick_pos = [x_pos]
    for chrom, length in positions:
        label_text = "\n" + chrom[3:]
        ax.text(
            x_pos + length / 2,
            ax.get_ylim()[0],
            label_text,
            ha="center",
            va="center_baseline",
            fontsize=font_size,
        )
        # Update the x position for the next chromosome
        x_pos += length
        tick_pos.append(x_pos)

    ax.set_xticks(tick_pos)
    ax.set_xticklabels([" "] * len(tick_pos))
    ax.set_xlim(0, x_pos)
    return tick_pos


def no_y_ticks(ax):
    ax.set_yticks([])
    ax.set_yticklabels([])
    

def plot_x_lines(ax, assembly=hg19, positions=None, width=1, alpha=.5):        
    if positions is None:
        positions = list(accumulate(assembly.chr_lens.values()))
    for pos in positions:
        ax.axvline(
            pos,
            color="black",
            linewidth=width,
            linestyle="--",
            alpha=alpha,
        )