
from matplotlib.patches import Rectangle
import matplotlib.patches as mpatches
from cns.utils.conversions import *

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


def plot_chr_bg(ax, assembly, y_min=0, y_max=1, colored=False):
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


def plot_cytobands(ax, assembly, y_min=0, y_max=1, alpha=.2, sel_chrom=None):
    for band in assembly.cytobands:
        chrom, start, end, name, gie_stain = band
        if (sel_chrom is not None) and (sel_chrom != chrom):
            continue
        x_pos = (start + assembly.chr_starts[chrom]) if sel_chrom is None else start
        color = cytoband_color(gie_stain)
        length = end - start + 1
        rect = Rectangle(
            (x_pos, y_min), length, y_max-y_min, color=color, alpha=alpha
        )
        ax.add_patch(rect)
   
    ax.set_ylim(y_min, y_max)     
    ax.set_xlim(0, assembly.gen_len if sel_chrom is None else assembly.chr_lens[sel_chrom])


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


def plot_gaps(ax, assembly, y_min=0, y_max=1, alpha=.2, sel_chrom=None):
    for band in assembly.gaps:
        chrom, start, end, gap_type, bridge = band        
        if (sel_chrom is not None) and (sel_chrom != chrom):
            continue
        x_pos = (start + assembly.chr_starts[chrom]) if sel_chrom is None else start
        color = gap_color(gap_type)
        length = end - start + 1
        rect = Rectangle(
            (x_pos, y_min), length, y_max-y_min, color=color, alpha=alpha
        )
        ax.add_patch(rect)
    
    ax.set_ylim(y_min, y_max)
    ax.set_xlim(0, assembly.gen_len if sel_chrom is None else assembly.chr_lens[sel_chrom])



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


def _plot_x_ticks(ax, chr_lens, fontsize):
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

    ax.set_xticks(tick_pos)
    ax.set_xticklabels([" "] * len(tick_pos))
    ax.set_xlim(0, x_pos)
    return tick_pos


def no_y_ticks(ax):
    ax.set_yticks([])
    ax.set_yticklabels([])
    

def plot_x_ticks(ax, assembly, fontsize=8):
    positions = list(assembly.chr_lens.items())
    return _plot_x_ticks(ax, positions, fontsize)


def _plot_x_lines(ax, positions, width, alpha):
    for pos in positions:
        ax.axvline(
            pos,
            color="black",
            linewidth=width,
            linestyle="--",
            alpha=alpha,
        )


def plot_x_lines(ax, assembly, width=1, alpha=.5):    
    positions = list(assembly.chr_starts.values())[1:]
    _plot_x_lines(ax, positions, width, alpha)