from cns.utils.conversions import cytobands_to_df
from cns.utils import hg19
import numpy as np


def create_step_breaks(reg_len, step_size, equidistant=True):
    if (step_size < 1) or (reg_len < step_size):
        return [0, reg_len]
    padding = reg_len % step_size
    if equidistant:
        step_count = reg_len // step_size
        if step_count < 1:
            return [0, reg_len]
        if padding > (step_size / 2):
            step_size -= (step_size - padding) / (step_count + 1)
        else:
            step_size += padding / step_count
        fracs = np.arange(0, reg_len + 1, step_size) ## 
        return [np.int32(np.floor(frac + .5)) for frac in fracs]
    else:
        if padding > (step_size / 2):
            padding = (step_size - padding) / 2
            start = step_size - padding
            end = reg_len - start
            fracs = np.arange(start, end + 1, step_size)
            return [0] + [np.int32(np.floor(frac + .5)) for frac in fracs] + [reg_len]
        else:
            if reg_len // step_size <= 1:
                return [0, reg_len]
            start = step_size + padding / 2
            end = reg_len - start
            fracs = np.arange(start, end + 1, step_size)
            return [0] + [np.int32(np.floor(frac + .5)) for frac in fracs] + [reg_len]
       

# Calculate breakpoints at the given resolution. The boundaries either use half
def calc_bin_breaks(step_size, equidistant=True, assembly=hg19):
    return [
        (chrom, create_step_breaks(length, step_size, equidistant))
        for chrom, length in assembly.chr_lens.items()
    ]


def calc_arm_breaks(assembly=hg19):
    cyto_df = cytobands_to_df(assembly.cytobands)
    acen = cyto_df.query("stain == 'acen' and name.str.contains('p')", engine="python")                                
    max_ends = cyto_df.groupby("chrom")["end"].max().to_dict()
    result = [
        (row["chrom"], [0, row["end"], max_ends[row["chrom"]]])
        for _, row in acen.iterrows()
    ]
    return result


# all the breakpoints around cytobands
def calc_cytoband_breaks(assembly=hg19):
    cyto_df = cytobands_to_df(assembly.cytobands)
    return [
        (chrom, [0] + [end for end in cyto_df.query(f"chrom == '{chrom}'")["end"]])
        for chrom in cyto_df["chrom"].unique()
    ]


def get_breakpoints(break_type, assembly=hg19):
    if break_type == "arms":
        return calc_arm_breaks(assembly)
    elif break_type == "cytobands":
        return calc_cytoband_breaks(assembly)
    else:
        try:
            bin_size = int(break_type)
        except:
            raise ValueError(
                "break_type must be 'arms', 'cytobands' or an integer, got "
                + break_type
            )
        return calc_bin_breaks(bin_size, equidistant=True, assembly=assembly)