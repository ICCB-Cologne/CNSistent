import numpy as np
import pandas as pd
from cns.utils.assemblies import hg19
from cns.utils.conversions import breaks_to_segments, cytobands_to_df


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
    return { chrom: create_step_breaks(length, step_size, equidistant) for chrom, length in assembly.chr_lens.items() }


def calc_arm_breaks(assembly=hg19):
    cyto_df = cytobands_to_df(assembly.cytobands)
    acen = cyto_df.query("stain == 'acen' and name.str.contains('p')", engine="python")                                
    max_ends = cyto_df.groupby("chrom")["end"].max().to_dict()
    result = { row["chrom"]: [0, row["end"], max_ends[row["chrom"]]] for _, row in acen.iterrows() }
    return result


# all the breakpoints around cytobands
def calc_cytoband_breaks(assembly=hg19):
    cyto_df = cytobands_to_df(assembly.cytobands)
    return { chrom: [0] + [end for end in cyto_df.query(f"chrom == '{chrom}'")["end"]] 
            for chrom in cyto_df["chrom"].unique() }

# TODO: Make breaks should be a dict
# Create breakpoints
def make_breaks(break_type, assembly=hg19):
    if break_type == "arms":
        return calc_arm_breaks(assembly)
    elif break_type == "cytobands":
        return calc_cytoband_breaks(assembly)
    else:
        try:
            bin_size = int(break_type)
        except:
            raise ValueError("break_type must be 'arms', 'cytobands' or an integer, got " + break_type)
        return calc_bin_breaks(bin_size, equidistant=True, assembly=assembly)
    

# Obtain breakpoints from a cns dataframe
def get_breaks(cns_df, keep_ends=True, assembly=hg19):
    breaks_start = cns_df[["chrom", "start"]].copy().drop_duplicates().rename(columns={"start": "break"})
    breaks_end = cns_df[["chrom", "end"]].copy().drop_duplicates().rename(columns={"end": "break"})
    breaks = pd.concat([breaks_start, breaks_end]).drop_duplicates()
    if not keep_ends:
        breaks = breaks[breaks["break"] != breaks["chrom"].map(assembly.chr_lens)]
        breaks = breaks[breaks["break"] > 0]
    breaks.sort_values(by=['chrom', 'break'], inplace=True)
    dict_start = breaks.groupby('chrom')['break'].agg(list).to_dict()
    # insert empty for chromosomes without breakpoints
    for chrom in assembly.chr_names:
        if chrom not in dict_start:
            dict_start[chrom] = []
    return dict_start


def get_breaks_in_segments(segments, breaks):
    res = { chrom: [] for chrom in breaks }
    for chrom, start, end in segments:
        if chrom not in breaks:
            continue
        for br in breaks[chrom]:
            if start <= br < end:
                res[chrom].append(br)
    return res


def insert_breaks_in_segments(segments, breaks):
    breaks = get_breaks_in_segments(segments, breaks)
    # add breaks from segments
    for chrom, start, end in segments:
        if chrom not in breaks:
            breaks[chrom] = []
        breaks[chrom].append(start)
        breaks[chrom].append(end)
    # unique and sort
    for chrom in breaks:
        breaks[chrom] = sorted(set(breaks[chrom]))
    # breaks to segments
    res = breaks_to_segments(breaks)
    return res