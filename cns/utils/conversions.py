import pandas as pd

from cns.utils.assemblies import hg19


def cytobands_to_df(cytobands):
    return pd.DataFrame(cytobands, columns=["chrom", "start", "end", "name", "stain"])


def gaps_to_df(gaps):
    return pd.DataFrame( gaps, columns=["chrom", "start", "end", "type", "bridge"])


def segs_to_df(segments):
    seg_list = []
    for chrom in sorted(segments.keys()):  # Sort the keys lexicographically
        for seg in segments[chrom]:
            seg_list.append((chrom, seg[0], seg[1]))
    return pd.DataFrame(seg_list, columns=["chrom", "start", "end"])


def df_to_segs(segs_df):
    res = {}
    for chrom, group in segs_df.groupby("chrom"):
        res[chrom] = list(zip(group["start"], group["end"]))
    return res


def chrom_to_sortable(chrom, aut_count = 22):
    if chrom == "chrX":
        return aut_count + 1  # Make 'chrX' sort last
    if chrom == "chrY":
        return aut_count + 2
    if chrom == "chrM":
        return aut_count + 3
    else:
        return int(chrom[3:])  # Remove 'chr' and convert to int
    

def sortable_to_chrom(sortable, aut_count = 22):
    if sortable <= aut_count:
        return "chr" + str(sortable)
    if sortable == aut_count + 1:
        return "chrX"
    if sortable == aut_count + 2:
        return "chrY"
    if sortable == aut_count + 3:
        return "chrM"


def segs_to_chrom_dict(segments):
    res = { }
    for seg in segments:
        chrom, start, end = seg
        if chrom not in res:
            res[chrom] = []
        res[chrom].append((start, end))
    return res


def tuples_to_segments(tuples):
    segs = {}
    if len(tuples) > 0 and len(tuples[0]) >= 3:
        for tuple in tuples:
            if tuple[0] not in segs:
                segs[tuple[0]] = []
            segs[tuple[0]].append((tuple[1], tuple[2]))
    return segs


def breaks_to_segments(breakpoints):
    segs = {}
    for chrom, breaks in breakpoints.items():
        segs[chrom] = []
        for i in range(len(breaks) - 1):
            segs[chrom].append((breaks[i], breaks[i + 1]))
    return segs


def genome_to_segments(assembly=hg19):
    segs = {}
    for chrom, len in assembly.chr_lens.items():
        segs[chrom] = [(0, len)]
    return segs
