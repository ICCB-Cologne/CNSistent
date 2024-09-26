import pandas as pd

from cns.utils.assemblies import hg19


def cytobands_to_df(cytobands):
    return pd.DataFrame(cytobands, columns=["chrom", "start", "end", "name", "stain"])


def gaps_to_df(gaps):
    return pd.DataFrame(gaps, columns=["chrom", "start", "end", "type", "bridge"])


def segs_to_df(segments):
    seg_list = []
    for chrom in sorted(segments.keys()):  # Sort the keys lexicographically
        for seg in segments[chrom]:
            if len(seg) > 2:
                seg_list.append((chrom, seg[0], seg[1], seg[2]))
            else:
                seg_list.append((chrom, seg[0], seg[1]))
    res_df = pd.DataFrame(seg_list)
    if len(res_df.columns) == 3:
        return res_df.rename(columns={0: "chrom", 1: "start", 2: "end"})
    else:
        return res_df.rename(columns={0: "chrom", 1: "start", 2: "end", 3: "name"})


def df_to_segs(segs_df):
    res = {}
    for chrom, group in segs_df.groupby("chrom"):
        res[chrom] = list(zip(group["start"], group["end"], group["name"]))
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
            if len(tuple) == 3:
                segs[tuple[0]].append((tuple[1], tuple[2]))
            else:
                segs[tuple[0]].append((tuple[1], tuple[2], tuple[3]))
    return segs


def breaks_to_segments(breakpoints):
    segs = {}
    seg_i = 0
    for chrom, breaks in breakpoints.items():
        segs[chrom] = []
        for i in range(len(breaks) - 1):
            segs[chrom].append((breaks[i], breaks[i + 1], seg_i))
            seg_i += 1
    return segs


def genome_to_segments(assembly=hg19):
    segs = {}
    for chrom, len in assembly.chr_lens.items():
        segs[chrom] = [(0, len, chrom)]
    return segs
