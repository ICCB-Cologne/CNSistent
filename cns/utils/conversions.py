import numpy as np
import pandas as pd


def cytobands_to_df(cytobands):
    return pd.DataFrame(cytobands, columns=["chrom", "start", "end", "name", "stain"])


def gaps_to_df(gaps):
    return pd.DataFrame( gaps, columns=["chrom", "start", "end", "type", "bridge"])


def sampleid_to_numeric_df(cna_df):
    cna_df = cna_df.copy()
    samples = cna_df["sample_id"].unique().tolist()
    # replace sampleid with numeric values
    sample_dict = {samples[i]: i for i in range(len(samples))}
    cna_df["sample_id"] = cna_df["sample_id"].replace(sample_dict)
    return samples, cna_df


def numeric_to_sampleid_df(samples, cna_df):
    cna_df = cna_df.copy()
    sample_dict = {i: samples[i] for i in range(len(samples))}
    cna_df["sample_id"] = cna_df["sample_id"].replace(sample_dict)
    return cna_df


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


def column_to_label(cn_column):
    if cn_column == "total_cn":
        return "Total CN"
    elif cn_column == "major_cn":
        return "Major CN"
    elif cn_column == "minor_cn":
        return "Minor CN"
    else:
        return cn_column


def segs_to_chrom_dict(segments):
    res = { }
    for seg in segments:
        chrom, start, end = seg
        if chrom not in res:
            res[chrom] = []
        res[chrom].append((start, end))
    return res
