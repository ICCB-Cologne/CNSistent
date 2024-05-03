from cns.utils.files import load_cns
import pandas as pd
from os.path import join as pjoin


class Dataset:
    def __init__(self, path, name, color):
        self.name = name
        self.cns = load_cns(pjoin(path, f"{name}_cns_imp.tsv"))
        self.samples = pd.read_csv(pjoin(path, f"{name}_samples.tsv"), sep="\t")
        self.color = color


def load_data(path):
    data = {
        "PCAWG" : Dataset(path, "PCAWG", "#1f77b4"),
        "TCGA": Dataset(path, "TCGA_hg19", "#ff7f0e"),
        "TRACERx": Dataset(path, "TRACERx", "#2ca02c")
    }
    return data