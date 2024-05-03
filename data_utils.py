from cns.utils.files import load_cns
import pandas as pd


class Dataset:
    def __init__(self, name, color):
        self.name = name
        self.cns = load_cns(f"./out/{name}_cna_imp.tsv")
        self.samples = pd.read_csv(f"./out/{name}_samples.tsv", sep="\t")
        self.color = color


def load_data():
    data = {
        "PCAWG" : Dataset("PCAWG", "#1f77b4"),
        "TCGA": Dataset("TCGA_hg19", "#ff7f0e"),
        "TRACERx": Dataset("TRACERx", "#2ca02c")
    }
    return data