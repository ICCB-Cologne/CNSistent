from cns.utils.files import load_cns, load_samples
import pandas as pd
from os.path import join as pjoin


out_path = "../../out"
data_path = "../../data"
docs_path = "../../docs"


def load_cns_out(filename):
    return load_cns(pjoin(out_path, filename))


def load_samples_out(filename):
    return load_samples(pjoin(out_path, filename))


class Dataset:
    def __init__(self, name, color):
        self.name = name
        self.cns = load_cns_out(f"{name}_cns_imp.tsv")
        self.samples = load_samples_out(f"{name}_samples.tsv")
        self.color = color


def load_data():
    data = {
        "PCAWG" : Dataset("PCAWG", "#1f77b4"),
        "TCGA": Dataset("TCGA_hg19", "#ff7f0e"),
        "TRACERx": Dataset("TRACERx", "#2ca02c")
    }
    return data