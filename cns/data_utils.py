from pathlib import Path
from cns.utils.files import load_cns, load_samples
import pandas as pd
from os.path import join as pjoin, abspath, dirname


def get_root_path():
    return abspath(pjoin(dirname(__file__), ".."))


out_path = pjoin(get_root_path(), "out")
data_path = pjoin(get_root_path(), "data")
docs_path = pjoin(get_root_path(), "docs")


def load_cns_out(filename):
    return load_cns(pjoin(out_path, filename))


def load_samples_out(filename):
    return load_samples(pjoin(out_path, filename))


def load_bins(dataset, bin_type):
    return load_cns(pjoin(out_path, f"{dataset}_bin_{bin_type}.tsv"))


class Dataset:
    def __init__(self, name, color):
        self.name = name
        self.cns = load_cns_out(f"{name}_cns_imp.tsv")
        self.samples = load_samples_out(f"{name}_samples.tsv")
        self.color = color


def load_pcawg():
    return Dataset("PCAWG", "#1f77b4")


def load_tcga():
    return Dataset("TCGA_hg19", "#ff7f0e")


def load_tracerx():
    return Dataset("TRACERx", "#2ca02c")


def load_data():
    return {
        "PCAWG" : load_pcawg(),
        "TCGA": load_tcga(),
        "TRACERx": load_tracerx()
    }


def load_data_file(filename):
    sep = "\t" if filename.endswith(".tsv") else ","
    return pd.read_csv(pjoin(data_path, filename), sep=sep)