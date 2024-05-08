from pathlib import Path
from cns.utils.files import load_cns, load_samples
import pandas as pd
from os.path import join as pjoin, abspath, dirname

from cns.utils.selection import filter_samples, select_CNS_samples


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


def load_filter_samples():
    samples = {
        "PCAWG": load_samples_out("PCAWG_samples.tsv"),
        "TRACERx": load_samples_out("TRACERx_samples.tsv"),
        "TCGA": load_samples_out("TCGA_hg19_samples.tsv")
    }
    for k, v in samples.items():
        print(k)
        min_frac = 0.95 if k != "TRACERx" else 0.85
        samples[k] = filter_samples(v, cover_min_frac=min_frac, whitelist=k=="PCAWG")
    return samples


def load_filter_bins(samples, bin_size):
    cns = {
        "PCAWG": load_cns_out(f"PCAWG_bin_{bin_size}.tsv"),
        "TRACERx": load_cns_out(f"TRACERx_bin_{bin_size}.tsv"),
        "TCGA": load_cns_out(f"TCGA_hg19_bin_{bin_size}.tsv")
    }
    for k, v in cns.items():
        cns[k] = select_CNS_samples(v, samples[k])
    return cns


def get_cns_for_type(cns, samples, type):
	query = f"type == '{type}' | TCGA_type == '{type}'" if "TCGA_type" in samples.columns else f"type == '{type}'"
	ids = samples.query(query).index
	select_cns = cns.set_index("sample_id").loc[ids].reset_index()
	return select_cns