import pandas as pd
from os.path import join as pjoin, abspath, dirname

from cns.process.binning import add_cns_loc, sum_cns
from cns.utils.selection import select_CNS_samples
from cns.utils.files import load_cns, load_samples


def get_root_path():
    return abspath(pjoin(dirname(__file__), ".."))


out_path = pjoin(get_root_path(), "out")
data_path = pjoin(get_root_path(), "data")
docs_path = pjoin(get_root_path(), "docs")


def load_cns_out(filename):
    return sum_cns(add_cns_loc(load_cns(pjoin(out_path, filename))))


def load_samples_out(filename):
    return load_samples(pjoin(out_path, filename))


def load_bins(dataset, bin_type):
    return sum_cns(add_cns_loc(load_cns(pjoin(out_path, f"{dataset}_bin_{bin_type}.tsv"))))


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
    data = {
        "PCAWG" : load_pcawg(),
        "TCGA": load_tcga(),
        "TRACERx": load_tracerx()
    }
    return data


def load_data_file(filename):
    sep = "\t" if filename.endswith(".tsv") else ","
    return pd.read_csv(pjoin(data_path, filename), sep=sep)


def filter_samples(samples, ane_min_frac = 0.001, cover_min_frac = 0.95, whitelist = False, remove_uncertain = False, print_info = False):
    cn_neutral = samples.query(f"ane_major_cn_frac_aut < {ane_min_frac} & ane_major_cn_frac_aut < {ane_min_frac}").index
    if print_info:
        print(len(cn_neutral), "samples are CN neutral")
    filtered = samples.query("(index not in @cn_neutral)")

    # Find samples with low coverage (below 95% in autosomes)
    low_coverage = samples.query(f"cover_frac_aut < {cover_min_frac}").index
    if print_info:
        print(len(low_coverage), "samples have low coverage")
    filtered = filtered.query("(index not in @low_coverage)")

    # Filter out CN neutral and low coverage samples 
    if whitelist:
        blacklisted = samples.query("whitelist == False").index
        if print_info:
            print(len(blacklisted), "samples are blacklisted")
        filtered = filtered.query("(index not in @blacklisted)")

    if remove_uncertain:
        samples["type"] = samples["type"].replace({"LUADx2": "LUAD"})
        samples["type"] = samples["type"].replace({"LUADx3": "LUAD"})
        untyped = samples[samples["type"].fillna('').apply(lambda x: any(not c.isupper() for c in x))].index
        if print_info:
            print(len(untyped), "samples do not have exact type")
        filtered = filtered.query("(index not in @untyped)")

    if print_info:
        print("Filtered samples:", len(filtered))

    return filtered.copy()


def load_filter_samples(print_info=False):
    samples = {
        "PCAWG": load_samples_out("PCAWG_samples.tsv"),
        "TRACERx": load_samples_out("TRACERx_samples.tsv"),
        "TCGA": load_samples_out("TCGA_hg19_samples.tsv")
    }
    for k, v in samples.items():
        if print_info:
            print(k)
        min_frac = 0.95 if k != "TRACERx" else 0.85
        samples[k] = filter_samples(v, cover_min_frac=min_frac, whitelist=k=="PCAWG", remove_uncertain=k=="TRACERx", print_info=print_info)
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

# TODO: Double check if works
def get_cns_for_type(cns, samples, type):
	query = f"type == '{type}' | TCGA_type == '{type}'" if "TCGA_type" in samples.columns else f"type == '{type}'"
	ids = samples.query(query).index
	select_cns = cns.set_index("sample_id").loc[ids].reset_index()
	return select_cns


def load_merged_samples(print_info=False):
    samples = load_filter_samples(print_info)
    for k, v in samples.items():
        v["source"] = k
    samples["PCAWG"]["type"] = samples["PCAWG"]["TCGA_type"]    
    samples["TRACERx"]["type"] = samples["TRACERx"]["type"].replace({"LUADx2": "LUAD"})
    samples["TRACERx"]["type"] = samples["TRACERx"]["type"].replace({"LUADx3": "LUAD"})
    # drop where TCGA_id is != NaN
    overlap_with_tcga = samples["PCAWG"].index[samples["PCAWG"]["TCGA_id"].notna()]
    if print_info:
        print(f"Overlapping samples with TCGA: {len(overlap_with_tcga)}")
    samples["PCAWG"] = samples["PCAWG"].drop(overlap_with_tcga)
    # drop columns TCGA_id and TCGA_type
    samples["PCAWG"] = samples["PCAWG"].drop(columns=["TCGA_id", "TCGA_type", "whitelist"])
    samples["PCAWG"].head()
    all_samp = pd.concat(samples.values())
    if print_info:
        print("Total samples:", len(all_samp))
    return all_samp


def load_merged_bins(samples, bin_size):
    cns = {
        "PCAWG": load_cns_out(f"PCAWG_bin_{bin_size}.tsv"),
        "TRACERx": load_cns_out(f"TRACERx_bin_{bin_size}.tsv"),
        "TCGA": load_cns_out(f"TCGA_hg19_bin_{bin_size}.tsv")
    }
    all_cns = pd.concat(cns.values())
    all_cns = select_CNS_samples(all_cns, samples)
    return all_cns


def load_merged_cns(samples):
    cns = {
        "PCAWG": load_cns_out("PCAWG_cns_imp.tsv"),
        "TRACERx": load_cns_out("TRACERx_cns_imp.tsv"),
        "TCGA": load_cns_out("TCGA_hg19_cns_imp.tsv")
    }
    all_cns = pd.concat(cns.values())
    all_cns = select_CNS_samples(all_cns, samples)
    return all_cns


def main_load_data():
    samples = load_merged_samples()
    cns = load_merged_cns(samples)
    return samples, cns


def load_COSMIC(change_coords=True):
    res = pd.read_csv(pjoin(data_path, "COSMIC_consensus_genes.tsv"), sep="\t")
    if change_coords:
        res.loc[:, "start"] -= 1
    return res


def load_ENSEMBL(change_coords=True):
    res = pd.read_csv(pjoin(data_path, "ENSEMBL_coding_genes.tsv"), sep="\t")
    if change_coords:
        res.loc[:, "start"] -= 1
    return res