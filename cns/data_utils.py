import pandas as pd
from os.path import join as pjoin, abspath, dirname
from cns.process.binning import add_cns_loc, sum_cns
from cns.utils.canonization import find_cn_cols_if_none
from cns.utils.logging import log_info
from cns.utils.selection import select_CNS_samples
from cns.utils.files import load_cns, load_samples
from cns.utils.kneepoint import find_bends
import matplotlib.pyplot as plt


def get_root_path():
    return abspath(pjoin(dirname(__file__), ".."))


img_path = pjoin(get_root_path(), "img")
out_path = pjoin(get_root_path(), "out")
data_path = pjoin(get_root_path(), "data")
docs_path = pjoin(get_root_path(), "docs")


def load_cns_out(filename, raw=False):
    cns_df = load_cns(pjoin(out_path, filename))
    if raw:
        return cns_df
    return sum_cns(add_cns_loc(rename_cns_columns(cns_df)))


def load_samples_out(filename):
    return load_samples(pjoin(out_path, filename))


def load_bins(dataset, bin_type):
    return sum_cns(add_cns_loc(load_cns(pjoin(out_path, f"{dataset}_bin_{bin_type}.tsv"))))


def filter_samples(samples, ane_min_frac=0.001, cover_min_frac=0.95, whitelist=False, filter_types=False, print_info=False):
    log_info(print_info, f"Total samples: {len(samples)}")
    
    cn_neutral = samples.query(f"ane_het_aut < {ane_min_frac}").index
    log_info(print_info, f"{len(cn_neutral)} samples are CN neutral (below {ane_min_frac:.5f})")
    filtered = samples.query("(index not in @cn_neutral)")

    # Find samples with low coverage (below 95% in autosomes)
    low_coverage = samples.query(f"cover_het_aut < {cover_min_frac}").index
    log_info(print_info, f"{len(low_coverage)} samples have low coverage (below {cover_min_frac:.5f})")
    filtered = filtered.query("(index not in @low_coverage)")

    # Filter out CN neutral and low coverage samples 
    if whitelist:
        blacklisted = samples.query("whitelist == False").index
        log_info(print_info, f"{len(blacklisted)} samples are blacklisted")
        filtered = filtered.query("(index not in @blacklisted)")

    if filter_types:
        samples["type"] = samples["type"].replace({"LUADx2": "LUAD"}).replace({"LUADx3": "LUAD"})
        untyped = samples[samples["type"].fillna('').apply(lambda x: any(not c.isupper() for c in x))].index
        log_info(print_info, f"{len(untyped)} samples do not have exact type")
        filtered = filtered.query("(index not in @untyped)")

    log_info(print_info, f"Filtered samples: {len(filtered)}")

    return filtered.copy()


def load_all_samples(filter=True, retype=True, drop_tcga=True, print_info=False):
    samples = {
        "PCAWG": load_samples_out("PCAWG_samples.tsv"),
        "TRACERx": load_samples_out("TRACERx_samples.tsv"),
        "TCGA_hg19": load_samples_out("TCGA_hg19_samples.tsv")
    }
    
    overlap_with_tcga = samples["PCAWG"]["TCGA_id"].dropna().unique()
    
    if filter:
        for k, v in samples.items():
            log_info(print_info, k)
            ane_vals = v["ane_het_aut"]
            ane_bends = find_bends(ane_vals)
            ane_min_frac = ane_bends[0][ane_bends[2]]
            cover_vals = v["cover_het_aut"]
            cover_bends = find_bends(cover_vals)
            cover_min_frac = cover_bends[0][cover_bends[4]]
            whitelist = k=="PCAWG"
            filter_types = k=="TRACERx"
            samples[k] = filter_samples(v, ane_min_frac, cover_min_frac, whitelist, filter_types, print_info)
    
    if drop_tcga:
        # drop where TCGA_id is != NaN
        samples["TCGA_hg19"] = samples["TCGA_hg19"].query("sample_id not in @overlap_with_tcga") 
        log_info(print_info, f"Overlapping samples with PCAWG: {len(overlap_with_tcga)}")
        log_info(print_info, f"After overlap removal: {len(samples['TCGA_hg19'])}")

    if retype:
        samples["PCAWG"]["type"] = samples["PCAWG"]["TCGA_type"]    
        samples["TRACERx"]["type"] = samples["TRACERx"]["type"].replace({"LUADx2": "LUAD"}).replace({"LUADx3": "LUAD"})
    
    samples["PCAWG"] = samples["PCAWG"].drop(columns=["TCGA_id", "TCGA_type", "whitelist"])

    return samples


def get_cns_for_type(cns, samples, type):
    query = f"type == '{type}'"
    ids = samples.query(query).index
    select_cns = cns.set_index("sample_id").loc[ids].reset_index()
    return select_cns


def load_merged_samples(filter=True, retype=True, print_info=False):
    samples = load_all_samples(filter, retype, drop_tcga=True, print_info=print_info)
    for k, v in samples.items():
        v["source"] = k
    all_samples = pd.concat(samples.values())        
    log_info(print_info, f"Total samples: {len(all_samples)}")
    return all_samples


def rename_cns_columns(cns_df, cn_columns=None):
    cn_columns = find_cn_cols_if_none(cns_df, cn_columns)
    return cns_df.rename(columns={cn_columns[0]: "major_cn", cn_columns[1]: "minor_cn"})


def load_merged_bins(select_samples, bin_size):
    cns = {
        "PCAWG": load_cns_out(f"PCAWG_bin_{bin_size}.tsv"),
        "TRACERx": rename_cns_columns(load_cns_out(f"TRACERx_bin_{bin_size}.tsv")),
        "TCGA_hg19": load_cns_out(f"TCGA_hg19_bin_{bin_size}.tsv")
    }
    all_cns = pd.concat(cns.values())
    if select_samples is not None:
        all_cns = select_CNS_samples(all_cns, select_samples)
    return all_cns


def load_merged_cns(select_samples=None):
    cns = {
        "PCAWG": load_cns_out("PCAWG_cns_imp.tsv"),
        "TRACERx": rename_cns_columns(load_cns_out("TRACERx_cns_imp.tsv")),
        "TCGA_hg19": load_cns_out("TCGA_hg19_cns_imp.tsv")
    }
    all_cns = pd.concat(cns.values())
    cns = None
    if select_samples is not None:
        all_cns = select_CNS_samples(all_cns, select_samples)
    return all_cns


def main_load_data(bins = None):
    samples = load_merged_samples()
    if bins == None:
        cns = load_merged_cns(samples)
    else:
        cns = load_merged_bins(samples, bins)
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


def save_cns_fig(fig_name, fig = None):
    if fig == None:
        fig = plt.gcf()
    fig.savefig(f"{img_path}/{fig_name}.png", bbox_inches="tight", transparent=True, dpi=300)
    fig.savefig(f"{img_path}/{fig_name}.pdf", bbox_inches="tight", transparent=True)


def save_doc_fig(fig_name, fig = None):
    if fig == None:
        fig = plt.gcf()
    fig.savefig(f"{docs_path}/{fig_name}.png", bbox_inches="tight", transparent=True, dpi=300)
    fig.savefig(f"{img_path}/{fig_name}.pdf", bbox_inches="tight", transparent=True)