import pandas as pd
from os.path import join as pjoin, abspath, dirname
from cns.utils.canonization import find_cn_cols_if_none
from cns.utils.logging import log_info
from cns.utils.selection import select_CNS_samples
from cns.utils.files import load_cns, load_samples, load_segments
from cns.utils.kneepoint import find_bends, z_score_filter
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
    return rename_cns_columns(cns_df)


def load_samples_out(filename):
    return load_samples(pjoin(out_path, filename))


def load_bins(dataset, segment_type):
    return load_cns(pjoin(out_path, f"{dataset}_bin_{segment_type}.tsv"))


def filter_samples(samples, ane_min_frac=0.001, cover_min_frac=0.95, filter_types=False, print_info=False):
    log_info(print_info, f"Total samples: {len(samples)}")
    
    cn_neutral = samples.query(f"ane_het_aut < {ane_min_frac}").index
    log_info(print_info, f"{len(cn_neutral)} samples are CN neutral (below {ane_min_frac:.5f})")
    filtered = samples.query("(index not in @cn_neutral)")

    # Find samples with low coverage (below 95% in autosomes)
    low_coverage = samples.query(f"cover_het_aut < {cover_min_frac}").index
    log_info(print_info, f"{len(low_coverage)} samples have low coverage (below {cover_min_frac:.5f})")
    filtered = filtered.query("(index not in @low_coverage)")

    if filter_types:
        samples["type"] = samples["type"].replace({"LUADx2": "LUAD"}).replace({"LUADx3": "LUAD"})
        untyped = samples[samples["type"].fillna('').apply(lambda x: any(not c.isupper() for c in x))].index
        log_info(print_info, f"{len(untyped)} samples do not have exact type")
        filtered = filtered.query("(index not in @untyped)")

    log_info(print_info, f"Filtered samples: {len(filtered)}")

    return filtered.copy()


def load_all_samples(filter=True, retype=True, print_info=False):
    samples = {
        "PCAWG": load_samples_out("PCAWG_samples.tsv"),
        "TRACERx": load_samples_out("TRACERx_samples.tsv"),
        "TCGA_hg19": load_samples_out("TCGA_hg19_samples.tsv")
    }
    total_count = sum([len(v) for v in samples.values()])
    log_info(print_info, f"Total samples: {total_count}")

    if filter:
        for k, v in samples.items():
            log_info(print_info, k)

            # calculate bend for aneuploidy
            ane_bends = find_bends(v["ane_het_aut"])
            ane_min_frac = ane_bends[0][ane_bends[2]]

            # calculate the z-score for coverage
            cover_filtered = z_score_filter(v["cover_het_aut"])
            cover_min_frac = cover_filtered.min()

            filter_types = k=="TRACERx" and retype
            samples[k] = filter_samples(v, ane_min_frac, cover_min_frac, filter_types, print_info)

    if retype:
        samples["PCAWG"]["type"] = samples["PCAWG"]["TCGA_type"]    
        samples["PCAWG"] = samples["PCAWG"].drop(columns=["TCGA_type"])
        samples["TRACERx"]["type"] = samples["TRACERx"]["type"].replace({"LUADx2": "LUAD"}).replace({"LUADx3": "LUAD"})

    return samples


def get_cns_for_type(cns, samples, type):
    query = f"type == '{type}'"
    ids = samples.query(query).index
    select_cns = cns.set_index("sample_id").loc[ids].reset_index()
    return select_cns


def load_merged_samples(filter=True, retype=True, print_info=False):
    samples = load_all_samples(filter, retype, print_info=print_info)
    for k, v in samples.items():
        v["source"] = k
    all_samples = pd.concat(samples.values())        
    log_info(print_info, f"Total samples: {len(all_samples)}")
    return all_samples


def rename_cns_columns(cns_df, cn_columns=None):
    cn_columns = find_cn_cols_if_none(cns_df, cn_columns)
    return cns_df.rename(columns={cn_columns[0]: "major_cn", cn_columns[1]: "minor_cn"})


def load_merged_bins(select_samples, segment_size):
    cns = {
        "PCAWG": load_cns_out(f"PCAWG_bin_{segment_size}.tsv"),
        "TRACERx": load_cns_out(f"TRACERx_bin_{segment_size}.tsv"),
        "TCGA_hg19": load_cns_out(f"TCGA_hg19_bin_{segment_size}.tsv")
    }
    all_cns = pd.concat(cns.values())
    if select_samples is not None:
        all_cns = select_CNS_samples(all_cns, select_samples)
    return all_cns.reset_index(drop=True)


def load_merged_cns(select_samples=None):
    cns = {
        "PCAWG": load_cns_out("PCAWG_cns_imp.tsv"),
        "TRACERx": load_cns_out("TRACERx_cns_imp.tsv"),
        "TCGA_hg19": load_cns_out("TCGA_hg19_cns_imp.tsv")
    }
    all_cns = pd.concat(cns.values())
    cns = None
    if select_samples is not None:
        all_cns = select_CNS_samples(all_cns, select_samples)
    return all_cns.reset_index(drop=True)


def main_load_data(segment_type = None):
    samples = load_merged_samples()
    if segment_type == None:
        cns = load_merged_cns(samples)
    else:
        cns = load_merged_bins(samples, segment_type)
    return samples, cns


def load_COSMIC():
    return load_segments(pjoin(data_path, "COSMIC_consensus_genes.bed"))


def load_ENSEMBL():
    return load_segments(pjoin(data_path, "ENSEMBL_coding_genes.bed"))


def save_cns_fig(fig_name, fig = None):
    if fig == None:
        fig = plt.gcf()
    fig.savefig(f"{img_path}/{fig_name}.png", bbox_inches="tight", transparent=True, dpi=300)
    fig.savefig(f"{img_path}/{fig_name}.pdf", bbox_inches="tight", transparent=True)


def save_doc_fig(fig_name, fig = None):
    if fig == None:
        fig = plt.gcf()
    fig.savefig(f"{docs_path}/{fig_name}.png", bbox_inches="tight", transparent=True, dpi=300)
    fig.savefig(f"{docs_path}/{fig_name}.pdf", bbox_inches="tight", transparent=True)