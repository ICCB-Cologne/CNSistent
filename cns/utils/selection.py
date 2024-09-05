import numpy as np
from cns.utils.assemblies import hg19

def cns_head(cns_df, n=5):
    samples = np.sort(cns_df["sample_id"].unique())[:n]
    cns_head = cns_df.query('sample_id in @samples')
    return cns_head.copy()


def cns_tail(cns_df, n=5):
    samples = np.sort(cns_df["sample_id"].unique())[-n:]
    cns_tail = cns_df.query('sample_id in @samples')
    return cns_tail.copy()


def cns_random(cns_df, n=5, seed=0):
    np.random.seed(seed)
    samples = np.random.choice(cns_df["sample_id"].unique(), n, replace=False)
    cns_random = cns_df.query('sample_id in @samples')
    return cns_random.copy()


def sample_head(samples_df, n=5):
    samples = np.sort(samples_df.index.unique())[:n]
    sample_head = samples_df.query('sample_id in @samples')
    return sample_head.copy()


def sample_tail(samples_df, n=5):
    samples = np.sort(samples_df.index.unique())[:n]
    sample_tail = samples.query('sample_id in @samples')
    return sample_tail.copy()


def sample_random(samples_df, n=5, seed=0):
    np.random.seed(seed)
    samples = np.random.choice(samples_df["sample_id"].unique(), n, replace=False)
    sample_random = samples_df.query('sample_id in @samples')
    return sample_random.copy()


def only_aut(cns_df, assembly=hg19):
    return cns_df.query(f"chrom != '{assembly.chr_x}' and chrom != '{assembly.chr_y}'").copy()


def only_sex(cns_df, assembly=hg19):
    return cns_df.query(f"chrom == '{assembly.chr_x}' or chrom == '{assembly.chr_y}'").copy()


def drop_Y(cns_df, assembly=hg19):
    return cns_df.query(f"chrom != '{assembly.chr_y}'").copy()


def select_CNS_samples(cns_df, samples):
    return cns_df.query("sample_id in @samples.index")


def get_cns_for_type(cns_df, samples, type):
    query = f"type == '{type}'"
    ids = samples.query(query).index
    select_cns = cns_df.set_index("sample_id").loc[ids].reset_index()
    return select_cns