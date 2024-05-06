import numpy as np

def cns_head(cns, n=5):
    samples = np.sort(cns["sample_id"].unique())[:n]
    cns_head = cns.query('sample_id in @samples')
    return cns_head.copy()


def cns_tail(cns, n=5):
    samples = np.sort(cns["sample_id"].unique())[-n:]
    cns_tail = cns.query('sample_id in @samples')
    return cns_tail.copy()


def cns_random(cns, n=5, seed=0):
    np.random.seed(seed)
    samples = np.random.choice(cns["sample_id"].unique(), n, replace=False)
    cns_random = cns.query('sample_id in @samples')
    return cns_random.copy()


def sample_head(cns, n=5):
    samples = np.sort(cns["sample_id"].unique())[:n]
    sample_head = cns.query('sample_id in @samples')
    return sample_head.copy()


def sample_tail(cns, n=5):
    samples = np.sort(cns["sample_id"].unique())[-n:]
    sample_tail = cns.query('sample_id in @samples')
    return sample_tail.copy()


def sample_random(cns, n=5, seed=0):
    np.random.seed(seed)
    samples = np.random.choice(cns["sample_id"].unique(), n, replace=False)
    sample_random = cns.query('sample_id in @samples')
    return sample_random.copy()


def get_autosomes(cns):
    return cns.query("chrom != 'chrX' and chrom != 'chrY'")


def get_sex_chroms(cns):
    return cns.query("chrom == 'chrX' or chrom == 'chrY'")