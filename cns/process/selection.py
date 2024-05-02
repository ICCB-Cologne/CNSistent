import numpy as np

def cna_head(cns, n=5):
    samples = np.sort(cns["sample_id"].unique())[:n]
    cna_head = cns.query('sample_id in @samples')
    return cna_head.copy()


def cna_tail(cns, n=5):
    samples = np.sort(cns["sample_id"].unique())[-n:]
    cna_tail = cns.query('sample_id in @samples')
    return cna_tail.copy()


def cna_random(cns, n=5, seed=0):
    np.random.seed(seed)
    samples = np.random.choice(cns["sample_id"].unique(), n, replace=False)
    cna_random = cns.query('sample_id in @samples')
    return cna_random.copy()


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