import numpy as np

def cns_head(cns, n=5):
    samples = np.sort(cns.index.unique())[:n]
    cns_head = cns.query('sample_id in @samples')
    return cns_head.copy()


def cns_tail(cns, n=5):
    samples = np.sort(cns.index.unique())[-n:]
    cns_tail = cns.query('sample_id in @samples')
    return cns_tail.copy()


def cns_random(cns, n=5, seed=0):
    np.random.seed(seed)
    samples = np.random.choice(cns.index.unique(), n, replace=False)
    cns_random = cns.query('sample_id in @samples')
    return cns_random.copy()


def sample_head(cns, n=5):
    samples = np.sort(cns.index.unique())[:n]
    sample_head = cns.query('sample_id in @samples')
    return sample_head.copy()


def sample_tail(cns, n=5):
    samples = np.sort(cns.index.unique())[-n:]
    sample_tail = cns.query('sample_id in @samples')
    return sample_tail.copy()


def sample_random(cns, n=5, seed=0):
    np.random.seed(seed)
    samples = np.random.choice(cns.index.unique(), n, replace=False)
    sample_random = cns.query('sample_id in @samples')
    return sample_random.copy()


def only_aut(cns):
    return cns.query("chrom != 'chrX' and chrom != 'chrY'")


def only_sex(cns):
    return cns.query("chrom == 'chrX' or chrom == 'chrY'")


def drop_Y(cns):
    return cns.query("chrom != 'chrY'")


def filter_samples(samples, ane_min_frac = 0.001, cover_min_frac = 0.95, whitelist = False, print_info = False):
    cn_netural = samples.query(f"ane_major_cn_frac_aut < {ane_min_frac} & ane_major_cn_frac_aut < {ane_min_frac}").index
    if print_info:
        print(len(cn_netural), "samples are CN neutral")
    filtered = samples.query("(index not in @cn_netural)")

    # Find samples with low coverage (below 95% in autosomes)
    low_coverage = samples.query(f"cover_frac_aut < {cover_min_frac}").index
    if print_info:
        print(len(low_coverage), "samples have low coverage")
    filtered = samples.query("(index not in @low_coverage)")

    # Filter out CN neutral and low coverage samples 
    if whitelist:
        blacklisted = samples.query("whitelist == False").index
        if print_info:
            print(len(blacklisted), "samples are blacklisted")        
        filtered = samples.query("(index not in @blacklisted)")

    return filtered.copy()


def select_CNS_samples(cns, samples):
    return cns.query("sample_id in @samples.index")


def get_cns_for_type(cns, samples, type):
    query = f"type == '{type}'"
    ids = samples.query(query).index
    select_cns = cns.set_index("sample_id").loc[ids].reset_index()
    return select_cns