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


def only_aut(cns_df, assembly=hg19, inplace=False):
    if inplace:
        if not cns_df.index.is_unique:
            cns_df.reset_index(drop=True, inplace=True)
        cns_df.drop(cns_df[(cns_df['chrom'] == assembly.chr_y) | (cns_df['chrom'] == assembly.chr_x)].index, inplace=True)
        return cns_df
    return cns_df.query(f"chrom != '{assembly.chr_x}' and chrom != '{assembly.chr_y}'")


def only_sex(cns_df, assembly=hg19, inplace=False):
    if inplace:
        if not cns_df.index.is_unique:
            cns_df.reset_index(drop=True, inplace=True)
        cns_df.drop(cns_df[(cns_df['chrom'] != assembly.chr_y) & (cns_df['chrom'] != assembly.chr_x)].index, inplace=True)
        return cns_df
    return cns_df.query(f"chrom == '{assembly.chr_x}' or chrom == '{assembly.chr_y}'")


def drop_Y(cns_df, assembly=hg19, inplace=False):
    if inplace:
        cns_df.drop(cns_df[cns_df['chrom'] == assembly.chr_y].index, inplace=True)
        return cns_df
    return cns_df.query(f"chrom != '{assembly.chr_y}'").copy()


def select_CNS_samples(cns_df, samples_df):
    return cns_df.query("sample_id in @samples_df.index")


def select_cns_by_type(cns_df, samples, type_val, type_col="type"):
    query = f"{type_col} == '{type_val}'"
    ids = samples.query(query).index
    cns_ids = cns_df["sample_id"].unique()
    intersect = np.intersect1d(ids, cns_ids)
    select_cns = cns_df.set_index("sample_id").loc[intersect].reset_index()
    return select_cns


def cn_not_nan(cns_df, cn_columns, het):
    nan_vals = cns_df[cn_columns].isna()
    nan_filter = ~nan_vals.all(axis=1) if het else ~nan_vals.any(axis=1)
    non_nan_df = cns_df.loc[nan_filter]
    return non_nan_df


def get_chr_sets(cns_df, assembly=hg19):
    chroms = cns_df["chrom"].unique().tolist()
    aut_selection = [chrom for chrom in chroms if chrom in assembly.aut_names]
    if len(aut_selection) == 0:
        print(cns_df)
        raise ValueError("No autosomes found in the input segments.")
    res_dict = { "aut": aut_selection}
    sex_selection = [chrom for chrom in chroms if chrom in assembly.sex_names]
    if len(sex_selection) != 0:
        res_dict["sex"] = sex_selection
        res_dict["all"] = aut_selection + sex_selection
    return res_dict


def dataframe_array_split(samples_df, n_splits):
    """
    Splits a DataFrame into n_splits parts as equally as possible.

    Parameters:
    - samples_df: The pandas DataFrame to split.
    - n_splits: The number of parts to split the DataFrame into.

    Returns:
    - A list of pandas DataFrame objects.
    """
    # Ensure n_splits is a positive integer
    n_splits = max(int(n_splits), 1)

    if n_splits == 1:
        return [samples_df]

    # Calculate the number of rows in each split
    total_rows = len(samples_df)
    rows_per_split = total_rows // n_splits
    remainder = total_rows % n_splits

    # Initialize variables to keep track of the current row and the list of splits
    current_row = 0
    splits = []

    for i in range(n_splits):
        # Calculate the number of rows for this split
        # Add one to some of the splits to distribute the remainder
        rows_in_split = rows_per_split + (1 if i < remainder else 0)

        # Slice the DataFrame for this split and append to the list
        split_df = samples_df.iloc[current_row:current_row + rows_in_split]
        if rows_in_split > 0:
            splits.append(split_df)
            # Update the current row index
            current_row += rows_in_split

    return splits