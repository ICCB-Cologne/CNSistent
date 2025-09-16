import pandas as pd
import re
from cns.utils.assemblies import hg19
from cns.utils.logging import log_info

def _requires_rename(cn_columns):
    if len(cn_columns) > 2:
        raise ValueError(f"""Discovery of CN columns failed.\n
                         Only one (total) or two (major, minor) or (hap1, hap2) CN columns are allowed.\n
                         Found {cn_columns} instead.""")
    elif len(cn_columns) == 2:
        for cn_col in cn_columns:
            if not cn_col in ["major_cn", "minor_cn", "hap1_cn", "hap2_cn"]:
                return True

    elif len(cn_columns) == 1:
        return cn_columns[0] != "total_cn"

    else:
        raise ValueError("Discovery of CN columns failed. No CN columns found.")

    return False


# Either obtains CN columns from the DataFrame or checks if the provided columns are valid
def get_cn_cols(cns_df, cn_cols=None):
    """Gets or validates copy number columns from DataFrame.
    
    Args:
        cns_df (pd.DataFrame): DataFrame containing copy number data
        cn_cols (str|list, optional): Column name(s) to validate. If None, discovers CN columns
    
    Returns:
        list: Valid copy number column names
        
    Raises:
        ValueError: If specified columns not found or invalid number of columns
    """
    if cn_cols is None:
        cn_cols = [col for col in cns_df.columns if _is_cn_column(col)]
    # check if cn_cols is a string
    elif isinstance(cn_cols, str):
        if cn_cols in cns_df.columns:
            return [cn_cols]
        else:
            raise ValueError(f"Column {cn_cols} not found in the CNS DataFrame.")
    # check if cn_cols is a list
    elif isinstance(cn_cols, list):
        for col in cn_cols:
            if col not in cns_df.columns:
                raise ValueError(f"Column {col} not found in the CNS DataFrame.")
        if len(cn_cols) > 2:
            raise ValueError(f"Only one (total) or two (major, minor) CN columns are allowed. Found {cn_cols} instead.")
        if len(cn_cols) == 0:
            raise ValueError("No CN columns found.")
        return cn_cols
    return cn_cols


def rename_cn_cols(cns_df, cn_columns=None, print_info=False):
    """
    Renames copy number columns in the DataFrame.

    Parameters:
    cns_df (pd.DataFrame): The DataFrame containing copy number data.
    cn_columns (list): List of column names to be renamed.
    print_info (bool): Flag to indicate whether to print information.

    Returns:
    pd.DataFrame: The DataFrame with renamed columns.
    list: The updated list of column names.
    """
    cn_columns = get_cn_cols(cns_df, cn_columns)
    if not _requires_rename(cn_columns):
        return cns_df, list(cn_columns)

    if len(cn_columns) == 2:
        rename_map = _get_major_minor_cols(cns_df, cn_columns)
        if len(rename_map) == 0: # Set haplotype specific
            rename_map={cn_columns[0]: "hap1_cn", cn_columns[1]: "hap2_cn"}
    elif len(cn_columns) == 1:
        rename_map = {cn_columns[0]: "total_cn"}

    cns_df.rename(columns=rename_map, inplace=True)
    log_info(print_info, f"Renamed CN columns: {rename_map}")
    return cns_df, list(rename_map.values())


def _get_major_minor_cols(cns_df, cn_columns):
    col1 = cn_columns[0]
    col2 = cn_columns[1]
    if (cns_df[col1] >= cns_df[col2]).all():
        return {col1: "major_cn", col2: "minor_cn"}
    elif (cns_df[col2] >= cns_df[col1]).all():
        return {col1: "minor_cn", col2: "major_cn"}
    else:
        return {}  


def _find_column(cns_df, patterns):
    # Find matching column
    matching_column = None
    for col in cns_df.columns:
        if any(pd.Series(col).str.contains(pattern, case=False, regex=True).any() for pattern in patterns):
            matching_column = col
            break
    return matching_column


def canonize_sample_id(df, print_info=False):
    """Ensures DataFrame has standardized sample_id column.
    
    Args:
        df (pd.DataFrame): DataFrame to canonize
        print_info (bool): Whether to print info messages
        
    Returns:
        pd.DataFrame: DataFrame with canonized sample_id column
    """
    # if the column sample_id does not exist, rename the first column to sample_id
    if "sample_id" not in df.columns:
        sample_col = _find_column(df, ['sample', 'id', 'sampleId', 'sample_id', 'sample-id', 'sample_name', 'sampleName', 'sample-name'])	
        if sample_col is None:
            df.columns = ["sample_id"] + df.columns[1:].tolist()
            log_info(print_info, f"Column sample_id not found, renamed first column to sample_id.")
        else:
            df.rename(columns={sample_col: "sample_id"}, inplace=True)
            log_info(print_info, f"Renamed column {sample_col} to sample_id.")

    return df


def canonize_chroms(cns_df, assembly=hg19, print_info=False):
    """Standardizes chromosome names and validates against assembly.
    
    Args:
        cns_df (pd.DataFrame): Copy number DataFrame
        assembly: Reference assembly defining valid chromosomes
        print_info (bool): Whether to print info messages
        
    Returns:
        pd.DataFrame: DataFrame with canonized chromosome names
        
    Raises:
        ValueError: If no valid chromosomes found
    """
        # if the column chrom does not exist, rename the second column to chrom
    if "chrom" not in cns_df.columns:
        chrom_col = _find_column(cns_df, ['chrom', 'chr', 'chromosome'])
        if chrom_col is None:
            cns_df.columns = cns_df.columns[:1].tolist() + ["chrom"] + cns_df.columns[2:].tolist()
            log_info(print_info, f"Column chrom not found, renamed second column to chrom.")
        else:
            cns_df.rename(columns={chrom_col: "chrom"}, inplace=True)
            log_info(print_info, f"Renamed column {chrom_col} to chrom.")

    chrom_vals = cns_df["chrom"].unique()
    # if the chromosomes values are all either digits or single characters, convert to chrX format
    if all([str(chrom).isdigit() or len(chrom) == 1 for chrom in chrom_vals]):
        cns_df["chrom"] = "chr" + cns_df["chrom"].astype(str)
        chrom_vals = cns_df["chrom"].unique()
        log_info(print_info, "Chromosome values converted to chr[1-Y] format.")
    # if the first 3 letters of the chromosome values are not lower case, convert these 3 letters to lower case
    if not all([chrom[:3].islower() for chrom in chrom_vals]):
        cns_df["chrom"] = cns_df["chrom"].apply(lambda x: x[:3].lower() + x[3:])
        chrom_vals = cns_df["chrom"].unique()
        log_info(print_info, "Chromosome values converted to lower case.")
    
    if not any([chrom in assembly.chr_names for chrom in chrom_vals]):
        raise ValueError(f"No chrom found. Chromosome values must be in {assembly.chr_names}, got {chrom_vals}.")
    not_known = [chrom for chrom in chrom_vals if chrom not in assembly.chr_names]
    if len(not_known) > 0:
        log_info(print_info, f"Found chromosomes not in assembly: {not_known}, these will be dropped.")
        rows_to_drop = cns_df[cns_df["chrom"].isin(not_known)].index
        cns_df.drop(rows_to_drop, inplace=True) 
    return cns_df


def canonize_positions(cns_df, print_info=False):
    """Ensures standard start/end position columns.
    
    Args:
        cns_df (pd.DataFrame): Copy number DataFrame
        print_info (bool): Whether to print info messages
        
    Returns:
        pd.DataFrame: DataFrame with canonized position columns
    """

        # if the column start does not exist, rename the third column to start
    if "start" not in cns_df.columns:
        start_col = _find_column(cns_df, ['start', 'begin', 'chromstart', 'chrom-start', 'chrom_start', 'startpos', 'start-pos', 'start_pos'])
        if start_col is None:
            cns_df.columns = cns_df.columns[:2].tolist() + ["start"] + cns_df.columns[3:].tolist()
            log_info(print_info, f"Column start not found, renamed third column to start.")
        else:
            cns_df.rename(columns={start_col: "start"}, inplace=True)
            log_info(print_info, f"Renamed column {start_col} to start.")
    cns_df["start"] = cns_df["start"].astype(int)

    # if the column end does not exist, rename the fourth column to end
    if "end" not in cns_df.columns:
        end_col = _find_column(cns_df, ['end', 'stop', 'endpos', 'end-pos', 'end_pos' 'chromend', 'chrom-end', 'chrom_end'])
        if end_col is None:
            cns_df.columns = cns_df.columns[:3].tolist() + ["end"] + cns_df.columns[4:].tolist()
            log_info(print_info, f"Column end not found, renamed fourth column to end.")
        else:
            cns_df.rename(columns={end_col: "end"}, inplace=True)
            log_info(print_info, f"Renamed column {end_col} to end.")
    return cns_df


def canonize_name(cns_df, print_info=False):
    """Standardizes optional name column if present.
    
    Args:
        cns_df (pd.DataFrame): Copy number DataFrame
        print_info (bool): Whether to print info messages
        
    Returns:
        pd.DataFrame: DataFrame with canonized name column
    """
    name_col = _find_column(cns_df, ['name'])
    if name_col is not None and name_col != "name":
        cns_df.rename(columns={name_col: "name"}, inplace=True)
        log_info(print_info, f"Renamed column {name_col} to name.")
    return cns_df

def canonize_cns_df(cns_df, input_cn_columns=None, order_columns=False, assembly=hg19, print_info=False):
    """Canonizes copy number DataFrame by standardizing all columns.
    
    Applies standard canonization steps in order:
    1. Sample ID column
    2. Chromosome names and validation
    3. Position columns (start/end)
    4. Copy number columns
    5. Optional name column
    
    Args:
        cns_df (pd.DataFrame): Copy number DataFrame to canonize
        assembly: Reference assembly for chromosome validation
        input_cn_columns (str|list, optional): Copy number column(s) to validate
        print_info (bool): Whether to print info messages
        
    Returns:
        pd.DataFrame: Canonized copy number DataFrame
        list: Names of canonized copy number columns
        
    Raises:
        ValueError: If required columns missing or invalid
    """
    # convert columns to strings
    cns_df.columns = cns_df.columns.astype(str)
    cn_columns = get_cn_cols(cns_df, input_cn_columns)

    cns_df = canonize_sample_id(cns_df, print_info)
    cns_df = canonize_chroms(cns_df, assembly, print_info)
    cns_df = canonize_positions(cns_df, print_info)
    cns_df = canonize_name(cns_df, print_info)

    log_info(print_info, f"Using CN columns: {cn_columns}")
    if len(cn_columns) == 2 and order_columns:
        major_cn = cns_df[[cn_columns[0], cn_columns[1]]].max(axis=1)
        minor_cn = cns_df[[cn_columns[0], cn_columns[1]]].min(axis=1)
        cns_df.drop(columns=cn_columns, inplace=True)
        cns_df["major_cn"] = major_cn
        cns_df["minor_cn"] = minor_cn
        cn_columns = ["major_cn", "minor_cn"]
        log_info(print_info, f"Converted columns to ordered")
    elif input_cn_columns == None:
        cns_df, cn_columns = rename_cn_cols(cns_df, cn_columns, print_info)

    select_cols = ["sample_id", "chrom", "start", "end"] + cn_columns
    if "name" in cns_df.columns:
        select_cols.append("name")

    cns_df = cns_df[select_cols]
    # set dtypes
    col_types = {"sample_id": str, "chrom": str, "start": int, "end": int}
    for cn_columns in cn_columns:
        col_types[cn_columns] = float
    cns_df = cns_df.astype(col_types)
    return cns_df


def _is_cn_column(column):
    if not isinstance(column, str):
        return False
    pattern = re.compile(r'^(cn).*|.*(cn)$|.*(major|minor|hap|total|allele).*', re.IGNORECASE)
    return bool(re.search(pattern, column))

