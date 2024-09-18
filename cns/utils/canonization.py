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


def find_cn_cols_if_none(cns_df, cn_cols=None):
    if cn_cols is None:
        cn_cols = [col for col in cns_df.columns if is_cn_column(col)]
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
            raise ValueError("Only one (total) or two (major, minor) CN columns are allowed.")
        if len(cn_cols) == 0:
            raise ValueError("No CN columns found.")
        return cn_cols
    return cn_cols


def rename_cn_cols(cns_df, cn_columns=None, print_info=False):
    cn_columns = find_cn_cols_if_none(cns_df, cn_columns)
    if not _requires_rename(cn_columns):
        return cns_df, list(cn_columns)

    if len(cn_columns) == 2:
        rename_map = _get_major_minor_cols(cns_df, cn_columns)
        if len(rename_map) == 0: # Try haplotype specific
            y_only = cns_df.query("chrom == 'chrY'")
            if (y_only[cn_columns[1]] == 0).all():
                rename_map={cn_columns[0]: "hap1_cn", cn_columns[1]: "hap2_cn"}
            elif (y_only[cn_columns[0]] == 0).all():
                rename_map={cn_columns[0]: "hap2_cn", cn_columns[1]: "hap1_cn"}
            else:
                raise ValueError("Haplotypes are of mixed sex. Either of the two CN columns")
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
        return {col1: "major_cn", col2: "minor_cn"}
    else:
        return {}  


def find_column(cns_df, patterns):
    # Find matching column
    matching_column = None
    for col in cns_df.columns:
        if any(pd.Series(col).str.contains(pattern, case=False, regex=True).any() for pattern in patterns):
            matching_column = col
            break
    return matching_column


def canonize_cns_df(cns_df, cn_columns=None, order_columns=False, assembly=hg19, print_info=False):
    cn_columns = find_cn_cols_if_none(cns_df, cn_columns)

    # if the column sample_id does not exist, rename the first column to sample_id
    if "sample_id" not in cns_df.columns:
        sample_col = find_column(cns_df, ['sample', 'Id', 'sampleId', 'sample_id'])
        if sample_col is None:
            cns_df.columns = ["sample_id"] + cns_df.columns[1:].tolist()
            log_info(print_info, f"Column sample_id not found, renamed first column to sample_id.")
        else:
            cns_df.rename(columns={sample_col: "sample_id"}, inplace=True)
            log_info(print_info, f"Renamed column {sample_col} to sample_id.")

    # if the column chrom does not exist, rename the second column to chrom
    if "chrom" not in cns_df.columns:
        chrom_col = find_column(cns_df, ['chrom', 'chr', 'chromosome'])
        if chrom_col is None:
            cns_df.columns = cns_df.columns[:1].tolist() + ["chrom"] + cns_df.columns[2:].tolist()
            log_info(print_info, f"Column chrom not found, renamed second column to chrom.")
        else:
            cns_df.rename(columns={chrom_col: "chrom"}, inplace=True)
            log_info(print_info, f"Renamed column {chrom_col} to chrom.")

    chrom_vals = cns_df["chrom"].unique()
    # if the chromosomes values are all either digits or single characters, convert to chrX format
    if all([chrom.isdigit() or len(chrom) == 1 for chrom in chrom_vals]):
        cns_df["chrom"] = "chr" + cns_df["chrom"]
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
        rows_to_drop = cns_df[cns_df["chrom"].isin(chrom_vals[not_known])].index
        cns_df.drop(rows_to_drop, inplace=True) 

    # if the column start does not exist, rename the third column to start
    if "start" not in cns_df.columns:
        start_col = find_column(cns_df, ['start', 'begin'])
        if start_col is None:
            cns_df.columns = cns_df.columns[:2].tolist() + ["start"] + cns_df.columns[3:].tolist()
            log_info(print_info, f"Column start not found, renamed third column to start.")
        else:
            cns_df.rename(columns={start_col: "start"}, inplace=True)
            log_info(print_info, f"Renamed column {start_col} to start.")
    cns_df["start"] = cns_df["start"].astype(int)

    # if the column end does not exist, rename the fourth column to end
    if "end" not in cns_df.columns:
        end_col = find_column(cns_df, ['end', 'stop'])
        if end_col is None:
            cns_df.columns = cns_df.columns[:3].tolist() + ["end"] + cns_df.columns[4:].tolist()
            log_info(print_info, f"Column end not found, renamed fourth column to end.")
        else:
            cns_df.rename(columns={end_col: "end"}, inplace=True)
            log_info(print_info, f"Renamed column {end_col} to end.")

    cn_columns = find_cn_cols_if_none(cns_df, cn_columns)  
    log_info(print_info, f"Using CN columns: {cn_columns}")
    if len(cn_columns) == 2 and order_columns:
        major_cn = cns_df[[cn_columns[0], cn_columns[1]]].max(axis=1)
        minor_cn = cns_df[[cn_columns[0], cn_columns[1]]].min(axis=1)
        cns_df.drop(columns=cn_columns, inplace=True)
        cns_df["major_cn"] = major_cn
        cns_df["minor_cn"] = minor_cn
        cn_columns = ["major_cn", "minor_cn"]
        log_info(print_info, f"Converted columns to ordered")
    else:
        cns_df, cn_columns = rename_cn_cols(cns_df, cn_columns)

    cns_df = cns_df[["sample_id", "chrom", "start", "end"] + cn_columns]
    return cns_df


def is_hap_spec(cn_columns):
    return len(cn_columns) == 2 and all([cn_col in ["hap1_cn", "hap2_cn"] for cn_col in cn_columns])


def is_cn_column(column):
    if not isinstance(column, str):
        return False
    pattern = re.compile(r'^(cn|CN).*|.*(cn|CN)$|.*(major|minor|hap).*', re.IGNORECASE)
    return bool(re.search(pattern, column))


def is_canonical_cns_df(cns_df):
    return "sample_id" in cns_df.columns and "chrom" in cns_df.columns and "start" in cns_df.columns and "end" in cns_df.columns and any([is_cn_column(col) for col in cns_df.columns])