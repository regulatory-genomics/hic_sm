# Copied from https://github.com/mirnylab/distiller-sm/blob/master/_distiller_common.py
import os, pathlib
import numpy as np
import pandas as pd
import shlex

def argstring_to_dict(argstring):
    """
    Convert a command line argument string into a dictionary.
    
    Args:
        argstring (str): The command line argument string.
        
    Returns:
        dict: A dictionary with argument names as keys and their values.
        
    Note:
        - Arguments that start with '-' are considered keys.
        - If an argument has no value, it is set to True.
        - Values can be separated by spaces and will be joined into a single string.
        - If the final argument is a value without a preceding key, it will be included as a value for the last key (issue in case of specified input as last argument).
    """
    args = shlex.split(argstring)
    keys = np.where([arg.startswith('-') for arg in args])[0]
    vals = np.where([not arg.startswith('-') for arg in args])[0]
    args_arrs = [arr for arr in np.split(args, keys) if arr.size > 0]
    argdict = {str(arr[0]): (' '.join(arr[1:]) if arr.size > 1 else True) for arr in args_arrs }
    return argdict


def parse_fastq_dataframe(df):
    """
    Convert a pandas DataFrame to library_run_fastqs dictionary structure.
    
    Args:
        df (pd.DataFrame): DataFrame with columns: sample_id, lane, fastq1, fastq2
        
    Returns:
        dict: A nested dictionary with structure {library: {run: [fastq1, fastq2]}}
    """
    # Accept either fastq1/fastq2 or R1/R2 column naming
    if {"fastq1", "fastq2"}.issubset(df.columns):
        fastq1_col, fastq2_col = "fastq1", "fastq2"
    elif {"R1", "R2"}.issubset(df.columns):
        fastq1_col, fastq2_col = "R1", "R2"
    else:
        raise Exception(
            "DataFrame is missing required FASTQ columns. Expected either "
            "{fastq1, fastq2} or {R1, R2}."
        )

    required_columns = {"sample_id", "lane", fastq1_col, fastq2_col}
    if not required_columns.issubset(df.columns):
        missing = required_columns - set(df.columns)
        raise Exception(
            f"DataFrame is missing required columns: {missing}. "
            f"Required columns are: {required_columns}"
        )
    
    library_run_fastqs = {}
    
    # Iterate through DataFrame rows and build the nested dictionary
    for _, row in df.iterrows():
        library = str(row['sample_id'])
        run = str(row['lane'])
        fastq1 = str(row[fastq1_col])
        fastq2 = str(row[fastq2_col])
        
        # Initialize nested dictionaries if they don't exist
        if library not in library_run_fastqs:
            library_run_fastqs[library] = {}
        if run not in library_run_fastqs[library]:
            library_run_fastqs[library][run] = []
        
        # Add fastq files
        library_run_fastqs[library][run] = [fastq1, fastq2]
    
    return library_run_fastqs


def organize_fastqs(config):
    load_csv = config["input"]["load_csv"]
    if load_csv:
        fastq_csv = config["input"]["fastq_csv"]
        raw_reads_input = pd.read_csv(fastq_csv)
    else:
        raw_reads_input = config["input"]["raw_reads_paths"]
    
    # Case 2: pandas DataFrame
    elif isinstance(raw_reads_input, pd.DataFrame):
        library_run_fastqs = parse_fastq_dataframe(raw_reads_input)
    
    # Case 3: Dictionary structure
    elif isinstance(raw_reads_input, dict):
        if not check_fastq_dict_structure(raw_reads_input):
            raise Exception(
                "An unknown format for library_fastqs! Please provide it as either "
                'a path to the folder structured as "library/run/fastqs", '
                "a dictionary specifying the project structure, or "
                "a pandas DataFrame with columns: sample_id, lane, fastq1, fastq2."
            )
        library_run_fastqs = raw_reads_input
    
    else:
        raise Exception(
            "An unknown format for library_fastqs! Please provide it as either "
            "a path to the folder with the structure library/run/fastqs, "
            "a dictionary specifying the project structure, or "
            "a pandas DataFrame with columns: sample_id, lane, fastq1, fastq2."
        )

    return library_run_fastqs



def get_units_fastqs(w):
    """
    Resolve input FASTQs for a given sample_run wildcard.
    Expects sample_run to be formatted as '{library}.{run}' if library/run are not present.
    """
    if hasattr(w, "library") and hasattr(w, "run"):
        return get_raw_fastqs(w)
    if hasattr(w, "sample_run"):
        library, run = w.sample_run.split(".")
        fastqs = LIBRARY_RUN_FASTQS[library][run]
        return [
            f"{downloaded_fastqs_folder}/{library}.{run}.1.fastq.gz"
            if needs_downloading(fastqs, 0)
            else fastqs[0],
            f"{downloaded_fastqs_folder}/{library}.{run}.2.fastq.gz"
            if needs_downloading(fastqs, 1)
            else fastqs[1],
        ]
    raise ValueError("sample_run wildcard is required to resolve FASTQs.")

