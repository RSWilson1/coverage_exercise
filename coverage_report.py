"""
This script generates a report listing any genes that have less than threshold coverage at 30x.
It takes the sambamba output as input and
amalgamates the coverage by exon to determine the coverage for each gene.

Usage: python coverage_report.py -i <sambamba_input_txt> -t 100
Alternatively: python coverage_report.py -D <directory of sambamba files> -t 100

Roadmap for the script:
# TODO: Additional feature in HGNC IDs as gene symbols are not unique.
# TODO: Additional feature to add a column for the number of exons in the gene.
# TODO: Additional feature to add a column for the number of exons with coverage less than threshold.
# TODO: Make more memory efficient by removing unnecessary columns and minimum dtypes.
# TODO: Handling for a list of files.
# TODO: Add testing for functional parts of script
 - parse_args
 - read_sambamba_input
 - single vs directory
 - errors raised.
Written 07/05/2024 by Robert Wilson
"""

import os
import argparse
import pandas as pd


def parse_args():
    """
    Parses the command line arguments.

    Returns
    -------
        args (argparse object): The argparse object containing user inputs.
    """
    parser = argparse.ArgumentParser(
        description="Generates a report listing any genes that have less than threshold coverage at 30x.")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-i", "--sambamba_input_file",
                       help="The path to the sambamba input file")
    group.add_argument("-D", "--input_directory",
                       help="The directory containing all the tsvs for input")
    parser.add_argument("-o", "--output_file_prefix",
                        help="The output file prefix that will be prepended to the output file name", required=False)
    parser.add_argument("-t", "--threshold", type=int,
                        help="The threshold for coverage, default is 100%",
                        required=False, default=100)
    args = parser.parse_args()
    return args


def find_files(directory):
    """
    Finds all the tsv files in the directory.

    parameters
    ----------
        directory (str): The path to the directory.

    Returns
    -------
        files_dict (dict): A dictionary of all the sambamba input files
            in the directory with their associated prefix.
    """
    files_dict = {}

    # Loop through all the files in the directory
    for file in os.listdir(directory):
        if file.endswith("sambamba_output.txt") or file.endswith("sambamba_output.tsv"):
            # Extract the file prefix
            file_prefix = file.split('.')[0]
            # Store the file prefix in the dictionary with the file name as the key
            files_dict[file] = file_prefix
    return files_dict


def read_sambamba_input(sambamba_input_file):
    """
    Reads the sambamba input file and returns a dictionary containing the coverage by exon for each gene.

    parameters
    ----------
        sambamba_input_file (str): The path to the sambamba input file
        containing the coverage by exon for each gene.

    Returns
    -------
        sambamba_df (df): A dataframe containing the coverage by exon for each gene.

    +-------------+---------------+-------------+---------------------+
    | #chromosome | StartPosition | EndPosition |    FullPosition     |
    +-------------+---------------+-------------+---------------------+
    |           1 |      26126711 |    26126914 | 1-26126711-26126914 |
    |           1 |      26127523 |    26127661 | 1-26127523-26127661 |
    |           1 |      26128496 |    26128618 | 1-26128496-26128618 |
    +-------------+---------------+-------------+---------------------+...
    +---------+-----------+----------------------+-------+-----------+
    | NotUsed | NotUsed.1 | GeneSymbol;Accession | Size  | readCount |
    +---------+-----------+----------------------+-------+-----------+
    |       0 | +         | SELENON;NM_020451.2  | 57190 |       142 |
    |       0 | +         | SELENON;NM_020451.2  | 57190 |      5748 |
    |       0 | +         | SELENON;NM_020451.2  | 57190 |      2637 |
    +---------+-----------+----------------------+-------+-----------+...
    +--------------+--------------+------------+------------+-------------+
    | meanCoverage | percentage30 | sampleName | GeneSymbol | Accession   |
    +--------------+--------------+------------+------------+-------------+
    |      39.5222 |      49.2611 |          1 | SELENON    | NM_020451.2 |
    |      3501.21 |          100 |          1 | SELENON    | NM_020451.2 |
    |      1105.45 |          100 |          1 | SELENON    | NM_020451.2 |
    +--------------+--------------+------------+------------+-------------+

    nb. percentage30 is the percentage of the exon covered at 30x.

    """
    # Read the sambamba input file, with tab and whitespace as separators
    try:
        sambamba_df = pd.read_csv(
            sambamba_input_file, sep=r"\s+", header=0, index_col=False)
    except Exception as e:
        raise RuntimeError(f"Error reading the file, see: {e}")
    # Set dtypes for the columns
    dtypes = {
        'StartPosition': "int64",
        'EndPosition': "int64",
        'FullPosition': str,
        'GeneSymbol;Accession': str,
        'Size': "int64",
        'readCount': "int64",
        'meanCoverage': "float64",
        'percentage30': "float64",
        'sampleName': str,
    }
    sambamba_df = sambamba_df.astype(dtypes)

    # Split 'GeneSymbol;Accession' into 'GeneSymbol' and 'Accession'
    sambamba_df[['GeneSymbol', 'Accession']] = sambamba_df[
        "GeneSymbol;Accession"
    ].str.split(';', expand=True)

    return sambamba_df


def check_coverage(sambamba_df, threshold):
    """
    Checks the coverage of all genes that have less than threshold,
    (default 100%) coverage at 30x.
    Generates a subset df containing the genes with less than threshold coverage at 30x.

    parameters
    ----------
        sambamba_df (dataframe): A dataframe containing
        the coverage by exon for each gene.

    Returns
    -------
        low_coverage_genes_df (dataframe): A dataframe containing the genes
            with less than threshold coverage at 30x.
        low_coverage_genes_list (list):
            A list of all the genes with less than threshold coverage at 30x.

    STOUT:
        Print statement for all genes with less than threshold coverage at 30x.
    """
    # Filter genes with less than threshold coverage at 30x
    low_coverage_exons_df = sambamba_df[sambamba_df['percentage30'] < threshold]

    # Create a dataframe for all genes with under threshold coverage
    low_coverage_exons_df = low_coverage_exons_df[[
        'GeneSymbol', 'percentage30']]

    # Print the genes with less than threshold coverage at 30x
    low_coverage_genes_list = low_coverage_exons_df['GeneSymbol'].unique(
    )
    print(f"Genes with less than threshold ({threshold}% coverage) at 30x: {
          ", ".join(low_coverage_genes_list)}.")

    low_coverage_genes_df = sambamba_df[sambamba_df['GeneSymbol'].isin(
        low_coverage_genes_list
    )]

    return low_coverage_genes_df, low_coverage_genes_list


def write_output_excels(low_coverage_genes_df, output_file_prefix):
    """
    Writes the output to an excel file.

    parameters
    ----------
        low_coverage_genes_df (dataframe): A dataframe containing the genes with less than threshold coverage at 30x.
        output_file_prefix (str): The path to the output file.

    Returns
    -------
        None

    Outputs:
        "*_report.xlsx" (excel): An excel file containing
            the genes with less than threshold coverage at 30x.

        "*_summary_report.xlsx" (excel): An excel file containing
            the genes with less than threshold coverage at 30x with summary metrics.
    """
    # generate excel with gene name header, coverage per exon, and accession

    # Group by 'GeneSymbol' and aggregate 'percentage30'
    low_coverage_genes_summary_df = \
        low_coverage_genes_df.groupby('GeneSymbol;Accession')['percentage30'].agg(
            LowestCoverage='min',
            HighestCoverage='max',
            MeanCoveragePerExon='mean',
            MedianCoveragePerExon='median'
        ).reset_index()

    # Split 'GeneSymbol;Accession' into 'GeneSymbol' and 'Accession'
    low_coverage_genes_summary_df[['GeneSymbol', 'Accession']] = \
        low_coverage_genes_summary_df["GeneSymbol;Accession"].str.split(
            ';', expand=True)

    #reorder columns
    low_coverage_genes_summary_df = low_coverage_genes_summary_df[[
        'GeneSymbol', 'Accession', 'LowestCoverage', 'HighestCoverage',
        'MeanCoveragePerExon', 'MedianCoveragePerExon'
    ]]
    low_coverage_genes_summary_df.to_excel(
        f"{output_file_prefix}_summary_report.xlsx", index=False
    )
    # write a detailed excel
    low_coverage_genes_df.to_excel(
        f"{output_file_prefix}_report.xlsx",
        index=False
    )


def generate_single_report(sambamba_input_file, output_file_prefix, threshold):
    """
    Function to generate a report listing any genes that have less than threshold coverage at 30x.

    parameters
    ----------
        args.sambamba_input_file (str): The path to the sambamba input file.

        threshold (int): The threshold for coverage, default is threshold.

    Returns
    -------
        None

    Outputs:
        "*_report.xlsx" (excel): An excel file containing
            the genes with less than threshold coverage at 30x.

        "*_summary_report.xlsx" (excel): An excel file containing
            the genes with less than threshold coverage at 30x with summary metrics.
    """
    sambamba_df = read_sambamba_input(sambamba_input_file)
    low_coverage_df, low_coverage_genes_list = check_coverage(
        sambamba_df, threshold)
    write_output_excels(low_coverage_df, f"{
        output_file_prefix}_report.xlsx")


def main(args):
    """
    Main function to generate a report listing any genes that have less than threshold coverage at 30x.

    parameters
    ----------
        args (argparse object): The argparse object containing user inputs.

    Returns
    -------
        None

    Errors
    ------
        RuntimeError: If the user does not provide either a single file or a directory of files.
    """
    # logic to check which inputs and logic for what to run.
    # If directory
    if args.input_directory:
        files = find_files(args.input_directory)
        for file, file_prefix in files.items():
            generate_single_report(file,
                                   f"{args.output_file_prefix}{file_prefix}",
                                   args.threshold
                                   )
    # If single file
    elif args.sambamba_input_file:
        if args.output_file_prefix is None:
            args.output_file_prefix = args.sambamba_input_file.split(".")[0]
        elif args.output_file_prefix:
            input_file_name = args.sambamba_input_file.split(".")[0]
            args.output_file_prefix = f"{
                args.output_file_prefix}{input_file_name}"
        generate_single_report(args.sambamba_input_file,
                               args.output_file_prefix,
                               args.threshold)
    # If multiple files
    # needs more logic to handle multiple files
    else:
        raise RuntimeError(
            "Please provide either a single file or a directory of files.")


if __name__ == "__main__":
    args = parse_args()

    # run main function
    main(args)
