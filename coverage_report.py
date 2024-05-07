"""
This script generates a report listing any genes that have less than threshold coverage at 30x.
It takes the sambamba output as input and
amalgamates the coverage by exon to determine the coverage for each gene.

Usage: python coverage_report.py <sambamba_input_tsv>
Alternatively: python coverage_report.py <directory of tsvs>

# TODO: Additional feature in HGNC IDs for all gene symbols.
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
        files (list): A list of all the tsv files in the directory.
    """
    files = []
    print(directory)
    print(os.listdir(directory))
    # Loop through all the files in the directory
    for file in os.listdir(directory):
        print(file)
        if file.endswith("sambamba_output.txt") or file.endswith("sambamba_output.tsv"):
            files.append(file)
    return files


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
    """
    # Read the sambamba input file, with tab and whitespace as separators
    sambamba_df = pd.read_csv(sambamba_input_file, sep=r"\s+", header=0, index_col=False)
    # Set dtypes for the columns
    dtypes = {
        'StartPosition': "int64",
        'EndPosition': "int64",
        'FullPosition': str,
        'GeneSymbol;Accession': str,
        'Size': "int64",
        'readCount': "int64",
        'meanCoverage': float,
        'percentage30': float,
        'sampleName': str,
    }
    sambamba_df = sambamba_df.astype(dtypes)
    print(sambamba_df.head())
    print(sambamba_df.columns)
    print(sambamba_df.dtypes)
    # Manipulate the dataframe to get the gene symbol and coverage
    sambamba_df["GeneSymbol"] = sambamba_df["GeneSymbol;Accession"].str.split(
        ";").str[0]
    sambamba_df["Coverage"] = sambamba_df["percentage30"]
    # Get accession
    sambamba_df["Accession"] = sambamba_df["GeneSymbol;Accession"].str.split(
        ";").str[1]
    print(sambamba_df.head())

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
        genes_low_coverage_df (dataframe): A dataframe containing the genes
            with less than threshold coverage at 30x.
        genes_low_coverage_list (list):
            A list of all the genes with less than threshold coverage at 30x.

    STOUT:
        Print statement for all genes with less than threshold coverage at 30x.
    """
    # Filter genes with less than threshold coverage at 30x
    exons_low_coverage_df = sambamba_df[sambamba_df['percentage30'] < threshold]

    # Create a dataframe for all genes with under threshold coverage
    exons_low_coverage_df = exons_low_coverage_df[[
        'GeneSymbol', 'Coverage']]

    # Print the genes with less than threshold coverage at 30x
    genes_low_coverage_list = exons_low_coverage_df['GeneSymbol'].unique(
    )
    print(f"Genes with less than threshold {threshold} coverage at 30x: {
          ", ".join(genes_low_coverage_list)}.")

    genes_low_coverage_df = sambamba_df[sambamba_df['GeneSymbol'].isin(
        genes_low_coverage_list
    )]

    return genes_low_coverage_df, genes_low_coverage_list


def write_output_excel(genes_low_coverage_df, output_file_prefix):
    """
    Writes the output to an excel file.

    parameters
    ----------
        genes_low_coverage_df (dataframe): A dataframe containing the genes with less than threshold coverage at 30x.
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

    # Group by 'GeneSymbol' and aggregate 'Coverage'
    summary_genes_low_coverage_df = \
        genes_low_coverage_df.groupby('GeneSymbol')['Coverage'].agg(
            LowestCoverage='min',
            HighestCoverage='max',
            AverageCoverage='mean',
            MedianCoverage='median'
        ).reset_index()

    summary_genes_low_coverage_df.to_excel(
        f"{output_file_prefix}_summary_report.xlsx", index=False
    )
    # write a detailed excel
    genes_low_coverage_df.to_excel(
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
    low_coverage_df, genes_low_coverage_list = check_coverage(sambamba_df, threshold)
    write_output_excel(low_coverage_df, f"{
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
        for file in files:
            generate_single_report(file,
                                   args.output_file_prefix,
                                   args.threshold
                                   )
    # If single file
    elif args.sambamba_input_file:
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
