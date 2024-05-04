"""
This script generates a report listing any genes that have less than 100% coverage at 30x.
It takes the sambamba output as input and
amalgamates the coverage by exon to determine the coverage for each gene.

Usage: python coverage_report.py <sambamba_input_tsv>
Alternatively: python coverage_report.py <directory of tsvs>

"""

import sys
import os
import argparse
import pandas as pd


def parse_args():
    """
    Parses the command line arguments.

    Returns:
        args (argparse object): The argparse object containing user inputs.
    """
    parser = argparse.ArgumentParser(
        description="Generates a report listing any genes that have less than 100% coverage at 30x.")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-i", "--sambamba_input_file",
                       help="The path to the sambamba input file")
    group.add_argument("-D", "--input_directory",
                       help="The directory containing all the tsvs for input")
    parser.add_argument("-o", "--output_file_prefix",
                        help="The output file prefix that will be prepended to the output file name", required=False)
    args = parser.parse_args()
    return args


def find_files(directory):
    """
    Finds all the tsv files in the directory.

    Args:
        directory (str): The path to the directory.

    Returns:
        files (list): A list of all the tsv files in the directory.
    """
    files = []
    for file in os.listdir(directory):
        if file.endswith("sambamba_output.tsv"):
            files.append(file)
    return files


def read_sambamba_input(sambamba_input_file):
    """
    Reads the sambamba input file and returns a dictionary containing the coverage by exon for each gene.

    Args:
        sambamba_input_file (str): The path to the sambamba input file.

    Returns:
        sambamba_df (df): A dataframe containing the coverage by exon for each gene.
    """

    header_list = ["chromosome", "StartPosition", "EndPosition", "FullPosition", "NotUsed1", "NotUsed2",
                   "GeneSymbol;Accession", "Size", "readCount", "meanCoverage", "percentage30", "sampleName"]
    sambamba_df = pd.read_csv(sambamba_input_file, sep='\t', header=None,
                              names=header_list, comment="#")

    # Manipulate the dataframe to get the gene symbol and coverage
    sambamba_df["GeneSymbol"] = sambamba_df["GeneSymbol;Accession"].str.split(
        ";").str[0]
    sambamba_df["Coverage"] = sambamba_df["percentage30"]
    # Get accession
    sambamba_df["Accession"] = sambamba_df["GeneSymbol;Accession"].str.split(
        ";").str[1]
    print(sambamba_df.head())

    return sambamba_df


def check_coverage(sambamba_df):
    """
    Checks the coverage of all genes that have less than 100% coverage at 30x.
    Generates a subset df containing the genes with less than 100% coverage at 30x.

    Args:
        sambamba_df (dataframe): A dataframe containing
        the coverage by exon for each gene.

    Returns:
        under_100_coverage_df (dataframe): A dataframe containing the genes
        with less than 100% coverage at 30x.

    STOUT:
        Print statement for all genes with less than 100% coverage at 30x.
    """
    # Filter genes with less than 100% coverage at 30x
    filtered_df = sambamba_df[sambamba_df['percentage30'] < 100]

    # Create a dataframe for all genes with under 100% coverage
    under_100_coverage_df = filtered_df[['GeneSymbol', 'Coverage']]

    # Print the genes with less than 100% coverage at 30x
    under_100_coverage_list = under_100_coverage_df['GeneSymbol'].unique()
    print(f"Genes with less than 100% coverage at 30x: {
          ", ".join(under_100_coverage_list)}.")

    print(under_100_coverage_df)

    return under_100_coverage_df


def write_output_excel(under_100_coverage_df, output_file_prefix):
    """
    Writes the output to an excel file.

    Args:
        under_100_coverage_df (dataframe): A dataframe containing the genes with less than 100% coverage at 30x.
        output_file_prefix (str): The path to the output file.

    Returns:
        None

    Outputs:
        excel_report (excel): An excel file containing the genes with less than 100% coverage at 30x.
    """
    # generate excel with gene name header, coverage per exon, and accession
    # Add in HGNC ID if time.
    under_100_coverage_df.to_excel(output_file_prefix, index=False)


def generate_single_report(sambamba_input_file, output_file_prefix):
    """
    Function to generate a report listing any genes that have less than 100% coverage at 30x.

    Args:
        args.sambamba_input_file (str): The path to the sambamba input file.

    Returns:
        None

    Outputs:
        "*.xlsx" (excel): An excel file containing the genes with less than 100% coverage at 30x.
    """
    sambamba_df = read_sambamba_input(sambamba_input_file)
    under_100_coverage_df = check_coverage(sambamba_df)
    write_output_excel(under_100_coverage_df, f"{output_file_prefix}.xlsx")


def main(args):
    """
    Main function to generate a report listing any genes that have less than 100% coverage at 30x.

    Args:
        args (argparse object): The argparse object containing user inputs.

    Returns:
        None

    Outputs:
        "*.xlsx" (excel): An excel file containing the genes with less than 100% coverage at 30x.
    """
    # logic to check which inputs and logic for what to run.
    # If directory
    if args.input_directory:
        files = find_files(args.input_directory)
        for file in files:
            generate_single_report(file, args.output_file_prefix)
    # If single file
    elif args.sambamba_input_file:
        generate_single_report(args.sambamba_input_file,
                               args.output_file_prefix)
    # If multiple files
    # needs more logic to handle multiple files
    else:
        raise RuntimeError(
            "Please provide either a single file or a directory of files.")


if __name__ == "__main__":
    args = parse_args()

    # run main function
    main(args)
