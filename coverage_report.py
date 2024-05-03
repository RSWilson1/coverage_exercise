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
        df: A dataframe containing the coverage by exon for each gene.
    """

    header_list = ["chromosome", "StartPosition", "EndPosition", "FullPosition", "NotUsed1", "NotUsed2", "GeneSymbol;Accession", "Size", "readCount", "meanCoverage", "percentage30", "sampleName"]
    sambamba_df = pd.read_csv(sambamba_input_file, sep='\t', header=None,
                              names=header_list, comment="#")

    # Manipulate the dataframe to get the gene symbol and coverage
    sambamba_df["GeneSymbol"] = sambamba_df["GeneSymbol;Accession"].str.split(";").str[0]
    sambamba_df["Coverage"] = sambamba_df["percentage30"]
    # Get accession
    sambamba_df["Accession"] = sambamba_df["GeneSymbol;Accession"].str.split(";").str[1]
    print(sambamba_df.head())

    return sambamba_df


def generate_report(sambamba_df):
    """
    Generates a report listing any genes that have less than 100% coverage at 30x.

    Args:
        sambamba_df (dataframe): A dataframe containing the coverage by exon for each gene.

    Returns:
        None

    STOUT:
        Print statement for all genes with less than 100% coverage at 30x.
    """
    # Filter genes with less than 100% coverage at 30x
    filtered_df = sambamba_df[sambamba_df['percentage30'] < 100]

    # Create a dataframe for all genes with under 100% coverage
    under_100_coverage_df = filtered_df[['GeneSymbol', 'Coverage']]

    # Print the genes with less than 100% coverage at 30x
    under_100_coverage_list = under_100_coverage_df['GeneSymbol'].unique()
    print(f"Genes with less than 100% coverage at 30x: {", ".join(under_100_coverage_list)}.")


    print(under_100_coverage_df)
def main(sambamba_input_file):
    """
    Main function to generate a report listing any genes that have less than 100% coverage at 30x.

    Args:
        sambamba_input_file (str): The path to the sambamba input file.

    Returns:
        None
    """
    sambamba_df = read_sambamba_input(sambamba_input_file)
    generate_report(sambamba_df)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generates a report listing any genes that have less than 100% coverage at 30x.")
    parser.add_argument("sambamba_input_file", help="The path to the sambamba input file")
    parser.add_argument("input_directory", help="The directory containing all the tsvs for input")
    parser.add_argument("output_file", help="The path to the output file")
    args = parser.parse_args()

    main(args.sambamba_input_file)