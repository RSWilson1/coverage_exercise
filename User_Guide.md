# Running `coverage_report.py` with Command Line Arguments

## Introduction
This script generates a report listing any genes that have less than 100% coverage at 30x.
It takes the sambamba output as input and
amalgamates the coverage by exon to determine the coverage for each gene.

The script generates two excel files, a summary excel for all genes that failed the 30X in at least one exon and summary statistics across all the exons (min, max, mean and median).
The detailed excel report contains all the exons for each gene which failed the 30x coverage for further inspection.

## Prerequisites
Full requirements listed in requirements.txt
- Python 3.12
- Pandas
- openpyxl


## Installation
To install the required dependencies, run the following command:
```
pip install requirements.txt
```

## Usage
Example cmd for running.
```
python coverage_report.py -i NGS148_34_139558_CB_CMCMD_S33_R1_001.sambamba_output.txt -o NGS148
```

## Options
- List all available command line options and their descriptions.
- Include any default values or required arguments.

## Examples
- Provide additional examples of running the script with different combinations of command line arguments.


## References
- List any external resources or references used in creating the document.
