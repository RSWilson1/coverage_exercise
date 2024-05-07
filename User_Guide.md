## Generate a report of all genes with low coverage from sambamba files.

## Introduction
This script `coverage_report.py` generates a report listing any genes that have less than 100% coverage at 30x.
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

Or for a directory of sambamba input files.
```
python coverage_report.py -D ./bsambamba_files/ -o directory_example -t 100
```

## Options
#### Required
- `-i`, `--sambamba_input_file`: The path to the single sambamba input file
- `-D`, `--input_directory`: The directory containing all the sambamba input files (must end in `sambamba_output.txt` or `sambamba_output.tsv`)
#### Optional
- `-o`, `--output_file_prefix`: The output file prefix that will be prepended to the output file name (optional)
- `-t`, `--threshold`: The threshold for coverage, default is 100%
