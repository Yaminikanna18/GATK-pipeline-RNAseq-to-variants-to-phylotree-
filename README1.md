# FastQC Analysis Pipeline

This repository contains a SLURM script for running FastQC on FASTQ files.

## Script Details

- **fastqc_script.sh**: A SLURM script to run FastQC on FASTQ files. The script performs quality checks on files with suffixes `_1.fastq.gz` and `_2.fastq.gz`.

## Usage

1. **Update SLURM Details**: Before using the script, make sure to modify the placeholders in the `fastqc_script.sh` script to include the correct SLURM job parameters and directories.

2. **Modify Directories**: Replace the placeholders in the `fastqc_script.sh` script with the appropriate directory paths for your FASTQ files and output results.

3. **Submit the Script**: Submit the script to the SLURM scheduler using the command:

   ```bash
   sbatch fastqc_script.sh

