#!/bin/bash

# Load the FastQC module
module load apps/fastqc/0.11.9   # FastQC version 0.11.9

# Change to the directory where the FASTQ files are located
cd <path_to_trimmed_files_directory>

# Create the directory for FastQC results if it doesn't exist
mkdir -p <path_to_fastqc_results_directory>

# Run FastQC on each *_1_trimmed_paired.fastq.gz file and save results in the specified directory
for file in *_1_trimmed_paired.fastq.gz; do
    fastqc -o <path_to_fastqc_results_directory> "$file"
done

# Run FastQC on each *_2_trimmed_paired.fastq.gz file and save results in the specified directory
for file in *_2_trimmed_paired.fastq.gz; do
    fastqc -o <path_to_fastqc_results_directory> "$file"
done

