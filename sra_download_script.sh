#!/bin/bash

#SBATCH --job-name=<job_name>                      # Replace with your job name (e.g., "sra_download")
#SBATCH --output=<output_directory>/sra_download_%j.out  # Output directory placeholder
#SBATCH --error=<output_directory>/sra_download_%j.err   # Error directory placeholder
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=<time_limit>                        # Replace with desired time limit (e.g., "24:00:00")
#SBATCH --mem=<memory_limit>                       # Replace with memory requirement (e.g., "32G")
#SBATCH --account=<account_name>                   # Replace with your account name

# Load the SRA Toolkit module
module load apps/sratoolkit/<sra_toolkit_version>   # Replace with the appropriate version of SRA Toolkit

# Define your SRR accession numbers here as placeholders
SRR_ACC_LIST=(
<SRR_ACC1> <SRR_ACC2> <SRR_ACC3> ...               # Replace or allow users to insert their own SRR accession numbers
)

# Set output directory
OUTPUT_DIR=<output_directory>                      # Replace with a path or allow users to set their directory

# Change to the output directory
cd $OUTPUT_DIR

# Iterate over each accession ID
for SRR_ACC in "${SRR_ACC_LIST[@]}"; do
    # Download SRA files and convert to fastq format
    fasterq-dump --split-files $SRR_ACC

    # Gzip the fastq files
    gzip ${SRR_ACC}_1.fastq
    gzip ${SRR_ACC}_2.fastq
done

