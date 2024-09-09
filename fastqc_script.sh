#!/bin/bash

# Load the FastQC module
module load apps/fastqc/0.11.9

# Change to the directory where the FASTQ files are located
cd /user/home/dq23100/Scratch/PRJNA599452

# Create the directory for FastQC results if it doesn't exist
mkdir -p /user/home/dq23100/Scratch/PRJNA599452/fastqc_results

# Run FastQC on each FASTQ file ending in _1.fastq.gz and _2.fastq.gz separately
for file in *_1.fastq.gz; do
    fastqc -o /user/home/dq23100/Scratch/PRJNA599452/fastqc_results $file
done

for file in *_2.fastq.gz; do
    fastqc -o /user/home/dq23100/Scratch/PRJNA599452/fastqc_results $file
done

