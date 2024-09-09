#!/bin/bash

# Load Java and Trimmomatic modules
module load java/1.8.0_352
module add apps/trimmomatic/0.39

# Input directory containing fastq.gz files
fastq_dir="<input_directory>"

# Output directory to save results
output_dir="<output_directory>"

# Create the output directory if it doesn't exist
mkdir -p "${output_dir}"

# Loop through each fastq.gz file in the directory
for fastq_file in "${fastq_dir}"/*_1.fastq.gz; do
    # Extract sample name from the file name
    sample_name=$(basename "${fastq_file}" "_1.fastq.gz")
    # Run Trimmomatic
    srun --ntasks-per-node=1 java -jar /sw/apps/Trimmomatic-0.39/trimmomatic-0.39.jar PE \
        "${fastq_dir}/${sample_name}_1.fastq.gz" "${fastq_dir}/${sample_name}_2.fastq.gz" \
        "${output_dir}/${sample_name}_1_trimmed_paired.fastq.gz" "${output_dir}/${sample_name}_1_trimmed_unpaired.fastq.gz" \
        "${output_dir}/${sample_name}_2_trimmed_paired.fastq.gz" "${output_dir}/${sample_name}_2_trimmed_unpaired.fastq.gz" \
        ILLUMINACLIP:/sw/apps/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10:2:keepBothReads \
        LEADING:3 TRAILING:3 MINLEN:55
    echo "Finished processing ${sample_name}"
done

