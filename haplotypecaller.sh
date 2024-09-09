#!/bin/bash
#SBATCH --job-name=haplotypecaller                # Job name
#SBATCH --output=<output_directory>/haplotypecaller_%j.out  # Output file
#SBATCH --error=<output_directory>/haplotypecaller_%j.err   # Error log file
#SBATCH --nodes=1                               # Number of nodes
#SBATCH --ntasks-per-node=1                     # Number of tasks (processes) per node
#SBATCH --cpus-per-task=4                       # Number of CPU cores per task
#SBATCH --time=<time_limit>                     # Time limit hours:minutes:seconds
#SBATCH --mem=<memory_limit>                    # Memory limit per node
#SBATCH --account=<account_name>                # Specify your Slurm account here

# Load GATK module
module load gatk/<gatk_version>

# Define paths
REFERENCE="<reference_genome>"
BAM_FILE="<input_bam>"
OUTPUT_GVCF="<output_gvcf>"
INDEX_FILE="${BAM_FILE}.bai"

# Print paths and check files
echo "BAM file: ${BAM_FILE}"
echo "Index file: ${INDEX_FILE}"
echo "Output GVCF: ${OUTPUT_GVCF}"

# Ensure the BAM file is indexed
if [ -f "${INDEX_FILE}" ]; then
    echo "Index file found."
else
    echo "Index file not found. Creating index file."
    # Make sure samtools is available; it might be loaded by default or require a module
    module load samtools
    samtools index ${BAM_FILE}
fi

# Ensure the output directory exists
mkdir -p $(dirname ${OUTPUT_GVCF})

# Run HaplotypeCaller to generate GVCF
gatk HaplotypeCaller \
    -R ${REFERENCE} \
    -I ${BAM_FILE} \
    -O ${OUTPUT_GVCF} \
    --emit-ref-confidence GVCF \
    --annotation FisherStrand \
    --annotation QualByDepth

echo "HaplotypeCaller done"

