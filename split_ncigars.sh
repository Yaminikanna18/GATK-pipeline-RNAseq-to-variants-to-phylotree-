#!/bin/bash
#SBATCH --job-name=split_ncigars                # Job name
#SBATCH --output=<output_directory>/split_ncigars_%j.out  # Output file
#SBATCH --error=<output_directory>/split_ncigars_%j.err   # Error log file
#SBATCH --nodes=1                             # Number of nodes
#SBATCH --ntasks=1                            # Number of tasks (processes) per node
#SBATCH --cpus-per-task=4                     # Number of CPU cores per task
#SBATCH --mem=<memory_limit>                  # Memory limit per node
#SBATCH --time=<time_limit>                   # Time limit hours:minutes:seconds
#SBATCH --account=<account_name>              # Specify your Slurm account here

# Load GATK module
module load gatk/<gatk_version>

# Define paths to files
INPUT_BAM="$1"
OUTPUT_BAM="$2"
REFERENCE="$3"

# Check if the input BAM file exists
if [ ! -f ${INPUT_BAM} ]; then
    echo "Input BAM file not found: ${INPUT_BAM}"
    exit 1
fi

# Run GATK SplitNCigarReads
gatk --java-options "-Xmx28G" SplitNCigarReads \
    -R ${REFERENCE} \
    -I ${INPUT_BAM} \
    -O ${OUTPUT_BAM}

# Check if the output BAM was created successfully
if [ -f ${OUTPUT_BAM} ]; then
    echo "SplitNCigarReads completed successfully."
else
    echo "SplitNCigarReads failed. Output BAM file not found."
fi

