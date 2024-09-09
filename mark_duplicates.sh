#!/bin/bash
#SBATCH --job-name=mark_duplicates      # Job name
#SBATCH --output=<output_directory>/mark_duplicates_%j.out  # Output file
#SBATCH --error=<output_directory>/mark_duplicates_%j.err   # Error log file
#SBATCH --nodes=1                       # Number of nodes
#SBATCH --ntasks=1                      # Number of tasks (processes) per node
#SBATCH --cpus-per-task=4               # Number of CPU cores per task
#SBATCH --mem=<memory_limit>            # Memory limit per node
#SBATCH --time=<time_limit>             # Time limit hours:minutes:seconds
#SBATCH --account=<account_name>        # Specify your Slurm account here

# Load Conda environment
source ~/miniforge3/bin/activate <conda_environment>

# Define paths to files
INPUT_BAM="$1"
OUTPUT_BAM="$2"
METRICS_FILE="$3"

# Check if the input BAM file exists
if [ ! -f ${INPUT_BAM} ]; then
    echo "Input BAM file not found: ${INPUT_BAM}"
    exit 1
fi

# Run Picard MarkDuplicates
~/miniforge3/envs/<conda_environment>/bin/picard MarkDuplicates \
    INPUT=${INPUT_BAM} \
    OUTPUT=${OUTPUT_BAM} \
    METRICS_FILE=${METRICS_FILE} \
    CREATE_INDEX=true \
    VALIDATION_STRINGENCY=SILENT \
    REMOVE_DUPLICATES=false

# Check if the output BAM was created successfully
if [ -f ${OUTPUT_BAM} ]; then
    echo "MarkDuplicates completed successfully."

    # Remove the merged BAM file
    rm -f ${INPUT_BAM}
    echo "Removed the input BAM file: ${INPUT_BAM}"
else
    echo "MarkDuplicates failed. Output BAM file not found."
fi

