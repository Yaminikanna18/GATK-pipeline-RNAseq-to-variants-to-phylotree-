#!/bin/bash
#SBATCH --job-name=<job_name>  # Job name
#SBATCH --output=<job_name>_%j.out  # Output file
#SBATCH --error=<job_name>_%j.err  # Error log file
#SBATCH --nodes=1  # Number of nodes
#SBATCH --ntasks-per-node=1  # Number of tasks (processes) per node
#SBATCH --cpus-per-task=<num_cpus>  # Number of CPU cores per task
#SBATCH --time=<time_limit>  # Time limit hours:minutes:seconds
#SBATCH --mem=<memory_limit>  # Memory limit per node
#SBATCH --account=<account_name>  # Specify your Slurm account here

# Source Conda initialization script
source /path/to/conda.sh

# Activate Conda environment
conda activate <your_conda_environment>

# Print environment information for debugging
echo "Current PATH: $PATH"

# Verify STAR and Samtools are available
if ! command -v STAR &> /dev/null; then
    echo "STAR could not be found"
    exit 1
fi

if ! command -v samtools &> /dev/null; then
    echo "Samtools could not be found"
    exit 1
fi

# Directly test the STAR and Samtools executables
STAR --version
samtools --version

# Set working directory
WORK_DIR=/path/to/working_directory
INDEX_DIR=$WORK_DIR/star_index
READS_DIR=$WORK_DIR
OUTPUT_DIR=$WORK_DIR/star_output

# Ensure output directory exists and is writable
if [ ! -d "$OUTPUT_DIR" ]; then
    echo "Creating output directory: $OUTPUT_DIR"
    mkdir -p "$OUTPUT_DIR"
fi

if [ ! -w "$OUTPUT_DIR" ]; then
    echo "No write permissions for $OUTPUT_DIR"
    exit 1
fi

# Run STAR alignment for each pair of read files
READ1=$1
READ2=$2
BASENAME=$(basename $READ1 _1_trimmed_paired.fastq.gz)
OUTPUT_PREFIX=$OUTPUT_DIR/${BASENAME}_

STAR --runThreadN 8 \
     --genomeDir $INDEX_DIR \
     --readFilesIn <(zcat $READS_DIR/$READ1) <(zcat $READS_DIR/$READ2) \
     --outSAMtype BAM SortedByCoordinate \
     --outFileNamePrefix $OUTPUT_PREFIX

# Check if alignment was successful
if [ $? -ne 0 ]; then
    echo "STAR alignment failed for $READ1 and $READ2!"
    exit 1
fi

