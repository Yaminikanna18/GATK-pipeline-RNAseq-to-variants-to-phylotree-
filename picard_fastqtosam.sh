#!/bin/bash
#SBATCH --job-name=picard_fastqtosam    # Job name
#SBATCH --output=picard_fastqtosam_%j.out    # Output file (%j expands to job ID)
#SBATCH --error=picard_fastqtosam_%j.err     # Error log file (%j expands to job ID)
#SBATCH --nodes=1                   # Number of nodes
#SBATCH --ntasks-per-node=1          # Number of tasks (processes) per node
#SBATCH --cpus-per-task=4            # Number of CPU cores per task
#SBATCH --time=24:00:00               # Time limit hours:minutes:seconds
#SBATCH --mem=32G                    # Memory limit per node
#SBATCH --account=<account_name>     # Specify your Slurm account here

# Load or activate Conda environment
source /user/home/dq23100/miniforge3/bin/activate <your_conda_environment>

# Define variables
PICARD_JAR="/user/home/dq23100/miniforge3/envs/<your_conda_environment>/share/picard-3.2.0-0/picard.jar"  # Path to the Picard JAR file
INPUT_DIR="/path/to/your/input_directory"  # Path to input directory
OUTPUT_DIR="/path/to/your/output_directory"  # Path to output directory

# Array of base filenames (without _1 and _2)
base_filenames=(
    "sample1"
    "sample2"
    "sample3"
    "sample4"
    "sample5"
)

# Loop through each base filename and run Picard FastqToSam
for base in "${base_filenames[@]}"; do
    echo "Processing ${base}..."

    # Run Picard FastqToSam command
    java -jar ${PICARD_JAR} FastqToSam \
        FASTQ=${INPUT_DIR}/${base}_1_trimmed_paired.fastq.gz \
        FASTQ2=${INPUT_DIR}/${base}_2_trimmed_paired.fastq.gz \
        OUTPUT=${OUTPUT_DIR}/${base}_unmapped.bam \
        READ_GROUP_NAME=rg_${base} \
        SAMPLE_NAME=sample_${base} \
        LIBRARY_NAME=lib_${base} \
        PLATFORM_UNIT=unit_${base} \
        PLATFORM=ILLUMINA

    echo "Completed ${base}."
done

# Print completion message
echo "All jobs completed."

