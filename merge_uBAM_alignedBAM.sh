#!/bin/bash
#SBATCH --job-name=<job_name>
#SBATCH --output=<output_directory>/merge_%j.out
#SBATCH --error=<output_directory>/merge_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=<cpus_per_task>
#SBATCH --mem=<memory_limit>
#SBATCH --time=<time_limit>
#SBATCH --account=<account_name>

# Load Conda and activate environment
module load anaconda
source ~/miniforge3/bin/activate <conda_environment>

# Define paths to files
UNMAPPED_BAM=$1
ALIGNED_BAM=$2
OUTPUT_BAM=$3
REFERENCE_FASTA="<reference_fasta>"

# Run Picard MergeBamAlignment
picard MergeBamAlignment \
    ALIGNED=${ALIGNED_BAM} \
    UNMAPPED=${UNMAPPED_BAM} \
    O=${OUTPUT_BAM} \
    R=${REFERENCE_FASTA} \
    SORT_ORDER=coordinate CREATE_INDEX=true CLIP_ADAPTERS=true ADD_MATE_CIGAR=true

# Check if the output BAM was created successfully
if [ -f ${OUTPUT_BAM} ]; then
    echo "Merge completed successfully. Deleting original BAM files."
    rm -f ${UNMAPPED_BAM}
    rm -f ${ALIGNED_BAM}
else
    echo "Merge failed. Output BAM file not found."
fi

