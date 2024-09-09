#!/bin/bash
#SBATCH --job-name=star_indexing  # Job name
#SBATCH --output=star_indexing_%j.out  # Output file (%j expands to job ID)
#SBATCH --error=star_indexing_%j.err  # Error log file (%j expands to job ID)
#SBATCH --nodes=1  # Number of nodes
#SBATCH --ntasks-per-node=1  # Number of tasks (processes) per node
#SBATCH --cpus-per-task=4  # Number of CPU cores per task
#SBATCH --time=2:00:00  # Time limit hours:minutes:seconds
#SBATCH --mem=32G  # Memory limit per node
#SBATCH --account=<account_name>  # Specify your Slurm account here

# Source Conda initialization script
source /user/home/<your_username>/miniforge3/etc/profile.d/conda.sh

# Activate Conda environment
conda activate <your_conda_environment>

# Print environment information for debugging
echo "Current PATH: $PATH"

# Verify STAR is available
if ! command -v STAR &> /dev/null; then
    echo "STAR could not be found"
    exit 1
fi

# Directly test the STAR executable
STAR --version

# Set working directory
WORK_DIR=/path/to/your/working_directory
INDEX_DIR=$WORK_DIR/star_index
GENOME_FASTA=$WORK_DIR/genome.fasta
GTF_FILE=$WORK_DIR/annotations.gtf

# Ensure index directory exists
mkdir -p $INDEX_DIR

# Generate STAR index
STAR --runThreadN <number_of_threads> \
     --runMode genomeGenerate \
     --genomeDir $INDEX_DIR \
     --genomeFastaFiles $GENOME_FASTA \
     --sjdbGTFfile $GTF_FILE \
     --sjdbOverhang <sjdb_overhang> \
     --genomeSAindexNbases <genome_saindex_nbases>

# Check if indexing was successful
if [ $? -ne 0 ]; then
    echo "STAR genome indexing failed!"
    exit 1
fi

echo "STAR genome indexing completed successfully."

