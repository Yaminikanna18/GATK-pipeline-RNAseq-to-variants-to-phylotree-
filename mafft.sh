#!/bin/bash
#SBATCH --job-name=<job_name>            # Job name
#SBATCH --output=<output_directory>/mafft_%j.out  # Standard output log
#SBATCH --error=<output_directory>/mafft_%j.err   # Standard error log
#SBATCH --nodes=1                       # Number of nodes
#SBATCH --ntasks=1                      # Number of tasks (processes)
#SBATCH --cpus-per-task=4               # Number of CPU cores per task
#SBATCH --mem=<memory_limit>            # Memory per node
#SBATCH --time=<time_limit>             # Time limit (hh:mm:ss)
#SBATCH --account=<account_name>        # Project account

# Load the MAFFT module
module load mafft/<mafft_version>

# Specify input and output files
INPUT_FASTA="<path_to_input_fasta>"
OUTPUT_ALN="<path_to_output_aligned_fasta>"

# Run MAFFT with 8 threads and auto selection of alignment method
mafft --thread 8 --auto $INPUT_FASTA > $OUTPUT_ALN

# End of script

