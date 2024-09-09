#!/bin/bash
#SBATCH --job-name=<job_name>               # Job name
#SBATCH --output=<output_directory>/raxml_%j.out  # Standard output log
#SBATCH --error=<output_directory>/raxml_%j.err   # Standard error log
#SBATCH --nodes=1                         # Number of nodes
#SBATCH --ntasks=1                        # Number of tasks (processes)
#SBATCH --cpus-per-task=4                 # Number of CPU cores per task
#SBATCH --mem=16G                         # Memory per node
#SBATCH --time=4:00                       # Time limit
#SBATCH --account=<account_name>          # Project account

# Load necessary modules
module load openmpi/<openmpi_version>
module load raxml/<raxml_version>

# Run RAxML to construct a phylogenetic tree
raxmlHPC -s <input_fasta> -n <output_tree_name> -m GTRGAMMA -p <random_seed> -# <bootstrap_replicates> -x <random_seed>

