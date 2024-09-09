#!/bin/bash
#SBATCH --job-name=genotype_gvcfs                # Job name
#SBATCH --output=<output_directory>/genotype_gvcfs_%j.out  # Output file
#SBATCH --error=<output_directory>/genotype_gvcfs_%j.err   # Error log file
#SBATCH --nodes=1                               # Number of nodes
#SBATCH --ntasks=1                             # Number of tasks (processes) per node
#SBATCH --cpus-per-task=4                       # Number of CPU cores per task
#SBATCH --mem=<memory_limit>                    # Memory limit per node
#SBATCH --time=<time_limit>                     # Time limit hours:minutes:seconds
#SBATCH --account=<account_name>                # Specify your Slurm account here

# Load GATK module
module load gatk/<gatk_version>

# Define directories and files
REFERENCE_DIR="<reference_directory>"
OUTPUT_DIR="<output_directory>"

REFERENCE_FASTA="${REFERENCE_DIR}/<reference_genome>"
COMBINED_GVCF="${OUTPUT_DIR}/combined.g.vcf"
RAW_VCF="${OUTPUT_DIR}/raw_variants.vcf"

# Perform joint genotyping
echo "Performing joint genotyping..."
gatk GenotypeGVCFs \
  -R ${REFERENCE_FASTA} \
  -V ${COMBINED_GVCF} \
  -O ${RAW_VCF}

# Check if GenotypeGVCFs completed successfully
if [ $? -ne 0 ]; then
  echo "Error: GenotypeGVCFs failed." >&2
  exit 1
fi

# Success message
echo "GenotypeGVCFs step completed successfully."

