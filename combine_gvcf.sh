#!/bin/bash
#SBATCH --job-name=combine_gvcfs                # Job name
#SBATCH --output=<output_directory>/combine_gvcfs_%j.out  # Output file
#SBATCH --error=<output_directory>/combine_gvcfs_%j.err   # Error log file
#SBATCH --nodes=1                               # Number of nodes
#SBATCH --ntasks=1                             # Number of tasks (processes) per node
#SBATCH --cpus-per-task=4                       # Number of CPU cores per task
#SBATCH --mem=<memory_limit>                    # Memory limit per node
#SBATCH --time=<time_limit>                     # Time limit hours:minutes:seconds
#SBATCH --account=<account_name>                # Specify your Slurm account here

# Load GATK module
module load gatk/4.5.0.0

# Define directories and files
GVCF_DIR="<gvcf_directory>"
REFERENCE_DIR="<reference_directory>"
OUTPUT_DIR="<output_directory>"

REFERENCE_FASTA="${REFERENCE_DIR}/<reference_genome>"
COMBINED_GVCF="${OUTPUT_DIR}/combined.g.vcf"

# Create a file list for CombineGVCFs
FILE_LIST="${OUTPUT_DIR}/file_list.txt"
awk -v dir="${GVCF_DIR}" '{print "--variant " dir "/" $1 "_split_raw.g.vcf.gz"}' <sample_name_map.txt> ${FILE_LIST}

# Check if file_list.txt was created successfully
if [ ! -f ${FILE_LIST} ]; then
  echo "Error: File list ${FILE_LIST} was not created." >&2
  exit 1
fi

# Debugging: Print the content of the file list
echo "Contents of ${FILE_LIST}:"
cat ${FILE_LIST}

# Check if all GVCF files exist
while IFS= read -r line; do
  file_path="${line#--variant }"
  if [ ! -f "${file_path}" ]; then
    echo "Error: GVCF file ${file_path} does not exist." >&2
    exit 1
  fi
done < ${FILE_LIST}

# Combine GVCFs
echo "Combining GVCFs..."
gatk CombineGVCFs \
  -R ${REFERENCE_FASTA} \
  $(cat ${FILE_LIST}) \
  -O ${COMBINED_GVCF}

# Check if CombineGVCFs completed successfully
if [ $? -ne 0 ]; then
  echo "Error: CombineGVCFs failed." >&2
  exit 1
fi

# Success message
echo "CombineGVCFs step completed successfully."

