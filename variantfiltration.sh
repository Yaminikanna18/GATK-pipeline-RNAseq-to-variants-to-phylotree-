#!/bin/bash
#SBATCH --job-name=filter_variants                    # Job name
#SBATCH --output=<output_directory>/filter_variants_%j.out  # Output file
#SBATCH --error=<output_directory>/filter_variants_%j.err   # Error log file
#SBATCH --nodes=1                                   # Number of nodes
#SBATCH --ntasks=1                                 # Number of tasks (processes) per node
#SBATCH --cpus-per-task=4                           # Number of CPU cores per task
#SBATCH --mem=<memory_limit>                        # Memory limit per node
#SBATCH --time=<time_limit>                         # Time limit hours:minutes:seconds
#SBATCH --account=<account_name>                    # Specify your Slurm account here

# Load bcftools module
module load bcftools/<bcftools_version>

# Define input and intermediate output files with full paths
INPUT_VCF="<input_vcf>"
TEMP_VCF_INDEL_REMOVED="<temp_vcf_indel_removed>"
TEMP_VCF1="<temp_vcf_gq_filtered>"
TEMP_VCF2="<temp_vcf_dp_filtered>"
OUTPUT_VCF="<output_vcf>"

# Step 1: Remove indels
bcftools view -V indels -O z -o $TEMP_VCF_INDEL_REMOVED $INPUT_VCF

# Check if the indel-removed VCF file was created successfully
if [ ! -f "$TEMP_VCF_INDEL_REMOVED" ]; then
    echo "Error: Indel-removed VCF file was not created."
    exit 1
fi

# Step 2: Filter variants by Genotype Quality (GQ < 20)
bcftools filter -e 'FORMAT/GQ < 20' -O z -o $TEMP_VCF1 $TEMP_VCF_INDEL_REMOVED

# Check if the GQ filtered VCF file was created successfully
if [ ! -f "$TEMP_VCF1" ]; then
    echo "Error: GQ filtered VCF file was not created."
    exit 1
fi

# Step 3: Filter variants by Depth (DP < 3)
bcftools filter -e 'FORMAT/DP < 3' -O z -o $TEMP_VCF2 $TEMP_VCF1

# Check if the DP filtered VCF file was created successfully
if [ ! -f "$TEMP_VCF2" ]; then
    echo "Error: DP filtered VCF file was not created."
    exit 1
fi

# Step 4: Filter variants by Missing Data (F_MISSING > 0.05)
bcftools filter -e 'F_MISSING > 0.05' -O z -o $OUTPUT_VCF $TEMP_VCF2

# Check if the final filtered VCF file was created successfully
if [ ! -f "$OUTPUT_VCF" ]; then
    echo "Error: Final filtered VCF file was not created."
    exit 1
fi

# Index the final VCF file
bcftools index $OUTPUT_VCF

# Final message
echo "Filtering complete. Variants stored in $OUTPUT_VCF."

# Clean up temporary files
rm -f $TEMP_VCF_INDEL_REMOVED $TEMP_VCF1 $TEMP_VCF2

