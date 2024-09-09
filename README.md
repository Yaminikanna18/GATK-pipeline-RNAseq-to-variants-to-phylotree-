## SRA Download Script

This script downloads SRA files and converts them to FASTQ format using SRA Toolkit. Before using the script, replace the placeholders with your specific values.

### Replaceable Placeholders:

- `<job_name>`: Your SLURM job name.
- `<output_directory>`: The directory where you want to save output files.
- `<time_limit>`: The time limit for the job (e.g., `24:00:00`).
- `<memory_limit>`: The memory required (e.g., `32G`).
- `<account_name>`: Your SLURM account name.
- `<sra_toolkit_version>`: The version of SRA Toolkit you are using.
- `<SRR_ACC1>`, `<SRR_ACC2>`, etc.: Your SRR accession numbers.

### Running the Script

To submit the job to SLURM, use the following command:

```bash
sbatch sra_download_script.sh

# FastQC and multiqc Analysis Pipeline

This repository contains a SLURM script for running FastQC on FASTQ files.

## Script Details

- **fastqc_script.sh**: A SLURM script to run FastQC on FASTQ files. The script performs quality checks on files with suffixes `_1.f$
## Usage

1. **Update SLURM Details**: Before using the script, make sure to modify the placeholders in the `fastqc_script.sh` script to inclu$
2. **Modify Directories**: Replace the placeholders in the `fastqc_script.sh` script with the appropriate directory paths for your F$
3. **Submit the Script**: Submit the script to the SLURM scheduler using the command:

   ```bash
   sbatch fastqc_script.sh
   
   conda activate <your_conda_environment> # Running mulitiqc is optional
   multiqc /path/to/fastqc_results

# Trimmomatic Analysis Pipeline

This repository contains a SLURM script for running Trimmomatic on paired-end FASTQ files. 

## Script Details

- **trimmomatic_script.sh**: A SLURM script to run Trimmomatic on FASTQ files. The script performs trimming on files with suffixes `_1.fastq.gz` and `_2.fastq.gz`.

## Usage

1. Modify the placeholders in the `trimmomatic_script.sh` script to point to the appropriate directories and adjust SLURM parameters:
   - **<job_name>**: Job name.
   - **<output_directory>**: Directory for output and error files.
   - **<number_of_nodes>**: Number of nodes.
   - **<tasks_per_node>**: Number of tasks per node.
   - **<memory>**: Memory per node.
   - **<time_limit>**: Time limit.
   - **<account_name>**: Account name.
   - **<java_version>**: Java module version.
   - **<trimmomatic_version>**: Trimmomatic module version.
   - **<input_directory>**: Directory containing input FASTQ files.
   - **<output_directory>**: Directory to save the trimmed FASTQ files.

2. Submit the script to the SLURM scheduler using the command:

   ```bash
   sbatch trimmomatic_script.sh

## FastQC Analysis on Trimmed Files

- **fastqc_after_trimming_script.sh**: A script to run FastQC on trimmed FASTQ files. The script performs quality checks on files with suffixes `_1_trimmed_paired.fastq.gz` and `_2_trimmed_paired.fastq.gz`.

### Usage

1. **Modify Directories**: Replace the placeholders in the `fastqc_after_trimming_script.sh` script with the appropriate directory paths for your trimmed FASTQ files and desired output location.

2. **Run the Script**: After updating the script, you can submit it to your job scheduler or run it locally, depending on your setup.

   ```bash
   bash fastqc_after_trimming_script.sh


   conda activate <your_conda_environment>
   multiqc /path/to/trimmed_fastqc_results

In this pipeline, we concatenated the FASTQ files for different tissues of the same individual using the cat command. For example, forward and reverse reads for each sample across tissues (Brain, Gill, Liver, etc.) were merged into a single file for each individual.

cat <input_dir>/<sample>_Brain_1_trimmed_paired.fastq.gz \
    <input_dir>/<sample>_Gill_1_trimmed_paired.fastq.gz \
    <input_dir>/<sample>_Liver_1_trimmed_paired.fastq.gz \
> <output_dir>/<sample>_All_1_trimmed_paired.fastq.gz



Once the concatenated FASTQ files are ready, download the reference genome and its corresponding GTF file from the NCBI repository. These files are necessary for downstream analyses such as alignment or quantification.


# Picard FastqToSam conversion 

This repository contains a SLURM script for converting FASTQ files to unmapped BAM files using Picard's `FastqToSam` tool.

## Script Details

- **picard_fastqtosam.sh**: A SLURM script to run Picard's `FastqToSam` command on paired-end FASTQ files. This script converts FASTQ files into unmapped BAM files, which are useful for further processing and alignment.

## Usage

1. **Update SLURM Details**: Modify the placeholders in the `picard_fastqtosam.sh` script to fit your specific requirements:
   - **`<account_name>`**: Your SLURM account name.
   - **`<your_conda_environment>`**: The name of your Conda environment where Picard is installed.
   - **`/path/to/your/input_directory`**: Directory containing input FASTQ files.
   - **`/path/to/your/output_directory`**: Directory where you want to save the unmapped BAM files.

2. **Define Base Filenames**: Update the `base_filenames` array with the base names of your samples (excluding `_1` and `_2` suffixes). For example:
   - `"sample1"`
   - `"sample2"`

3. **Submit the Script**: After updating the script, use the following command to submit it to the SLURM scheduler:

    ```bash
    sbatch picard_fastqtosam.sh
    ```

Make sure to replace all placeholders with your actual values before running the script.


## STAR Genome Indexing Script

This script generates a STAR genome index from a genome FASTA file and an annotation GTF file.

### Script Details

- **star_indexing.sh**: A SLURM script to create a STAR genome index.

### Replaceable Placeholders:

- `<account_name>`: Your SLURM account name.
- `<your_conda_environment>`: The name of your Conda environment where STAR is installed.
- `<working_directory>`: The directory where you will store the STAR index and where your genome files are located.
- `<genome_fasta_file>`: The path to your genome FASTA file.
- `<gtf_file>`: The path to your genome GTF file.

### Running the Script

To submit the job to SLURM, use the following command:

```bash
sbatch star_indexing.sh


## STAR Alignment Script

This script performs alignment of paired-end FASTQ files using STAR aligner. Before using the script, replace the placeholders with your specific values.

### Script Details

- **star_align.sh**: A SLURM script to perform alignment of paired-end FASTQ files using STAR.

### Replaceable Placeholders:

- `<job_name>`: Your SLURM job name.
- `<num_cpus>`: Number of CPU cores per task.
- `<time_limit>`: Time limit for the job (e.g., `46:00:00`).
- `<memory_limit>`: Memory required (e.g., `68G`).
- `<account_name>`: Your SLURM account name.
- `<path_to_conda.sh>`: Path to your Conda initialization script.
- `<your_conda_environment>`: Name of your Conda environment.
- `<working_directory>`: Directory for working files.
- `<index_directory>`: Directory containing STAR index files.

### Running the Script

To submit the job to SLURM, use the following command:

```bash
sbatch star_align.sh <read1> <read2>


## BAM File Merging Script

This script merges unmapped and aligned BAM files using Picard's `MergeBamAlignment` tool.

### Script Details

- **merge_bam_alignment.sh**: A SLURM script to merge unmapped and aligned BAM files into a single BAM file.

### Replaceable Placeholders:

- `<job_name>`: Your SLURM job name.
- `<output_directory>`: Directory where you want to save output and error files.
- `<cpus_per_task>`: Number of CPU cores per task (e.g., `4`).
- `<memory_limit>`: Memory limit per node (e.g., `68G`).
- `<time_limit>`: Time limit for the job (e.g., `24:00:00`).
- `<account_name>`: Your SLURM account name.
- `<conda_environment>`: Name of the Conda environment to activate.
- `<reference_fasta>`: Path to the reference FASTA file.

### Running the Script

To submit the job to SLURM, use the following command:

```bash
sbatch merge_bam_alignment.sh <unmapped_bam> <aligned_bam> <output_bam>


## MarkDuplicates Script

This script marks duplicate reads in BAM files using Picard's `MarkDuplicates` tool.

### Script Details

- **mark_duplicates.sh**: A SLURM script to mark duplicate reads in BAM files.

### Replaceable Placeholders:

- `<output_directory>`: Directory where you want to save output and error files.
- `<memory_limit>`: Memory limit per node (e.g., `32G`).
- `<time_limit>`: Time limit for the job (e.g., `10:00:00`).
- `<account_name>`: Your SLURM account name.
- `<conda_environment>`: Name of the Conda environment to activate.

### Running the Script

To submit the job to SLURM, use the following command:

```bash
sbatch mark_duplicates.sh <input_bam> <output_bam> <metrics_file>

## SplitNCigarReads Script

This script uses GATK's `SplitNCigarReads` to split reads with N's in their CIGAR string.

### Script Details

- **split_ncigars.sh**: A SLURM script to run GATK's `SplitNCigarReads` on BAM files.

### Replaceable Placeholders:

- `<output_directory>`: Directory where you want to save output and error files.
- `<memory_limit>`: Memory limit per node (e.g., `32G`).
- `<time_limit>`: Time limit for the job (e.g., `10:00:00`).
- `<account_name>`: Your SLURM account name.
- `<gatk_version>`: The version of GATK you are using.

### Running the Script

To submit the job to SLURM, use the following command:

```bash
sbatch split_ncigars.sh <input_bam> <output_bam> <reference>

## HaplotypeCaller Script

This script uses GATK's `HaplotypeCaller` to call variants and generate a GVCF file.

### Script Details

- **haplotypecaller.sh**: A SLURM script to run GATK's `HaplotypeCaller` on BAM files.

### Replaceable Placeholders:

- `<output_directory>`: Directory where you want to save output and error files.
- `<time_limit>`: Time limit for the job (e.g., `2:00:00`).
- `<memory_limit>`: Memory limit per node (e.g., `32G`).
- `<account_name>`: Your SLURM account name.
- `<gatk_version>`: The version of GATK you are using.
- `<reference_genome>`: Path to the reference genome file.
- `<input_bam>`: Path to the input BAM file.
- `<output_gvcf>`: Path to save the output GVCF file.

### Running the Script

To submit the job to SLURM, use the following command:

```bash
sbatch haplotypecaller.sh <input_bam> <output_gvcf> <reference_genome>


## Combine GVCFs Script

This script combines multiple GVCF files into a single GVCF file using GATK's `CombineGVCFs`.

### Script Details

- **combine_gvcfs.sh**: A SLURM script to combine GVCF files into a single GVCF file using GATK's `CombineGVCFs`.

### Replaceable Placeholders:

- `<output_directory>`: Directory where you want to save output and error files.
- `<time_limit>`: Time limit for the job (e.g., `46:00:00`).
- `<memory_limit>`: Memory limit per node (e.g., `32G`).
- `<account_name>`: Your SLURM account name.
- `<gatk_version>`: The version of GATK you are using.
- `<gvcf_directory>`: Directory containing GVCF files to be combined.
- `<reference_directory>`: Directory containing the reference genome file.
- `<reference_genome>`: Name of the reference genome file.

### Running the Script

To submit the job to SLURM, use the following command:

```bash
sbatch combine_gvcfs.sh


## Joint Genotyping Script

This script performs joint genotyping on combined GVCF files to generate raw variant calls using GATK's `GenotypeGVCFs`.

### Script Details

- **jointgenotyping.sh**: A SLURM script to perform joint genotyping on GVCF files using GATK's `GenotypeGVCFs`.

### Replaceable Placeholders:

- `<output_directory>`: Directory where you want to save output and error files.
- `<time_limit>`: Time limit for the job (e.g., `46:00:00`).
- `<memory_limit>`: Memory limit per node (e.g., `32G`).
- `<account_name>`: Your SLURM account name.
- `<gatk_version>`: The version of GATK you are using.
- `<reference_directory>`: Directory containing the reference genome file.
- `<reference_genome>`: Name of the reference genome file.

### Running the Script

To submit the job to SLURM, use the following command:

```bash
sbatch jointgenotyping.sh

## Variant Filtration Script

This script filters VCF files to remove indels, apply genotype quality and depth filters, and manage missing data using `bcftools`.

### Script Details

- **variant_filtration.sh**: A SLURM script to filter variants from a VCF file using `bcftools`.

### Replaceable Placeholders:

- `<output_directory>`: Directory where you want to save output and error files.
- `<memory_limit>`: Memory limit per node (e.g., `32G`).
- `<time_limit>`: Time limit for the job (e.g., `24:00:00`).
- `<account_name>`: Your SLURM account name.
- `<bcftools_version>`: The version of `bcftools` you are using.
- `<input_vcf>`: Path to the input VCF file.
- `<temp_vcf_indel_removed>`: Path for the temporary VCF file with indels removed.
- `<temp_vcf_gq_filtered>`: Path for the temporary VCF file with genotype quality filtered.
- `<temp_vcf_dp_filtered>`: Path for the temporary VCF file with depth filtered.
- `<output_vcf>`: Path to the final filtered VCF file.

### Running the Script

To submit the job to SLURM, use the following command:

```bash
sbatch variant_filtration.sh


## VCF to Table Conversion Script

This script converts a VCF file into a tabular format using GATK.

### Script Details

- **vcf_to_tab.sh**: A SLURM script to convert a VCF file into a tabular format using GATK's `VariantsToTable` command.

### Replaceable Placeholders:

- `<output_directory>`: Directory where you want to save output and error files.
- `<memory_limit>`: Memory limit per node (e.g., `32G`).
- `<time_limit>`: Time limit for the job (e.g., `2:00:00`).
- `<account_name>`: Your SLURM account name.
- `<gatk_version>`: The version of GATK you are using.
- `<input_vcf>`: Path to the input VCF file.
- `<output_table>`: Path for the output tabular file.

### Running the Script

To submit the job to SLURM, use the following command:

```bash
sbatch vcf_to_tab.sh

Once you have converted the VCF file to a tabular format, you will need to process this table to generate a FASTA file. This involves several steps. First, you should clean the table by removing any unwanted symbols or special characters that may be present. Next, convert any ambiguous nucleotide representations to IUPAC codes, which are used to denote multiple possible nucleotides at a given position. Finally, use the cleaned table with IUPAC codes to generate a FASTA file. Ensure that you have the appropriate tools or scripts to perform these conversions and verify that your final FASTA file meets your requirements for subsequent analyses.

## MAFFT Alignment Script

This script performs multiple sequence alignment using the MAFFT tool. The alignment is performed on a FASTA file containing variant sequences, preparing the data for subsequent phylogenetic analysis.

### Script Details

- **mafft_alignment.sh**: A SLURM script to align sequences in a FASTA file using MAFFT.

### Replaceable Placeholders:

- `<job_name>`: Your SLURM job name.
- `<output_directory>`: The directory where you want to save output and error files.
- `<memory_limit>`: Memory limit per node (e.g., `32G`).
- `<time_limit>`: Time limit for the job (e.g., `2:00`).
- `<account_name>`: Your SLURM account name.
- `<mafft_version>`: The version of MAFFT you are using.
- `<input_fasta>`: Path to the input FASTA file.
- `<output_alignment>`: Path for the output aligned FASTA file.

### Running the Script

To submit the job to SLURM, use the following command:

```bash
sbatch mafft_alignment.sh


## RAxML Phylogenetic Tree Construction

This script performs phylogenetic tree construction using RAxML.

### Script Details

- **raxml_analysis.sh**: A SLURM script to run RAxML for constructing a phylogenetic tree.

### Replaceable Placeholders:

- `<job_name>`: Name of your SLURM job.
- `<output_directory>`: Directory where you want to save output files.
- `<account_name>`: Your SLURM account name.
- `<openmpi_version>`: The version of OpenMPI you are using.
- `<raxml_version>`: The version of RAxML you are using.
- `<input_fasta>`: The input FASTA file.
- `<output_tree_name>`: Name of the output tree file.
- `<random_seed>`: Random seed for RAxML.
- `<bootstrap_replicates>`: Number of bootstrap replicates.

### Running the Script

To submit the job to SLURM, use the following command:

```bash
sbatch raxml_analysis.sh

