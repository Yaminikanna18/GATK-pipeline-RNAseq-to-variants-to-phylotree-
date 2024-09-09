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

