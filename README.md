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

