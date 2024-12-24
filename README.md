

# Single Cell DNA Methylation Analysis Pipeline

This repository contains a Bash script for performing a DNA methylation analysis pipeline. The pipeline includes quality control, read trimming, alignment, duplicate removal, and methylation extraction. It leverages tools like FastQC, Trim Galore!, Bismark, and ParaFly for parallel processing.

> Paper
https://www.nature.com/articles/nmeth.3035
> 
## Prerequisites

Before running the script, ensure you have the following tools installed:

- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [Trim Galore!](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)
- [Bismark](https://www.bioinformatics.babraham.ac.uk/projects/bismark/)
- [ParaFly](https://github.com/itmat/ParaFly)

Additionally, ensure you have access to a reference genome for Bismark alignment.

## Usage

To run the pipeline, use the following command:

```bash
 bash SCBS/SCBS_Parafly_V2.sh -f <input_folder> -o <output_folder> -t <threads>
```

### Options

- `-f <input_folder>`: Directory containing input FASTQ files (compressed with `.fq.gz`).
- `-o <output_folder>`: Directory where output files will be saved.
- `-t <threads>`: Number of threads to use for parallel processing (default is 1).

## Pipeline Steps

1. **Quality Control**: Runs FastQC on raw FASTQ files to assess the quality of the sequencing data.
2. **Read Trimming**: Uses Trim Galore! to trim adapters and low-quality bases from reads.
3. **Alignment**: Aligns trimmed reads to a reference genome using Bismark in single-end, non-directional mode.
4. **Duplicate Removal**: Removes duplicate reads from the aligned BAM files.
5. **Methylation Extraction**: Extracts methylation information from deduplicated BAM files.

## Output

The pipeline generates the following output directories within the specified output folder:

- `FastQC`: Contains FastQC reports for raw reads.
- `Trimmed`: Contains trimmed FASTQ files and FastQC reports post-trimming.
- `Aligned`: Contains BAM files of aligned reads.
- `Deduplicated`: Contains BAM files with duplicates removed.
- `Methylation`: Contains methylation extraction results.


## Installation

1. **Clone the repository:**

   ```bash
   git clone https://github.com/thekingofall/SCBS.git
   cd SCBS
   ```

2. **Create a new environment using Mamba:**

   ```bash
   mamba create -n SCBS -c bioconda fastqc trim-galore bismark
   ```

3. **Activate the environment:**

   ```bash
   conda activate SCBS
   ```
   
## Notes

- Ensure that the reference genome path is correctly specified in the script for Bismark alignment.
- The script assumes paired-end reads with filenames ending in `_1.fq.gz` and `_2.fq.gz`. Adjust the script if using single-end reads or different naming conventions.



