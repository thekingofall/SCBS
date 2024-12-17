# Methylation Analysis Pipeline

This repository provides a Bash script for a methylation analysis pipeline, which includes quality control, trimming, alignment, deduplication, and methylation extraction. The pipeline utilizes tools such as FastQC, Trim Galore!, Bismark, and others.

## Prerequisites

Before running the pipeline, ensure you have the following software installed:

- [Mamba](https://github.com/mamba-org/mamba): A fast, robust package manager for conda.

## Installation

1. **Clone the repository:**

   ```bash
   git clone https://github.com/yourusername/methylation-analysis-pipeline.git
   cd methylation-analysis-pipeline
   ```

2. **Create a new environment using Mamba:**

   ```bash
   mamba create -n methylation-pipeline -c bioconda fastqc trim-galore bismark
   ```

3. **Activate the environment:**

   ```bash
   conda activate methylation-pipeline
   ```

## Usage

To run the pipeline, use the following command:

```bash
bash run_pipeline.sh <input_folder> <output_folder>
```

- `<input_folder>`: Directory containing input FASTQ files (compressed with `.fq.gz`).
- `<output_folder>`: Directory where results will be saved.

### Example

```bash
bash run_pipeline.sh /path/to/input /path/to/output
```

## Notes

- Ensure that the path to the reference genome in the script is correctly set.
- The script assumes paired-end reads are named with `_1` and `_2` suffixes.



