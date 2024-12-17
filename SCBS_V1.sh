#!/bin/bash

# Function to display usage information
usage() {
    echo "Usage: \$0 -f <input_folder> -o <output_folder>"
    exit 1
}

# Parse command-line options
while getopts ":f:o:" opt; do
    case $opt in
        f) INPUT_DIR="$OPTARG"
        ;;
        o) OUTPUT_DIR="$OPTARG"
        ;;
        \?) echo "Invalid option -$OPTARG" >&2
            usage
        ;;
        :) echo "Option -$OPTARG requires an argument." >&2
           usage
        ;;
    esac
done

# Check if both input and output directories are provided
if [ -z "$INPUT_DIR" ] || [ -z "$OUTPUT_DIR" ]; then
    usage
fi

# Create output directories for each step
mkdir -p ${OUTPUT_DIR}/FastQC ${OUTPUT_DIR}/Trimmed ${OUTPUT_DIR}/Aligned ${OUTPUT_DIR}/Deduplicated ${OUTPUT_DIR}/Methylation

# Step 1: Run FastQC to perform quality control
echo "Running FastQC..."
for fq in ${INPUT_DIR}/*.fq.gz; do
    fastqc $fq -o ${OUTPUT_DIR}/FastQC
done

# Step 2: Trim reads using Trim Galore! (including adapter trimming and quality trimming)
echo "Running Trim Galore!..."
for fq1 in ${INPUT_DIR}/*_1.fq.gz; do
    fq2="${fq1/_1.fq.gz/_2.fq.gz}"
    if [[ -f "$fq2" ]]; then
        trim_galore --fastqc --output_dir ${OUTPUT_DIR}/Trimmed $fq1 $fq2
    else
        trim_galore --fastqc --output_dir ${OUTPUT_DIR}/Trimmed $fq1
    fi
done

# Step 3: Align trimmed reads using Bismark (single-end, nondirectional)
echo "Running Bismark Alignment..."
for fq1 in ${OUTPUT_DIR}/Trimmed/*_1.fq.gz; do
    bismark --non_directional --genome /path/to/genome/ --output_dir ${OUTPUT_DIR}/Aligned $fq1
done

# Step 4: Remove duplicate reads
echo "Removing duplicates..."
for bam in ${OUTPUT_DIR}/Aligned/*.bam; do
    deduplicate_bismark --bam $bam --output_dir ${OUTPUT_DIR}/Deduplicated
done

# Step 5: Extract methylation information
echo "Extracting methylation calls..."
for dedup_bam in ${OUTPUT_DIR}/Deduplicated/*.bam; do
    bismark_methylation_extractor --gzip --bedGraph --output_dir ${OUTPUT_DIR}/Methylation $dedup_bam
done

echo "Pipeline completed successfully."
