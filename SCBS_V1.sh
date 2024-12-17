#!/bin/bash

if [ $# -ne 2 ]; then
    echo "Usage: \$0 <input_folder> <output_folder>"
    exit 1
fi

INPUT_DIR="\$1"  # 使用双引号
OUTPUT_DIR="\$2" # 使用双引号

# Create output directories for each step
mkdir -p "${OUTPUT_DIR}/FastQC" "${OUTPUT_DIR}/Trimmed" "${OUTPUT_DIR}/Aligned" "${OUTPUT_DIR}/Deduplicated" "${OUTPUT_DIR}/Methylation"

# Step 1: Run FastQC to perform quality control
echo "Running FastQC..."
for fq in "${INPUT_DIR}"/*.fq.gz; do
    fastqc "$fq" -o "${OUTPUT_DIR}/FastQC"
done

# Step 2: Trim reads using Trim Galore! (including adapter trimming and quality trimming)
echo "Running Trim Galore!..."
for fq1 in "${INPUT_DIR}"/*_1.fq.gz; do
    fq2="${fq1/_1.fq.gz/_2.fq.gz}"
    if [[ -f "$fq2" ]]; then
        trim_galore --fastqc --output_dir "${OUTPUT_DIR}/Trimmed" "$fq1" "$fq2"
    else
        trim_galore --fastqc --output_dir "${OUTPUT_DIR}/Trimmed" "$fq1"
    fi
done

# Step 3: Align trimmed reads using Bismark (single-end, nondirectional)
echo "Running Bismark Alignment..."
for fq1 in "${OUTPUT_DIR}/Trimmed"/*_1.fq.gz; do
    bismark --non_directional --genome /path/to/genome/ --output_dir "${OUTPUT_DIR}/Aligned" "$fq1"
done

# Step 4: Remove duplicate reads
echo "Removing duplicates..."
for bam in "${OUTPUT_DIR}/Aligned"/*.bam; do
    deduplicate_bismark --bam "$bam" --output_dir "${OUTPUT_DIR}/Deduplicated"
done

# Step 5: Extract methylation information
echo "Extracting methylation calls..."
for dedup_bam in "${OUTPUT_DIR}/Deduplicated"/*.bam; do
    bismark_methylation_extractor --gzip --bedGraph --output_dir "${OUTPUT_DIR}/Methylation" "$dedup_bam"
done

echo "Pipeline completed successfully."
