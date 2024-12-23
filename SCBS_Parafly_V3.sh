#!/bin/bash

# Function to display usage information
usage() {
    echo "Usage: \\$0 -f <input_folder> -o <output_folder> -t <threads> [-g <genome_path>]"
    exit 1
}

# Default values
THREADS=1
GENOME_PATH="/home/maolp/mao/Ref/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/HG38fa/"

# Parse command-line options
while getopts ":f:o:t:g:" opt; do
    case $opt in
        f) INPUT_DIR="$(pwd)/$OPTARG"
        ;;
        o) OUTPUT_DIR="$(pwd)/$OPTARG"
        ;;
        t) THREADS="$OPTARG"
        ;;
        g) GENOME_PATH="$(pwd)/$OPTARG"
        ;;
        \?) echo "Invalid option -$OPTARG" >&2
            usage
        ;;
        :) echo "Option -$OPTARG requires an argument." >&2
           usage
        ;;
    esac
done

# Check if input and output directories are provided
if [ -z "$INPUT_DIR" ] || [ -z "$OUTPUT_DIR" ]; then
    usage
fi

# Create output directories for each step
mkdir -p "${OUTPUT_DIR}/FastQC" "${OUTPUT_DIR}/Trimmed" "${OUTPUT_DIR}/Aligned" "${OUTPUT_DIR}/Deduplicated" "${OUTPUT_DIR}/Methylation" "${OUTPUT_DIR}/shell"

# Create a directory for shell scripts
SHELL_DIR="${OUTPUT_DIR}/shell"
mkdir -p "$SHELL_DIR"

# Step 1: Run FastQC to perform quality control
echo "Running FastQC..."
FASTQC_CMD_FILE="${SHELL_DIR}/fastqc_commands.txt"
> "$FASTQC_CMD_FILE" # Create an empty file, overwrite if exists
for fq in "${INPUT_DIR}"/*.fq.gz; do
    echo "fastqc \"$fq\" -o \"${OUTPUT_DIR}/FastQC\"" >> "$FASTQC_CMD_FILE"
done
ParaFly -c "$FASTQC_CMD_FILE" -CPU "$THREADS" -v

# Step 2: Trim reads using Trim Galore! (including adapter trimming and quality trimming)
echo "Running Trim Galore!..."
TRIM_CMD_FILE="${SHELL_DIR}/trim_commands.txt"
> "$TRIM_CMD_FILE"
for fq1 in "${INPUT_DIR}"/*_1.fq.gz; do
    fq2="${fq1/_1.fq.gz/_2.fq.gz}"
    
    # Check if the paired-end file exists
    if [[ -f "$fq2" ]]; then
        echo "trim_galore --paired --fastqc --output_dir \"${OUTPUT_DIR}/Trimmed\" \"$fq1\" \"$fq2\"" >> "$TRIM_CMD_FILE"
    else
        echo "echo \"Warning: Paired file for $fq1 not found. Skipping.\" >&2" >> "$TRIM_CMD_FILE"
    fi
done
ParaFly -c "$TRIM_CMD_FILE" -CPU "$THREADS" -v

# Step 3: Align trimmed reads using Bismark (paired-end, nondirectional)
echo "Running Bismark Alignment..."
BISMARK_CMD_FILE="${SHELL_DIR}/bismark_commands2.txt"
> "$BISMARK_CMD_FILE"
for fq1 in "${OUTPUT_DIR}/Trimmed/"*_val_1.fq.gz; do
    # 修改替换模式，确保生成正确的配对文件路径
    fq2="${fq1/_1_val_1.fq.gz/_2_val_2.fq.gz}"
    
    # 检查配对文件是否存在
    if [[ -f "$fq2" ]]; then
        echo "bismark --non_directional --genome \"$GENOME_PATH\" --output_dir \"${OUTPUT_DIR}/Aligned\" -1 \"$fq1\" -2 \"$fq2\"" >> "$BISMARK_CMD_FILE"
    else
        echo "echo \"Warning: Trimmed paired file for $fq1 not found. Skipping alignment.\" >&2" >> "$BISMARK_CMD_FILE"
    fi
done
ParaFly -c "$BISMARK_CMD_FILE" -CPU "$THREADS" -v

# Step 4: Remove duplicate reads
echo "Removing duplicates..."
DEDUP_CMD_FILE="${SHELL_DIR}/dedup_commands2.txt"
> "$DEDUP_CMD_FILE"
for bam in "${OUTPUT_DIR}/Aligned/"*.bam; do
    echo "deduplicate_bismark --bam \"$bam\" --paired --output_dir \"${OUTPUT_DIR}/Deduplicated\"" >> "$DEDUP_CMD_FILE"
done
ParaFly -c "$DEDUP_CMD_FILE" -CPU "$THREADS" -v

# Step 5: Extract methylation information
echo "Extracting methylation calls..."
METHYL_CMD_FILE="${SHELL_DIR}/methyl_commands2.txt"
> "$METHYL_CMD_FILE"
for dedup_bam in "${OUTPUT_DIR}/Deduplicated/"*.bam; do
    echo "bismark_methylation_extractor --gzip --bedGraph --paired-end --output_dir \"${OUTPUT_DIR}/Methylation\" \"$dedup_bam\"" >> "$METHYL_CMD_FILE"
done
ParaFly -c "$METHYL_CMD_FILE" -CPU "$THREADS" -v

# Step 6: Generate CX report
echo "Generating CX reports..."
CX_REPORT_CMD_FILE="${SHELL_DIR}/cx_report_commands.txt"
> "$CX_REPORT_CMD_FILE"
for dedup_bam in "${OUTPUT_DIR}/Deduplicated/"*.bam; do
    sample_name=$(basename "$dedup_bam" .bam)
    output_dir="${OUTPUT_DIR}/Methylation/CX_Reports"
    mkdir -p "$output_dir"
    
    echo "bismark_methylation_extractor \
        -p \
        --comprehensive \
        --no_overlap \
        --bedGraph \
        --cytosine_report \
        --CX_context \
        #--split_by_chromosome \
        --counts \
        --buffer_size 20G \
        --genome_folder \"$GENOME_PATH\" \
        \"$dedup_bam\" \
        --output_dir \"$output_dir\"" >> "$CX_REPORT_CMD_FILE"
done
ParaFly -c "$CX_REPORT_CMD_FILE" -CPU "$THREADS" -v

echo "Pipeline completed successfully."
