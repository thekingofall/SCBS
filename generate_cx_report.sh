#!/bin/bash

# Function to display usage information
usage() {
    echo "Usage: \$0 -i <input_bam_dir> -o <output_dir> [-g <genome_path>] [-t <threads>]"
    exit 1
}

# Default values
THREADS=1
GENOME_PATH="/home/maolp/mao/Ref/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/HG38fa/"

# Parse command-line options
while getopts ":i:o:g:t:" opt; do
    case $opt in
        i) INPUT_DIR="$OPTARG"
        ;;
        o) OUTPUT_DIR="$OPTARG"
        ;;
        g) GENOME_PATH="$OPTARG"
        ;;
        t) THREADS="$OPTARG"
        ;;
        \?) echo "Invalid option -$OPTARG" >&2
            usage
        ;;
        :) echo "Option -$OPTARG requires an argument." >&2
           usage
        ;;
    esac
done

# Check if required parameters are provided
if [ -z "$INPUT_DIR" ] || [ -z "$OUTPUT_DIR" ]; then
    usage
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Create command file
CMD_FILE="${OUTPUT_DIR}/cx_report_commands.txt"
> "$CMD_FILE"

# Generate commands for each BAM file
for bam in "${INPUT_DIR}"/*.bam; do
    echo "bismark_methylation_extractor \
        -p \
        --comprehensive \
        --no_overlap \
        --bedGraph \
        --cytosine_report \
        --CX_context \
        --split_by_chromosome \
        --counts \
        --buffer_size 20G \
        --genome_folder \"$GENOME_PATH\" \
        \"$bam\" \
        --output_dir \"$OUTPUT_DIR\"" >> "$CMD_FILE"
done

# Execute commands in parallel
ParaFly -c "$CMD_FILE" -CPU "$THREADS" -v

echo "CX report generation completed."
