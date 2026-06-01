#!/usr/bin/env bash
set -euo pipefail

REF_DIR="GRCh38_ZIKV"        # output dir from mkref
FASTQ_DIR="/path/to/fastqs"  # demultiplexed 10x FASTQs
THREADS=8
MEM_GB=64

# One entry per capture; names must match the FASTQ sample prefixes
SAMPLES=(cap1 cap2 cap3 cap4 cap5 cap6 cap7 cap8)

for s in "${SAMPLES[@]}"; do
  cellranger count \
    --id="${s}_count" \
    --transcriptome="$REF_DIR" \
    --fastqs="$FASTQ_DIR" \
    --sample="$s" \
    --chemistry=auto \
    --create-bam=true \
    --localcores="$THREADS" \
    --localmem="$MEM_GB"
done