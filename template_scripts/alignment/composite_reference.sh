#!/usr/bin/env bash
set -euo pipefail

# Inputs (set to matched-release host files and the viral genome FASTA)
HOST_FASTA="GRCh38.primary_assembly.genome.fa"
HOST_GTF="gencode.primary_assembly.annotation.gtf"
VIRUS_FASTA="ZIKV.fa"
VIRUS_GENE="ZIKV"
REF_NAME="GRCh38_ZIKV"
THREADS=8
MEM_GB=64

# Pull the viral contig name (col 1 of the GTF must match the FASTA header)
VIRUS_CONTIG=$(grep '^>' "$VIRUS_FASTA" | head -n1 | sed 's/^>//; s/ .*//')

# Total length of the viral record (assumes a single sequence in the FASTA)
VIRUS_LEN=$(awk '/^>/{next} {n+=length($0)} END{print n}' "$VIRUS_FASTA")

# Annotated feature span: defaults to the full genome.
# For primer-based flavivirus capture, restrict to the region upstream of the
# virus-specific primer where reads accumulate (e.g. 1-996 for PR 996R).
FEAT_START=1
FEAT_END="$VIRUS_LEN"

# Write a minimal viral GTF (tab-delimited) with one gene/transcript/exon
printf '%s\tcustom\tgene\t%s\t%s\t.\t+\t.\tgene_id "%s"; gene_name "%s"; gene_biotype "protein_coding";\n' \
  "$VIRUS_CONTIG" "$FEAT_START" "$FEAT_END" "$VIRUS_GENE" "$VIRUS_GENE" > virus.gtf
printf '%s\tcustom\ttranscript\t%s\t%s\t.\t+\t.\tgene_id "%s"; transcript_id "%s_tx"; gene_name "%s"; gene_biotype "protein_coding";\n' \
  "$VIRUS_CONTIG" "$FEAT_START" "$FEAT_END" "$VIRUS_GENE" "$VIRUS_GENE" "$VIRUS_GENE" >> virus.gtf
printf '%s\tcustom\texon\t%s\t%s\t.\t+\t.\tgene_id "%s"; transcript_id "%s_tx"; gene_name "%s"; gene_biotype "protein_coding";\n' \
  "$VIRUS_CONTIG" "$FEAT_START" "$FEAT_END" "$VIRUS_GENE" "$VIRUS_GENE" "$VIRUS_GENE" >> virus.gtf

# Filter the host GTF to standard biotypes
cellranger mkgtf "$HOST_GTF" host_filtered.gtf \
  --attribute=gene_biotype:protein_coding \
  --attribute=gene_biotype:lncRNA \
  --attribute=gene_biotype:IG_C_gene \
  --attribute=gene_biotype:IG_D_gene \
  --attribute=gene_biotype:IG_J_gene \
  --attribute=gene_biotype:IG_V_gene \
  --attribute=gene_biotype:TR_C_gene \
  --attribute=gene_biotype:TR_D_gene \
  --attribute=gene_biotype:TR_J_gene \
  --attribute=gene_biotype:TR_V_gene

# Concatenate host and virus into one reference
cat "$HOST_FASTA" "$VIRUS_FASTA" > combined.fa
cat host_filtered.gtf virus.gtf      > combined.gtf

# Build the Cell Ranger reference
cellranger mkref \
  --genome="$REF_NAME" \
  --fasta=combined.fa \
  --genes=combined.gtf \
  --nthreads="$THREADS" \
  --memgb="$MEM_GB"