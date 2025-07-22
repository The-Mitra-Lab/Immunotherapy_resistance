#!/bin/bash

# Input file (DESeq2 results with peak_id column)
INPUT_FILE="deseq2_results_down_significant.tsv"

# Intermediate and final BED filenames
RAW_BED="ifn_down_peaks.bed"
FINAL_BED="ifn_down_peaks_fixed.bed"

# Extract peak IDs and convert to BED format
cut -f1 "$INPUT_FILE" | grep -v peak_id > "$RAW_BED"
sed 's/_/\t/g' "$RAW_BED" > "$FINAL_BED"
rm "$RAW_BED"
