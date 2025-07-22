#!/bin/bash

# input file with DESeq2 results
INPUT_FILE="deseq2_results_up_significant.tsv"

# output filenames
RAW_BED="up_peaks_raw.bed"
FINAL_BED="up_peaks.bed"

# extract peak_id column and remove header
cut -f1 "$INPUT_FILE" | grep -v peak_id > "$RAW_BED"

# convert underscores to tab-delimited BED format
sed 's/_/\t/g' "$RAW_BED" > "$FINAL_BED"

# remove intermediate file
rm "$RAW_BED"
