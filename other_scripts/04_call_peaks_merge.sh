#!/bin/bash
#SBATCH --job-name=call_peaks_merg
#SBATCH --output=callPeaks_%A_%a.out
#SBATCH --error=callPeaks_%A_%a.err
#SBATCH --array=2-5
#SBATCH --mem=64000
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=XXX

set -e  # Exit on error
set -u  # Exit if variable is undefined

eval "$(conda shell.bash hook)"
conda activate atac_env

# load in sample info (contained in a manifest file with the following structure:
# filename,
SAMPLE_LINE=$(sed -n ${SLURM_ARRAY_TASK_ID}p path/to/call_peaks_merged_manifest.csv)

FILENAME=""
GENOME="mm10"  # Change to "hg38" if using human data

while IFS=, read -r filename; do
    FILENAME="$filename"
done <<< "$SAMPLE_LINE"

BASE=$(basename "$FILENAME" .bam)

# Tn5 coordinate shift
echo "Shifting Tn5 coordinates for $BASE..."
Tn5_SHIFTED="${BASE}_shifted.bam"
alignmentSieve --bam "$FILENAME" --outFile "$Tn5_SHIFTED" --ATACshift --numberOfProcessors 4

# peak calling
echo "Calling peaks for $BASE..."
macs2 callpeak \
    -t "$Tn5_SHIFTED" \
    -f BAM \
    -g mm \
    -n "$BASE" \
    --outdir "macs2_out_${BASE}" \
    --nomodel --shift -100 --extsize 200 -q 0.01 --keep-dup 1

# blacklist removal
echo "Filtering peaks with ENCODE blacklist..."

# define blacklist
if [[ $GENOME == "hg38" ]]; then
    BLACKLIST="/path/to/GRCh38_blacklist.bed"
elif [[ $GENOME == "mm10" ]]; then
    BLACKLIST="/path/to/mm10_blacklist.bed"
else
    echo "Unsupported genome: $GENOME"
    exit 1
fi

# apply blacklist filtering
for peakfile in "macs2_out_${BASE}/${BASE}"*.narrowPeak "macs2_out_${BASE}/${BASE}"*.broadPeak; do
    if [[ -f "$peakfile" ]]; then
        OUTFILE="${peakfile%.*}_noBL.${peakfile##*.}"
        bedtools intersect -v -a "$peakfile" -b "$BLACKLIST" > "$OUTFILE"
        echo "Blacklist filtered: $OUTFILE"
    else
        echo "Skipping: $peakfile (file not found)"
    fi
done

echo "Finished processing $BASE!"
