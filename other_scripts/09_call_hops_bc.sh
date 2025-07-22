#!/bin/bash
#
#SBATCH --job-name=bulkRNACallingCardsBarcodes
#SBATCH --output=logs/bulkRNACallingCardsBarcodes_%a.out
#SBATCH --error=logs/bulkRNACallingCardsBarcodes_%a.err
#SBATCH --array=2-7 # for parallel processing
#SBATCH --cpus-per-task=2
#SBATCH --mem=64000
#SBATCH --mail-type=ALL
#SBATCH --mail-user=xxx

#Assign these variables
GENOME=mm10 #mm10 or hg38 only
#Change the name of your experiment directory (the directory that contains scripts, raw, and output_and_analysis)
PROJECT_DIR=/path/to/project
BOWTIE2_INDEX_PATH_AND_PREFIX=/path/to/bowtie/mm10
TWO_BIT_PATH=/path/to/mm10.2bit

# load required packages with spack
eval $(spack load --sh py-cutadapt@2.9 )
eval $(spack load --sh py-umi-tools@1.0.0 )
eval $(spack load --sh bowtie2@2.3.5 )
eval $(spack load --sh samtools@1.13 )
eval $(spack load --sh py-numpy@1.21.3 /i7mcgz4 )
eval $(spack load --sh py-pandas@1.3.4 /erg2vc5 )
eval $(spack load --sh py-pysam@0.18.0 )
eval $(spack load --sh py-twobitreader@3.1.7 )
eval $(spack load --sh py-importlib-metadata )

# confirm software verions are correct
echo -e "\n\ncutadapt version:"
cutadapt --version

echo -e "\n\nUMI-tools version:"
umi_tools --version

echo -e "\n\nbowtie2 version"
bowtie2 --version

echo -e "\n\nsamtools version:"
samtools --version

echo -e "\n\nnumpy version:"
python3 -c "import numpy; print(numpy.version.version)"

echo -e "\n\npandas version:"
python3 -c "import pandas; print(pandas.__version__)"

echo -e "\n\npysam version:"
python3 -c "import pysam; print(pysam.__version__)"

#echo -e "\n\ntwobitreader version"
#python3 -c "from importlib_metadata import version; print(version('twobitreader'))"

# Initialize settings
PROJECT_OUT=$PROJECT_DIR"/output_and_analysis"
PROJECT_RAW=$PROJECT_DIR"/raw"

# read in sample manifest file for parallel processing
SAMPLE_LINE=$( sed -n ${SLURM_ARRAY_TASK_ID}p call_hops_manifest.csv )
FILENAME=
BARCODE=
INDEX1=
INDEX2=
while IFS=, read -r filename barcode index1 index2
do
    { FILENAME=$filename; }
    { BARCODE=$barcode; }
    { INDEX1=$index1; }
    { INDEX2=$index2; }
done <<< "$SAMPLE_LINE"

#Set up stem for naming output files.
BASE=`basename $FILENAME`
STEM=$BARCODE"_""${BASE%%.*}"

#Input list of SRT barcodes
safelist_full=$'srt_barcode_whitelist.txt'

# Prepare output variables
OUT_STEM=$PROJECT_OUT/$STEM
OUT_SAM=$OUT_STEM".sam"
OUT_BAM=$OUT_STEM".bam"
OUT_MAP_SORT_PREFIX=$OUT_STEM"_map_sort"
OUT_BAM_MAP=$OUT_MAP_SORT_PREFIX".bam"
OUT_BED=$OUT_MAP_SORT_PREFIX".bed"
OUT_BEDGRAPH=$OUT_MAP_SORT_PREFIX".bedgraph"

# Transposon TR sequence
TRANSPOSASE=PB
PB_LTR_1="CGTCAATTTTACGCAGACTATCTTT"
PB_LTR_2="GTTAA"

echo -e "\n\nSUMMARY FOR 1st ADAPTER TRIMMING"
#Extract reads with the OM-PB-3bp combo
cutadapt \
    -g "^"$BARCODE$PB_LTR_1 \
    -o $OUT_STEM".temp1.fastq.gz" \
    --minimum-length 1 \
    --discard-untrimmed \
    -e 0 \
    --no-indels \
    $PROJECT_RAW/$FILENAME


#Extract reads with the srt barcode and add srt barcode (and fake umi) to read id
echo -e "\n\nSUMMARY FOR EXTRACTING SRT BARCODES"
umi_tools extract \
    --stdin $OUT_STEM".temp1.fastq.gz" \
    --stdout $OUT_STEM".temp2.fastq.gz" \
    --extract-method=string \
    --bc-pattern=CCCCN \
    --whitelist=$safelist_full \
    --filter-cell-barcode;


echo -e "\n\nSUMMARY FOR 2nd ADAPTER TRIMMING"
cutadapt \
    -g "^"$PB_LTR_2 \
    -o $OUT_STEM".temp3.fastq.gz" \
    --minimum-length 1 \
    --discard-untrimmed \
    -e 0 \
    --no-indels \
    $OUT_STEM".temp2.fastq.gz"

#add barcode to read id
#while read bar; do
#   addon=_G_"$bar"_
#   sed -i "N;/\n.*$bar/s/ 1:N:0:\(.*[^\n]\)/$addon 1:N:0:\1/" $OUT_BC
#done < $safelist_full

# Trim and demultiplex by TR and insertion site (exact match)
#cutadapt \
#    -g "^"$BARCODE$Long_PB_LTR \
#    -o $OUT_STEM".temp.fastq.gz" \
#z    --minimum-length 1 \
#    --discard-untrimmed \
#    -e 0 \
#    --no-indels \
 #   $OUT_BC".gz"

# Trim any trailing Nextera adapters (allowing mismatches)
echo -e "\n\nSUMMARY FOR 3rd ADAPTER TRIMMING (percent trimmed not important)"
cutadapt \
    -a "CTGTCTCTTATACACATCTCCGAGCCCACGAGACTNNNNNNNNNNTCTCGTATGCCGTCTTCTGCTTG" \
    -o $OUT_STEM".fastq.gz" \
    --minimum-length 1 \
    $OUT_STEM".temp3.fastq.gz"

# Align the reads
bowtie2 \
    -x $BOWTIE2_INDEX_PATH_AND_PREFIX \
    -U $OUT_STEM".fastq.gz" \
    -S $OUT_SAM

# Convert to BAM
samtools view \
    -bS -h \
    $OUT_SAM \
    -o $OUT_BAM

# Filter only mapped reads, convert to BAM, and sort
samtools view \
    -bS -h -F 260 \
    $OUT_SAM | \
samtools sort - -o $OUT_BAM_MAP
samtools index $OUT_BAM_MAP

# Tag reads with barcodes
python3 TagBam.py \
    --tag XP:Z:$BARCODE \
    $OUT_BAM_MAP \
    $OUT_MAP_SORT_PREFIX"_tagged.bam"
samtools index $OUT_MAP_SORT_PREFIX"_tagged.bam"

# Tag reads wih i7 indexes
python3 TagBam.py \
    --tag XJ:Z:$INDEX1 \
    $OUT_MAP_SORT_PREFIX"_tagged.bam" \
    $OUT_MAP_SORT_PREFIX"_tagged2.bam"
samtools index $OUT_MAP_SORT_PREFIX"_tagged2.bam"

# Tag reads with i5 indexes
python3 TagBam.py \
    --tag XK:Z:$INDEX2 \
    $OUT_MAP_SORT_PREFIX"_tagged2.bam" \
    $OUT_MAP_SORT_PREFIX"_tagged3.bam"
samtools index $OUT_MAP_SORT_PREFIX"_tagged3.bam"

# Tag reads with SRT barcode
python3 TagBamWithSrtBC.py \
    --tag XQ:Z \
    $OUT_MAP_SORT_PREFIX"_tagged3.bam" \
    $OUT_MAP_SORT_PREFIX"_tagged4.bam"
samtools index $OUT_MAP_SORT_PREFIX"_tagged4.bam"

# Tag reads with transposon insertions
python3 AnnotateInsertionSites.py \
    --transposase $TRANSPOSASE \
    -f \
    $OUT_MAP_SORT_PREFIX"_tagged4.bam" \
    $TWO_BIT_PATH \
    $OUT_MAP_SORT_PREFIX"_final.bam"
samtools index $OUT_MAP_SORT_PREFIX"_final.bam"

# Get an (unsorted) list of unique insertions
python3 BamToCallingCard.py -b XP XJ XQ -i $OUT_MAP_SORT_PREFIX"_final.bam" -o $OUT_MAP_SORT_PREFIX"_unsorted.qbed"

# Sort the qbed file
cat $OUT_MAP_SORT_PREFIX"_unsorted.qbed" | sort -k 1,1 -k2,2n > $OUT_MAP_SORT_PREFIX"_final.qbed"

# Clean up
rm $OUT_SAM
rm $OUT_BAM
rm $OUT_STEM".temp1.fastq.gz"
rm $OUT_STEM".temp2.fastq.gz"
rm $OUT_STEM".temp3.fastq.gz"
rm $OUT_BAM_MAP
rm $OUT_BAM_MAP".bai"
rm $OUT_MAP_SORT_PREFIX"_tagged.bam"
rm $OUT_MAP_SORT_PREFIX"_tagged2.bam"
rm $OUT_MAP_SORT_PREFIX"_tagged3.bam"
rm $OUT_MAP_SORT_PREFIX"_tagged.bam.bai"
rm $OUT_MAP_SORT_PREFIX"_tagged2.bam.bai"
rm $OUT_MAP_SORT_PREFIX"_tagged3.bam.bai"
rm $OUT_MAP_SORT_PREFIX"_unsorted.qbed"

#check that srt barcode completed successfully
echo -e "\n\nCHECK FILE "$OUT_STEM"_final_SRT_barcode_distribution.txt FOR SRT BARCODE DISTRIBUTIONS"
echo -e "# insertions\tSRT barcode" > $OUT_STEM"_final_SRT_barcode_distribution.txt" #make header line
awk -F '/' '{print $3}' $OUT_MAP_SORT_PREFIX"_final.qbed" | sort | uniq -c | sort -nr | awk -F" " -v OFS="\t" '{print $1, $2}' >> $OUT_STEM"_final_SRT_barcode_distribution.txt"
