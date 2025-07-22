#!/bin/bash
#
#SBATCH --job-name=atac_proc
#SBATCH --output=atac%a.out
#SBATCH --error=atac%a.err
#SBATCH --mem=64000
#SBATCH --mail-type=ALL
#SBATCH --mail-user=XXX
################################################################################################################
# ATAC-seq peak calling
# Reilly Sample, Mitra Lab, WashU
# This script is an adapted peak caller from Howard Chang Lab.
# Description: Performs Nextera cutadapt step, alignment, read de-duplication, flagstat QC, nucleosome-free DNA
# restriction, Tn5 coordinate shift, narrow and broak peak calling, bigwig file generation for browser viewing.
# Input files: .fastq files (Read 1+2)
# Output files: indexed bam files, begraph, bed, flagstat QC report, fragment length distribution plot and .txt
# report, log and error files.
# Instructions: Run this script with your .fastq.gz files in ./raw. Make sure you  have
# parallel installed in the atac_env conda environment. Change the parameters
# under "## Default settings:" when applicable. Change email address above.
# Usage: sbatch call_atac_peaks.sh
################################################################################################################

## Default settings:
MAXJOBS=32 # Maximum parallel jobs for alignment
FASTQFOLDER="/path/to/fastq/raw" path here
READ1EXTENSION="R1_001.fastq.gz"
INSERT_SIZE_LIMIT="100" # limit to 100 when interested in nucleosome-free regions
MACSQVALUE=".05" # Q-value for MACS2 peak calling
GENOME="mm10" # hg38 or mm10
MIN_LENGTH=36 # minimum length of reads to keep following cutadapt

# Environment set-up
echo "setting up env."
eval "$(conda shell.bash hook)"
echo "conda loaded"
conda activate atac_env
echo "env activated"

mkdir -p trimmed
cd trimmed

# Cutadapt step:
# -a removes 'CTGTCTCTTATACACATCT' adapter sequence from R1
# -A removes 'CTGTCTCTTATACACATCT' adapter sequence from R2
# Bases at the 3' end will be trimmed until a base with a Phred score ≥ 15 is encountered.
# Bases at the 5' end will be trimmed until a base with a Phred score ≥ 10 is found.
# -o is output file for R1; -p is output file for R2.

find $FASTQFOLDER -name "*${READ1EXTENSION}" | while read file; do
    xbase=$(basename $file)
    echo "cutadapt -a CTGTCTCTTATACACATCT -A CTGTCTCTTATACACATCT \
        --quality-cutoff=15,10 --minimum-length=$MIN_LENGTH -o Trimmed_"$xbase" \
        -p Trimmed_"${xbase/R1/R2}" "$file" "${file/R1/R2}" > "$xbase"_cutadapt.log" >> 1_cutadaptcommands.txt
done

echo "Trimming Reads"
parallel -j $MAXJOBS < 1_cutadaptcommands.txt 2> 1_cutadaptcommands.err &> 1_cutadaptcommands.log

# Output files: 'Trimmed_*' (trimmed .fastq files (.gz compressed))

# Run FastQC
find . -maxdepth 1 -name "*.gz" | while read file; do
    xbase=$(basename "$file")
    outdir="${xbase%.*}_fastqc"

    mkdir -p "$outdir"  # Ensure the directory is created
    echo "fastqc -o \"$outdir\" \"$file\"" >> 2_fastqc_commands.txt
done

parallel -j $MAXJOBS < 2_fastqc_commands.txt

echo "FastQC complete."

# Run alignment using Bowtie2
# Input: trimmed .fastq files

cd ..
mkdir -p aligned
cd aligned

if [[ $GENOME == "hg38" ]]; then
    INDEX="/path/to/GRCh38_noalt_as"
elif [[ $GENOME == "mm10" ]]; then
    INDEX="/path/to/mm10"
else
    echo "Error: Unsupported genome '$GENOME'"
    exit 1
fi

find ../trimmed -name "*${READ1EXTENSION}" | while read file; do
    xbase=$(basename "$file")
    bam_out="${xbase%.*}.bam"

    echo "bowtie2 --very-sensitive-local --no-unal --no-mixed --no-discordant --phred33 -I 10 -X 700 \
        -x \"$INDEX\" -1 \"$file\" -2 \"${file/R1/R2}\" | \
        samtools view -u - | samtools sort -o \"$bam_out\"; \
        samtools index \"$bam_out\"" >> 2_alignCommands.txt
done

parallel -j $MAXJOBS < 2_alignCommands.txt 2> 2_alignCommands.err &> 2_alignCommands.log

# Output: *.bam

# Use samtools to remove mitochondrial reads and duplicate reads

find . -name "*fastq.bam" | while read file; do
    echo "samtools view -b \"$file\" chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY | \
        samtools sort > \"${file%.*}\"_nochrM.bam ; samtools index \"${file%.*}\"_nochrM.bam" >> 3_remove_chrM.txt
done

parallel -j $MAXJOBS < 3_remove_chrM.txt

# Output: *_nochrM.bam

# Calculate flagstat statistics

find . -name "*_nochrM.bam" | while read file; do
    echo "samtools flagstat \"$file\" > \"${file}_flag.stats\"" >> 5_flagstat_calc.txt
done

parallel -j $MAXJOBS < 5_flagstat_calc.txt 2> 5_flagstat_calc.err &> 5_flagstat_calc.log

# Output: *.stats

# Extract properly paired and uniquely mapped reads

find . -name "*_nochrM.bam" | while read file; do
    echo "samtools view -h -b -q 10 -f 2 \"$file\" > \"${file%.*}\"_qffilter.bam"
done >> 6_extractProperlyPaired.txt

echo "(6/17) Removing Unpaired and Low-Quality Reads"
parallel -j $MAXJOBS < 6_extractProperlyPaired.txt \
    2> 6_extractProperlyPaired.err &> 6_extractProperlyPaired.log

# Output: *_qffilter.bam

# Clean Up Intermediate Files
rm *fastq.bam
rm *nochrM.bam
rm *bai

# Index the filtered file
# Creates a .bai file for each *_qffilter.bam
find . -name "*_qffilter.bam" | while read file; do
    echo "samtools index \"$file\""
done >> 7_indexCommands.txt

parallel -j $MAXJOBS < 7_indexCommands.txt 2> 7_indexCommands.err &> 7_indexCommands.log

# Calculate and plot the insert size distribution.
echo "starting to plot inserts"
find . -name "*_qffilter.bam" | while read file; do
    frag_length_file="${file}_frag_length.txt"
    echo "samtools view \"$file\" | awk '\$9>0' | cut -f 9 | sort | uniq -c | sort -b -k2,2n | sed -e 's/^[ \t]*//' > \"$frag_length_file\""
done >> 16_fragLengthCommands.txt

# Run fragment length count generation in parallel
parallel -j "$MAXJOBS" < 16_fragLengthCommands.txt 2> 16_fragLengthCommands.err &> 16_fragLengthCommands.log

cd ..

mkdir -p QC

find ./aligned -name "*_qffilter.bam_frag_length.txt" | while read frag_length_file; do
    echo "Rscript plot_insert_size.R \"$frag_length_file\""
done >> 16_plotCommands.txt

# Run R script plotting in parallel
parallel -j "$MAXJOBS" < 16_plotCommands.txt 2> 16_plotCommands.err &> 16_plotCommands.log

cd aligned

# Limit peak calling to certain fragment sizes.
# This value may dpeend on your data. Generally, <= 100 bp will correspond to nucleosome-free DNA.
# We isolate nucleosome-free DNA when we are interested in identifying enhancers, promoters, and TSS.

find . -name "*_qffilter.bam" | while read file; do
    # Filter based on insert size <= INSERT_SIZE_LIMIT and create a new filtered BAM file
    echo "samtools view -h \"$file\" | awk '\$9 <= $INSERT_SIZE_LIMIT' | samtools view -bS - > \"${file%.bam}_ins_size_filt.bam\""
done >> 17_filterInsertSizeCommands.txt

# Run the filtering commands in parallel
parallel -j $MAXJOBS < 17_filterInsertSizeCommands.txt 2> 17_filterInsertSizeCommands.err &> 17_filterInsertSizeCommands.log

# Shift coordinates (Tn5 creates a staggered cut):
# For "+" strand reads (forward reads), the actual cut site is 4 bp downstream (shift +4 bp).
# For "-" strand reads (reverse reads), the actual cut site is 5 bp upstream (shift -5 bp).

find . -name "*_ins_size_filt.bam" | while read file; do
    echo "samtools sort -n $file | \
        bedtools bamtobed -bedpe -i stdin | \
        awk -v OFS='\t' '{if(\$9==\"+\"){print \$1,\$2+4,\$6+4} \
                          else if(\$9==\"-\"){print \$1,\$2-5,\$6-5}}' \
        > ${file%.*}.tagAlign" >> 9_tagAlignGenerate.txt
done

parallel -j $MAXJOBS < 9_tagAlignGenerate.txt 2> 9_tagAlignGenerate.err &> 9_tagAlignGenerate.log

# Output: *.tagAlign

# Call peaks with "summits" (i.e. sharp peaks)
# Use the *.narrowPeak calls for identifying TF motifs and mapping regulatory elements.

find . -name "*tagAlign" | while read file; do
    echo "macs2 callpeak -t $file \
        -f BED \
        -n ${file%.*} \
        -q $MACSQVALUE \
        --nomodel \
        --shift 37 \
        --extsize 73 \
        -B \
        --keep-dup 1 \
        --call-summits" >> 10_macs2Call.txt
done

parallel -j $MAXJOBS < 10_macs2Call.txt 2> 10_macs2Call.err &> 10_macs2Call.log

# Outputs: *.narrowPeak, *_treat_pileup.bdg

# Call broader peaks for chromatin states.

find . -name "*tagAlign" | while read file; do
    echo "macs2 callpeak -t $file -f BED -n ${file%.*}_broadpeak -q $MACSQVALUE \
        --nomodel --shift 37 --extsize 73 -B --keep-dup 1" >> 10_macs2Call_broadpeak.txt
done

parallel -j $MAXJOBS < 10_macs2Call_broadpeak.txt 2> 10_macs2Call_broadpeak.err &> 10_macs2Call_broadpeak.log

# Output: *.broadpeak

# Remove blacklisted regions from narrowPeak files

if [[ $GENOME == "hg38" ]]; then
    find . -name "*.narrowPeak" | while read file ; do
        echo "bedtools intersect -v -a "$file" -b /path/to/GRCh38_blacklist.bed > "${file%.*}"_noBL.narrowPeak" >> 11_removeBlacklist.txt
    done
elif [[ $GENOME == "mm10" ]]; then
    find . -name "*.narrowPeak" | while read file ; do
        echo "bedtools intersect -v -a "$file" -b /path/to/mm10_blacklist.bed > "${file%.*}"_noBL.narrowPeak" >> 11_removeBlacklist.txt
    done
fi

parallel -j $MAXJOBS < 11_removeBlacklist.txt 2> 11_removeBlacklist.err &> 11_removeBlacklist.log

# Output: *_noBL.narrowPeak

# Generate BigWig file for genome browsing
if [[ $GENOME == "hg38" ]]; then
    find . -name "*_treat_pileup.bdg" | while read file; do
        echo "sort -k1,1 -k2,2n "$file" > "${file%.*}"_sorted.bdg ; \
        bedGraphToBigWig "${file%.*}"_sorted.bdg /path/to/hg38.chrom.sizes "${file%.*}"_sorted.bw" \
        >> 12_generateBW.txt
    done

elif [[ $GENOME == "mm10" ]]; then
    find . -name "*_treat_pileup.bdg" | while read file; do
        echo "sort -k1,1 -k2,2n "$file" > "${file%.*}"_sorted.bdg ; \
        bedGraphToBigWig "${file%.*}"_sorted.bdg /path/to/mm10.chrom.sizes "${file%.*}"_sorted.bw" \
        >> 12_generateBW.txt
    done
fi

parallel -j $MAXJOBS < 12_generateBW.txt 2> 12_generateBW.err &> 12_generateBW.log

# Output: *.bw

# Clean up directories by storing results and removing intermediate files

mkdir -p bw
mkdir -p bam
mkdir -p bdg
mkdir -p bed
mkdir -p narrowpeak
mkdir -p xls_output
mkdir -p flagstats
mkdir -p frag_lengths
mkdir -p tagAlign

mv *.bw bw
mv *.bam bam
mv *.bai bam
mv *.bdg bdg
mv *.bed bed
mv *.narrowPeak narrowpeak
mv *.xls xls_output
mv *.stats flagstats
mv *length.txt frag_lengths
mv *.tagAlign tagAlign

rm *.log
rm *.err
rm *.txt

cd ..
rm 16_plot*
cd trimmed

rm 1_cutadaptcommands*
mv *.log ..
cd ..
rm -r trimmed
