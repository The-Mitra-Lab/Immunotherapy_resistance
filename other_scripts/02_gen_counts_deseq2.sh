#!/bin/bash
#SBATCH --job-name=gen_counts
#SBATCH --output=macs_%j.log
#SBATCH --ntasks=1
#SBATCH -c 4
#SBATCH --mem=20000

# This script takes .bam files as input an outputs a counts matrix for DESeq2

# activate Conda environment
eval "$(conda shell.bash hook)"
conda activate atac_analysis

# define BAM files 
bamfile1="path/to/sample1_rep1.bam"
bamfile2="path/to/sample1_rep2.bam"
bamfile3="path/to/sample2_rep1.bam"
bamfile4="path/to/sample2_rep2.bam"
combinedBed="merged_peak_counts_replicates.bed"

# merge replicates
samtools merge sample2_merged.bam "$bamfile3" "$bamfile4"
samtools sort -o sample2_merged.sorted.bam sample2_merged.bam
samtools index sample2_merged.sorted.bam

samtools merge sample1_merged.bam "$bamfile1" "$bamfile2"
samtools sort -o sample1_merged.sorted.bam sample1_merged.bam
samtools index sample1_merged.sorted.bam

# switch to MACS2 environment
conda deactivate
conda activate macs2_env_clean

# peak calling 
macs2 callpeak \
  -t sample1_merged.sorted.bam \
  -f BAM -g mm \
  -n sample1 --outdir macs2_sample1 \
  --nomodel --shift -100 --extsize 200 -q 0.01

macs2 callpeak \
  -t sample2_merged.sorted.bam \
  -f BAM -g mm \
  -n sample2 --outdir macs2_sample2 \
  --nomodel --shift -100 --extsize 200 -q 0.01

conda deactivate
conda activate atac_analysis

# merge peak sets
cat macs2_sample1/sample1_peaks.narrowPeak macs2_sample2/sample2_peaks.narrowPeak \
  | sort -k1,1 -k2,2n | bedtools merge > merged_peaks.bed

# quantify reads in merged peaks
bedtools multicov \
  -bams "$bamfile1" "$bamfile2" "$bamfile3" "$bamfile4" \
  -bed merged_peaks.bed > "$combinedBed"
