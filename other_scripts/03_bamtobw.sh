#!/bin/bash
#
#SBATCH --job-name=bamToBrow
#SBATCH --output=bamToBrow%a.out
#SBATCH --error=bamToBrow%a.err
#SBATCH --mem=20000
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=XXX

# this file takes a set of .bam files, combines and normalizes them, then converts them to .bw files for genome browser.

# Environment set-up
echo "Setting up environment."
eval "$(conda shell.bash hook)"
echo "Conda initialized"
conda activate atac_env

# Merge replicates
samtools merge MOC1veh_merged.bam \
    path/to/Arep1.fastq.gz_ins_size_filt.bam \
    path/to/Arep2.fastq.gz_ins_size_filt.bam

samtools merge MOC1ifn_merged.bam \
    path/to/Brep1.fastq.gz_ins_size_filt.bam \
    path/to/Brep2.fastq.gz_ins_size_filt.bam

samtools merge KOveh_merged.bam \
    path/to/Crep1.fastq.gz_ins_size_filt.bam \
    path/to/Crep2.fastq.gz_ins_size_filt.bam

samtools merge KOifn_merged.bam \
    path/to/Drep1.fastq.gz_ins_size_filt.bam \
    path/to/Drep2.fastq.gz_ins_size_filt.bam

# Sort and index merged BAMs
samtools sort -o MOC1veh_merged.sorted.bam MOC1veh_merged.bam
samtools index MOC1veh_merged.sorted.bam

samtools sort -o MOC1ifn_merged.sorted.bam MOC1ifn_merged.bam
samtools index MOC1ifn_merged.sorted.bam

samtools sort -o KOveh_merged.sorted.bam KOveh_merged.bam
samtools index KOveh_merged.sorted.bam

samtools sort -o KOifn_merged.sorted.bam KOifn_merged.bam
samtools index KOifn_merged.sorted.bam

# Normalize by read depth and create BigWig files
bamCoverage -b MOC1veh_merged.sorted.bam -o MOC1veh_merged.bw --normalizeUsing CPM --binSize 25
bamCoverage -b MOC1ifn_merged.sorted.bam -o MOC1ifn_merged.bw --normalizeUsing CPM --binSize 25
bamCoverage -b KOveh_merged.sorted.bam -o KOveh_merged.bw --normalizeUsing CPM --binSize 25
bamCoverage -b KOifn_merged.sorted.bam -o KOifn_merged.bw --normalizeUsing CPM --binSize 25

mkdir -p norm_bw_files
mv *.bw ./norm_bw_files
