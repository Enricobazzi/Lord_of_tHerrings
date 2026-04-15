#!/bin/bash -l
#SBATCH -A naiss2025-5-565
#SBATCH -p shared
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH -t 0-01:00:00
#SBATCH --mem=12G

## load software
ml samtools

# arguments
sample=${1}

# filter out repeats and inversions
samtools view -hb -L data/diversity/beds/noreps_noinvs.bed \
    data/bams/${sample}.subsampled_3X.bam \
    -o data/diversity/bams/${sample}.subsampled_3X.noreps_noinvs.bam