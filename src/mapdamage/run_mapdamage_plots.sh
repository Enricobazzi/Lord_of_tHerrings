#!/bin/bash -l
#SBATCH -A naiss2025-5-565
#SBATCH -p shared
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH -t 0-01:00:00
#SBATCH --mem=32G

# load software
ml PDCOLD/23.12 R/4.1.1 mapdamage/2.2.3

# constants
OUT=data/mapdamage
IBAMS=data/bams
REF=Reference/GCF_900700415.2_Ch_v2.0.2_genomic.fna

# arg
sample=${1}

# vars
bam=${IBAMS}/${sample}.subsampled_3X.bam
out=${OUT}/${sample}

# run mapdamage
mapDamage -i ${bam} -r ${REF} -d ${out} -n 1000000 --merge-reference-sequences
