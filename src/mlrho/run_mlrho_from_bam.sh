#!/bin/bash -l
#SBATCH -A naiss2025-5-565
#SBATCH -p shared
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH -t 0-01:00:00
#SBATCH --mem=32G

# from the GenErode pipeline

# load modules
ml bioinfo-tools mlRho/2.9
ml samtools

# constants
REF=Reference/GCF_900700415.2_Ch_v2.0.2_genomic.fna
OUT=data/mlrho
callable_bed=Reference/GCF_900700415.2_Ch_v2.0.2_genomic.noreps_noinvs.bed
minDP=2
maxDP=8

# arguments
sample=${1}

# variables
input_bam=data/bams/${sample}.subsampled_3X.noreps_noinvs.bam
pro=${OUT}/pros/${sample}.pro
mlrho=${OUT}/output/${sample}.mlrho.txt

# generate PRO from BAM 
samtools mpileup -f ${REF} -q 30 -Q 30 -B -l ${callable_bed} ${input_bam} | \
    awk -v minDP="${minDP}" -v maxDP="${maxDP}" '$4 >=minDP && $4 <=maxDP' | \
    sam2pro -c 5 > ${pro}

# format PRO file
formatPro -c ${minDP} -n ${OUT}/pros/${sample} ${pro}

# run mlRho on formatted PRO file
mlRho -M 0 -I -n ${OUT}/pros/${sample} > ${mlrho}
