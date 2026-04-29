#!/bin/bash -l
#SBATCH -A naiss2025-5-565
#SBATCH -p shared
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH -t 0-06:00:00
#SBATCH --mem=32G

# from https://github.com/erikrfunk/PPM_env-correlates
# angsd -i ${id}_*.bam -minQ 20 -minmapQ 20 -C 50 -noTrans 1 -doSaf 1 -GL 1 -out ${id}_He_noTrans -P 16 -anc /home/centos/USS/erik/PPM/PPM_reference.fasta -ref /home/centos/USS/erik/PPM/PPM_reference.fasta
# realSFS ${id}_*.saf.idx -maxiter 2000 -tole 1e-16 > ${id}_He_noTrans.ml

# load software
ml angsd

# constants
THREADS=8
REF=/cfs/klemming/projects/supr/naiss2024-6-170/analyses/Herring/Enrico/Reference/GCF_900700415.2_Ch_v2.0.2_genomic.fna
OUT=data/diversity/output

# arg
sample=${1}

# vars
bam=data/diversity/bams/${sample}.subsampled_3X.noreps_noinvs.bam

# run angsd
angsd -P ${THREADS} -i ${bam} -ref ${REF} -anc ${REF} \
    -out ${OUT}/${sample}.het_notrans \
    -minQ 30 -minmapQ 30 -C 50 -noTrans 1 -doSaf 1 -GL 1

# run realSFS
realSFS ${OUT}/${sample}.het_notrans.saf.idx \
    -maxiter 2000 -tole 1e-16 > ${OUT}/${sample}.het_notrans.ml