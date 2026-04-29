#!/bin/bash -l
#SBATCH -A naiss2025-5-565
#SBATCH -p shared
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=3
#SBATCH -t 0-00:20:00
#SBATCH --mem=124G

pop=${1}
i=${2}

src/GONE2/gone2 \
    -t 1 \
    -x \
    -g 0 \
    -r 2.54 \
    -l 0.01 \
    -o data/GONE2/output/${pop}.${i} \
    data/GONE2/input/${pop}.lowmiss_snps.vcf

src/GONE2/currentNe2/currentne2 \
    -t 1 \
    -x \
    -r 2.54 \
    -o data/GONE2/output/${pop}.${i} \
    data/GONE2/input/${pop}.lowmiss_snps.vcf
