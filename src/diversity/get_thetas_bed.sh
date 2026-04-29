#!/bin/bash -l
#SBATCH -A naiss2025-5-565
#SBATCH -p shared
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH -t 0-08:00:00
#SBATCH --mem=64G

# load software
ml angsd

# constants
OUT=data/diversity/output

# arg
dataset=${1}

# run thetaStat print and convert to bed format
thetaStat print ${OUT}/${dataset}.folded.thetas.idx \
    > ${OUT}/${dataset}.folded.thetas.tmp

awk '{print $1, $2-1, $2, $3, $4}' ${OUT}/${dataset}.folded.thetas.tmp | \
    tr ' ' '\t' > ${OUT}/${dataset}.folded.thetas.bed

rm ${OUT}/${dataset}.folded.thetas.tmp