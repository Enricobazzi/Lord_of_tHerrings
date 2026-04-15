#!/bin/bash -l
#SBATCH -A naiss2025-5-565
#SBATCH -p shared
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH -t 0-06:00:00
#SBATCH --mem=64G

# load software
ml angsd

# constants
THREADS=8
OUT=data/diversity/output

# arg
dataset=${1}

# vars
saf=${OUT}/${dataset}.folded.saf.idx
sfs=${OUT}/${dataset}.folded.sfs

# run realSFS
realSFS ${saf} -cores ${THREADS} -fold 1 > ${sfs}
