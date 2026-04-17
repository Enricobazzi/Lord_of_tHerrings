#!/bin/bash -l
#SBATCH -A naiss2025-5-565
#SBATCH -p shared
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH -t 0-08:00:00
#SBATCH --mem=64G

# load software
ml angsd

# constants
THREADS=8
OUT=data/diversity/output
WINDOWSIZE=10000
WINDOWSTEP=5000

# arg
dataset=${1}

# vars
saf=${OUT}/${dataset}.folded.saf.idx
sfs=${OUT}/${dataset}.folded.sfs

# run realSFS saf2theta and thetaStat
realSFS saf2theta ${saf} -cores ${THREADS} -fold 1 -sfs ${sfs} -outname ${OUT}/${dataset}.folded
thetaStat do_stat ${OUT}/${dataset}.folded.thetas.idx \
    -win ${WINDOWSIZE} \
    -step ${WINDOWSTEP} \
    -outnames ${OUT}/${dataset}.folded.thetas.W${WINDOWSIZE}.S${WINDOWSTEP}

