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

thetaStat do_stat ${OUT}/${dataset}.folded.postprob.thetas.idx \
    -win ${WINDOWSIZE} \
    -step ${WINDOWSTEP} \
    -outnames ${OUT}/${dataset}.folded.postprob.thetas.W${WINDOWSIZE}.S${WINDOWSTEP}.gz

