#!/bin/bash -l
#SBATCH -A naiss2025-5-565
#SBATCH -p shared
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH -t 0-12:00:00
#SBATCH --mem=64G

ml bioinfo-tools NGSadmix/32

# constants
THR=10

# arguments
dataset=${1}
k=${2}
i=${3}

# variables
beaglefile=data/gtlike/${dataset}.beagle.gz
out=data/ngsadmix/${dataset}.admix.k_${k}.seed_${i}

# run NGSadmix
NGSadmix \
  -likes ${beaglefile} \
  -K ${k} \
  -P ${THR} \
  -seed ${i} \
  -o ${out}
