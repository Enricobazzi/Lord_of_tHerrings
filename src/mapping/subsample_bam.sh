#!/bin/bash -l
#SBATCH -A naiss2025-5-565
#SBATCH -p shared
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH -t 0-02:30:00
#SBATCH --mem=24G

## load software
ml samtools

## arguments
ibam=${1}
obam=${2}
fraction=${3}

if [[ ${fraction} < 1 ]]; then 
    samtools view -h -b -s ${fraction} -@ 1 -o ${obam} ${ibam}
else
    cp ${ibam} ${obam}
fi

samtools index ${obam}
