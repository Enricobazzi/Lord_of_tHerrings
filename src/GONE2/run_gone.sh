#!/bin/bash -l
#SBATCH -A naiss2025-5-565
#SBATCH -p shared
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=3
#SBATCH -t 0-00:20:00
#SBATCH --mem=124G

pop=${1}
i=${2}
g=${3}

if [ ! -f data/GONE2/${pop}.lowmiss_snps.ped ]; then
    echo "Input PED file data/GONE2/${pop}.lowmiss_snps.ped not found!"
    exit 1
fi

src/GONE2/GONE2-1.0.2/gone2 \
    -t 3 \
    -g ${g} \
    -r 2.54 \
    -s 150000 \
    -o data/GONE2/output/${pop}.${i} \
    data/GONE2/${pop}.lowmiss_snps.ped
