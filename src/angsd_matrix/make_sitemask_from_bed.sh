#!/bin/bash -l
#SBATCH -A naiss2025-5-565
#SBATCH -p shared
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH -t 0-12:00:00
#SBATCH --mem=128G

# modules
ml bedtools

# arguments
dataset=${1}
inbed=${2}

# constants
GTLIKE="data/angsd_matrix/gtlike"
SITES="data/angsd_matrix/sites"

# derived variables
sites=$(basename ${inbed} .bed)
maffile=${GTLIKE}/${dataset}.mafs.gz
databed=${SITES}/${dataset}.bed
sitemask=${SITES}/${dataset}.${sites}.sitemask

## transform MAF to BED
zcat ${maffile} | \
    cut -f1,2 | tail -n +2 | \
    awk '{print $1, $2 - 1, $2}' | tr ' ' '\t' \
    > ${databed}

## use intersect to create sitemask
bedtools intersect -c \
    -a ${databed} \
    -b ${inbed} \
    | cut -f4 > ${sitemask}
