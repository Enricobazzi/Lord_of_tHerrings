#!/bin/bash -l
#SBATCH -A naiss2025-5-565
#SBATCH -p shared
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH -t 0-12:00:00
#SBATCH --mem=128G

# arguments
dataset=${1}
sites=${2}

# constants
GTLIKE="data/angsd_matrix/gtlike"
SITES="data/angsd_matrix/sites"

# derived variables
beaglefile=${GTLIKE}/${dataset}.beagle.gz
sitemask=${SITES}/${dataset}.${sites}.sitemask
outbeagle=${GTLIKE}/${dataset}.${sites}.filtered.beagle

## filter beagle file:
# copy header
zcat ${beaglefile} | head -1 > ${outbeagle}
# paste beagle file with sitemask and filter out sites that are not in the bed file
paste <(zcat ${beaglefile} | tail -n +2) <(cat ${sitemask}) | \
    awk '$NF == 1 {NF--; print}' >> ${outbeagle}

