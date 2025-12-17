#!/bin/bash -l
#SBATCH -A naiss2025-5-565
#SBATCH -p shared
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH -t 1-12:00:00
#SBATCH --mem=128G

## load software
ml angsd

## arguments
# dataset
dataset=${1}
# sitesfile
sitesfile=${2}

angsd \
    -bam data/angsd_matrix/bamlists/${dataset}.bamlist \
    -sites ${sitesfile} \
    -out data/angsd_matrix/output/${dataset} \
    -minMapQ 30 -minQ 20 -GL 1 -doMajorMinor 1 -doMaf 1 -SNP_pval 1e-6 \
    -doIBS 1 -doCounts 1 -doCov 1 -makeMatrix 1 -minMaf 0.001 \
    -remove_bads 1 -only_proper_pairs 1 -uniqueOnly 1 -skipTriallelic 1 \
    -nThreads 8