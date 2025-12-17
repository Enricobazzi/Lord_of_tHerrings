#!/bin/bash -l
#SBATCH -A naiss2025-5-565
#SBATCH -p shared
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH -t 3-00:00:00
#SBATCH --mem=128G

## load software
ml angsd

## arguments
dataset=${1}
# if there is a second argument, it's the sites file
sitesfile=${2}
sitename=$(basename ${sitesfile} .sites)

if [[ -n ${sitesfile} ]]; then
    angsd -GL 2 -nThreads 8 -doGlf 2 -doMajorMinor 1 -doMaf 2 -SNP_pval 1e-6 \
        -minMapQ 30 -minQ 20 -remove_bads 1 -only_proper_pairs 1 -uniqueOnly 1 -skipTriallelic 1 \
        -sites ${sitesfile} \
        -bam data/angsd_matrix/bamlists/${dataset}.bamlist \
        -out data/angsd_matrix/gtlike/${dataset}.${sitename}
else
    echo "No sites file provided, running without -sites option"
    angsd -GL 2 -nThreads 8 -doGlf 2 -doMajorMinor 1 -doMaf 2 -SNP_pval 1e-6 \
        -minMapQ 30 -minQ 20 -remove_bads 1 -only_proper_pairs 1 -uniqueOnly 1 -skipTriallelic 1 \
        -bam data/angsd_matrix/bamlists/${dataset}.bamlist \
        -out data/angsd_matrix/gtlike/${dataset}.wg
fi
