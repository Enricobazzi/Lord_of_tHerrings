#!/bin/bash -l
#SBATCH -A naiss2025-5-565
#SBATCH -p shared
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH -t 0-24:00:00
#SBATCH --mem=64G

# load software
ml angsd

# constants
THREADS=8
REF=/cfs/klemming/projects/supr/naiss2024-6-170/analyses/Herring/Enrico/Reference/GCF_900700415.2_Ch_v2.0.2_genomic.fna
OUT=data/diversity/output

# arg
dataset=${1}

# vars
bamlist=data/diversity/bamlists/${dataset}.bamlist
n_ind=$(wc -l < $bamlist)
min_ind=$((n_ind / 2))
maxD=$((n_ind * 3 * 5))
minD=$((n_ind * 3 / 3))
sfs=${OUT}/${dataset}.folded.sfs

# run angsd
angsd -P ${THREADS} -b ${bamlist} -ref ${REF} -anc ${REF} \
    -out ${OUT}/${dataset}.folded.postprob \
    -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -baq 1 -C 50 \
    -minMapQ 30 -minQ 20 -doCounts 1 \
    -GL 1 -doSaf 1 \
    -minInd ${min_ind} -setMaxDepth ${maxD} -setMinDepth ${minD} \
    -doThetas 1 -pest ${sfs}
