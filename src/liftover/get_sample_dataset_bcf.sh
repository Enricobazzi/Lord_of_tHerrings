#!/bin/bash -l
#SBATCH -A naiss2025-22-1002
#SBATCH -p shared
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH -t 0-00:20:00
#SBATCH --mem=24G

module load bcftools

BCFS=/cfs/klemming/projects/supr/naiss2024-6-170/analyses/Bcfs_Herring
THR=5

sample=${1}
dataset=${2}

bed=beds/${dataset}_snps.v3.bed
input_bcf=${BCFS}/${sample}.step3.bcf
output_bcf=bcfs/${sample}.${dataset}.bcf

bcftools view -Ob --threads ${THR} \
    ${input_bcf} \
    -T ${bed} \
    -o ${output_bcf}

bcftools index -o ${output_bcf}.csi ${output_bcf}
