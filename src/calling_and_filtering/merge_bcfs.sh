#!/bin/bash -l
#SBATCH -A naiss2025-22-1002
#SBATCH -p shared
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH -t 1-00:00:00
#SBATCH --mem=50G

module load bcftools

THR=20
input_dataset=${1}
output_bcf=${2}

bcftools merge -m snps -Ob -o ${output_bcf} --file-list ${input_dataset}
bcftools index -o ${output_bcf}.csi ${output_bcf}
bcftools stats ${output_bcf} > ${output_bcf}.stats
