#!/bin/bash -l
#SBATCH -A naiss2025-22-1002
#SBATCH -p shared
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH -t 0-3:00:00
#SBATCH --mem=24G

module load bcftools

REPMA=/cfs/klemming/projects/supr/naiss2024-6-170/analyses/Generode_Herring/Reference/GCA_040183275.1_Ch_v3.0_genomic.repma.bed
BCFS=/cfs/klemming/projects/supr/naiss2024-6-170/analyses/Bcfs_Herring
THR=5

sample=${1}
input_bcf=${BCFS}/${sample}.step1.bcf
output_bcf=${BCFS}/${sample}.step2.bcf

if [[ ! -f ${input_bcf} ]]; then
  echo "Input BCF file for sample ${sample} not found!"
  exit 1
fi

bcftools view -Ob ${input_bcf} -T ${REPMA} -o ${output_bcf} --threads ${THR}
bcftools index -o ${output_bcf}.csi ${output_bcf}
bcftools stats ${output_bcf} > ${output_bcf}.stats