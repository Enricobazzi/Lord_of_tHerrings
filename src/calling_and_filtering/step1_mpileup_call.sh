#!/bin/bash -l
#SBATCH -A naiss2025-22-1002
#SBATCH -p shared
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH -t 1-12:00:00
#SBATCH --mem=24G

module load bcftools

REF=/cfs/klemming/projects/supr/naiss2024-6-170/analyses/Generode_Herring/Reference/GCA_040183275.1_Ch_v3.0_genomic.fna
THR=10
BAMS=/cfs/klemming/projects/supr/naiss2024-6-170/analyses/Bams_Herring
BCFS=/cfs/klemming/projects/supr/naiss2024-6-170/analyses/Bcfs_Herring

sample=${1}

if [[ -f ${BAMS}/${sample}.merged.rmdup.merged.realn.rescaled.bam ]]; then
  input_bam=${BAMS}/${sample}.merged.rmdup.merged.realn.rescaled.bam
elif [[ -f ${BAMS}/${sample}.merged.rmdup.merged.realn.bam ]]; then
  input_bam=${BAMS}/${sample}.merged.rmdup.merged.realn.bam
else
  echo "BAM file for sample ${sample} not found!"
  exit 1
fi

output_bcf=${BCFS}/${sample}.step1.bcf

bcftools mpileup -Ou -Q 30 -q 30 -B -f ${REF} ${input_bam} | bcftools call -c -M -O b --threads ${THR} -o ${output_bcf}
bcftools index -o ${output_bcf}.csi ${output_bcf}
bcftools stats ${output_bcf} > ${output_bcf}.stats
