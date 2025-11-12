#!/bin/bash -l
#SBATCH -A naiss2025-5-565
#SBATCH -p shared
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH -t 0-3:30:00
#SBATCH --mem=24G

module load bcftools

BCFS=/cfs/klemming/projects/supr/naiss2024-6-170/analyses/Bcfs_Herring
REPMA=/cfs/klemming/projects/supr/naiss2024-6-170/analyses/Generode_Herring/Reference/GCA_040183275.1_Ch_v3.0_genomic.repma.bed
THR=5
QUAL=30
INDEL_DIST=5

sample=${1}
input_bcf=${BCFS}/${sample}.step1.bcf
output_bcf=${BCFS}/${sample}.step3.bcf
if [[ ! -f ${input_bcf} ]]; then
  echo "Input BCF file for sample ${sample} not found!"
  exit 1
fi

# min is 2, max is 5 times the average depth calculated from the dp field in the bcf
avgDP=$(bcftools view ${input_bcf} | grep -v "#" | cut -f8 | cut -d';' -f1 | cut -d'=' -f2 | awk '{sum += $1} END {print sum/NR}')
# minDP=$(echo "$avgDP / 3" | bc)
minDP=2
maxDP=$(echo "$avgDP * 5" | bc)
echo "Average DP: $avgDP" > ${BCFS}/${sample}.DP.info.txt
echo "Min DP: $minDP" >> ${BCFS}/${sample}.DP.info.txt
echo "Max DP: $maxDP" >> ${BCFS}/${sample}.DP.info.txt

bcftools view -Ob --threads ${THR} \
    -T ${REPMA} ${input_bcf} | \
bcftools filter -Ob --threads ${THR} \
    -g ${INDEL_DIST} | \
bcftools filter -Ob --threads ${THR} \
    -i "(DP4[0]+DP4[1]+DP4[2]+DP4[3])>${minDP} & \
    (DP4[0]+DP4[1]+DP4[2]+DP4[3])<${maxDP} & \
    QUAL>=${QUAL} & INDEL=0" | \
bcftools view -Ob --threads ${THR} \
    -e 'GT="0/1" & (DP4[2]+DP4[3])/(DP4[0]+DP4[1]+DP4[2]+DP4[3]) < 0.2' | \
bcftools view -Ob --threads ${THR} \
    -e 'GT="0/1" & (DP4[2]+DP4[3])/(DP4[0]+DP4[1]+DP4[2]+DP4[3]) > 0.8' \
    -o ${output_bcf}

bcftools index -o ${output_bcf}.csi ${output_bcf}
bcftools stats ${output_bcf} > ${output_bcf}.stats