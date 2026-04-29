#!/bin/bash -l
#SBATCH -A naiss2025-5-565
#SBATCH -p shared
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH -t 0-01:30:00
#SBATCH --mem=32G

# load software
ml bcftools

# arg
sample=${1}

# calculate heterozygosity from bcf
sites=$(grep -E "SN.*number of records" data/bcfs/${sample}.step3.bcf.stats | cut -f4)
hets=$(bcftools view -H data/bcfs/${sample}.step3.bcf | grep 0/1 | wc -l)
het=$(echo "scale=6; $hets / $sites" | bc)

# save output
echo "$hets $sites $het" > data/diversity/output/${sample}.het_from_bcf.txt