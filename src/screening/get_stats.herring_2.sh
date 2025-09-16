#!/bin/bash -l
#SBATCH -A naiss2024-5-520
#SBATCH -p shared
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH -t 0-01:00:00

## load software
ml PDC
ml bioinfo-tools
ml samtools

# directory containing raw sequencing stats
raw_reads='data/raw_reads_symlinks/historical/stats'
# directory containing trimmed reads
trimming='results/historical/trimming/stats'
# directory containing mapped reads stats
mapped='results/historical/mapping/GCA_040183275.1_Ch_v3.0_genomic'
# path to reference repeat masked bed file
bed='//cfs/klemming/projects/snic/naiss2024-6-170/analyses/Generode_Herring/Reference/GCA_040183275.1_Ch_v3.0_genomic.repma.bed'

sample=${1}

# Get total sequences from fastqc data of fastq files
total_seqs=$(grep "Total Sequences" $raw_reads/${sample}_*_R1_fastqc/fastqc_data.txt | awk '{sum += $3} END {print sum}')

# get the total number of sequences after trimming and merging
total_merged=$(grep "Total Sequences" $trimming/${sample}_*_trimmed_merged_fastqc/fastqc_data.txt | awk '{sum += $3} END {print sum}')

# Get mapped reads from sorted_bam stats
total_mapped=$(grep -m1 mapped $mapped/stats/bams_sorted/${sample}_*.sorted.bam.stats.txt | cut -d':' -f2 | cut -d' ' -f1 | awk '{sum += $1} END {print sum}')

# Get unique reads from rmdup_bam stats
total_uniq=$(grep -m1 mapped $mapped/stats/bams_rmdup/${sample}_*.merged.rmdup.bam.stats.txt | cut -d':' -f2 | cut -d' ' -f1 | awk '{sum += $1} END {print sum}')

# get mq 25 reads
total_MQ25=$(samtools view -c -F 4 -q 25 $mapped/${sample}.merged.rmdup.merged.realn.rescaled.bam)

# percentage of merged ND031_1_L004.sorted.bam
proportion_merged=$(awk "BEGIN {printf \"%.5f\", $total_merged / $total_seqs}")

# get endogenous
endogenous=$(awk "BEGIN {printf \"%.5f\", $total_mapped / $total_seqs}")

# merged endogenous
merged_endogenous=$(awk "BEGIN {printf \"%.5f\", $total_mapped / $total_merged}")

# get complexity
complexity=$(awk "BEGIN {printf \"%.5f\", $total_uniq / $total_mapped}")

# get genome recovery rate
grr=$(awk "BEGIN {printf \"%.5f\", $total_MQ25 / $total_seqs}")

# get total coverage
total_cov=$(samtools depth -a -q 30 -Q 25 -b $bed $mapped/${sample}.merged.rmdup.merged.realn.rescaled.bam | awk '{{sum+=$3}} END {{ print sum/NR }}' | awk '{print$1}')

echo ${sample} $total_seqs $total_merged $total_mapped $total_uniq $total_MQ25 $proportion_merged $endogenous $merged_endogenous $complexity $grr $total_cov > seqstats/${sample}_seqstats.txt
