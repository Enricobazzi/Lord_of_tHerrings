#!/bin/bash -l
#SBATCH -A naiss2024-5-520
#SBATCH -p shared
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH -t 0-00:20:00
#SBATCH -J mapping_stats

## information

# this script collects and calculates summary statistics from historical samples in the GenErode pipeline
# the script is run per sample, and adds output to the file "historical_sequencing_stats.txt"
# in order to get a sample list ("samp_list")from the metadata table, you can use sed as below:
# awk '{print$1}' ../GenErode/config/historical_samples.txt | sed -r -e 's/(_[0-9]+_L[0-9]+)//g' > samp_list
# the script is run within the main GenErode directory
# stats based on mapping quality (MQ) 30.

## edits to make before running:
# change the {ref} to match your reference genome name
# change the path to your reference bed file
# change the mitogenome scaffold
# change the path to the scripts directory

## usage of the script
# cd /path/to/GenErode/
# for i in $(cat samp_list); do sbatch collect_stats.sh $i; done

# alternative usage for a single sample
# sbatch collect_stats.sh {samp}

## directories and files

# directory containing raw sequencing stats
raw_reads='data/raw_reads_symlinks/historical/stats'
# directory containing trimmed reads
trimming='results/historical/trimming/stats'
# directory containing mapped reads stats
mapped='results/historical/mapping/GCA_040183275.1_Ch_v3.0_genomic'
# path to reference repeat masked bed file
bed='//cfs/klemming/projects/snic/naiss2024-6-170/analyses/Generode_Herring/Reference/GCA_040183275.1_Ch_v3.0_genomic.repma.bed'
# mitogenome scaffold
mito=''
# path to directory containing calculate.awk
scripts='seqstats'

## load software
ml PDC
ml bioinfo-tools
ml samtools

## summary statistics

# Get total sequences from fastqc data of fastq files

total_seqs=$(grep "Total Sequences" $raw_reads/${1}_*_R1_fastqc/fastqc_data.txt | awk '{print$3}')

# get the total number of sequences after trimming and merging

total_merged=$(grep "Total Sequences" $trimming/${1}_*_trimmed_merged_fastqc/fastqc_data.txt | awk '{print$3}')

# Get mapped reads from sorted_bam stats

total_mapped=$(awk 'NR==7 {print$1}' $mapped/stats/bams_sorted/${1}_*.sorted.bam.stats.txt)

# Get unique reads from rmdup_bam stats

total_uniq=$(awk 'NR==7 {print$1}' $mapped/stats/bams_rmdup/${1}_*.merged.rmdup.bam.stats.txt)

# get mq 25 reads

total_MQ25=$(samtools view -c -F 4 -q 25 $mapped/${1}.merged.rmdup.merged.realn.bam)

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

total_cov=$(samtools depth -a -q 30 -Q 25 -b $bed $mapped/${1}.merged.rmdup.merged.realn.bam | awk '{{sum+=$3}} END {{ print sum/NR }}' | awk '{print$1}')

# get mito coverage

#mito_cov=$(samtools depth -a -q 30 -Q 25 -r $mito $mapped/${1}.merged.rmdup.merged.realn.bam | awk '{{sum+=$3}} END {{ print sum/NR }}' | awk '{print$1}')

# get read lengths
# prints min max median mean

read_length=$(samtools view -q 25 $mapped/${1}.merged.rmdup.merged.realn.bam | awk '{print length($10)}' | 
    awk '{a[NR]=$1; sum+=$1; if(min==""){min=max=$1}; if($1>max){max=$1}; if($1<min){min=$1}} END { 
        n=asort(a); 
        median = (n % 2 ? a[(n+1)/2] : (a[n/2] + a[n/2+1])/2); 
        mean = sum / n; 
        print "min="min,"max="max,"median="median,"mean="mean 
    }' | tr ' ' ',')

# print file (add $mito_cov after $total_cov if you want to include mito coverage - uncomment above)

echo ${1} $total_seqs $total_merged $total_mapped $total_uniq $total_MQ25 $proportion_merged $endogenous $merged_endogenous $complexity $grr $total_cov $read_length >> seqstats/${1}_ancient_screening_stats.txt
