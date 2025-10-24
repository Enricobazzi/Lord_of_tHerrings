#!/bin/bash -l
#SBATCH -A naiss2025-22-1002
#SBATCH -p shared
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH -t 0-01:30:00
#SBATCH --mem=24G

## load software
ml samtools

REPMA='/cfs/klemming/projects/snic/naiss2024-6-170/analyses/Generode_Herring/Reference/GCA_040183275.1_Ch_v3.0_genomic.repma.bed'
ibam=${1}
bamfolder=$(dirname ${ibam})
sample=$(basename ${ibam} | cut -d'.' -f1)

samtools depth -a -q 30 -Q 25 -b ${REPMA} ${ibam} | \
    awk ' {
        v = $3
        count[v]++
        sum += v
        sumsq += v*v
        n++
    } END {
        mean = sum/n
        stddev = sqrt(sumsq/n - mean*mean)
        printf("# mean=%f median=NA stddev=%f\n", mean, stddev)
        PROCINFO["sorted_in"] = "@ind_num_asc"
        for (v in count) print v, count[v]
    }' > ${bamfolder}/${sample}.depth
# Note: median calculation is not implemented here