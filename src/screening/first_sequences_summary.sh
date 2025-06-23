#!/bin/bash

wd=/cfs/klemming/projects/supr/naiss2024-6-170/analyses/Generode_Herring
cd $wd
for sample in $(grep "ND" config/historical_samples_paths.txt | cut -d'_' -f1 | sort -u); do
    mq25=$(
        samtools view -c -F 4 -q 25 \
            results/historical/mapping/GCA_040183275.1_Ch_v3.0_genomic/${sample}_1_L004.sorted.bam
        )
    tot_seq=$(
        grep "Total Sequences" \
            results/historical/trimming/stats/${sample}_1_L004_trimmed_merged_fastqc/fastqc_data.txt | cut -f2
        )
    prop=$(echo "scale=4; $mq25 / $tot_seq" | bc)
    echo "$sample $mq25 $tot_seq $prop"
done

# ND001 8407414 84766400 .0991
# ND006 4015874 104119309 .0385
# ND009 523448 137643471 .0038
# ND010 1062479 108414805 .0098
# ND011 7379716 45304729 .1628
# ND021 17257474 32060333 .5382
# ND022 27418048 40650606 .6744
# ND023 8785531 15575870 .5640
# ND024 7581207 14072037 .5387
# ND026 11099336 18946841 .5858
# ND028 9491707 29263135 .3243
# ND029 11654545 18790811 .6202
# ND030 12106655 19382871 .6246
# ND031 14406811 23506013 .6128
# ND032 6462360 10637917 .6074
# ND034 6091194 20199046 .3015
# ND036 21558474 32208877 .6693
# ND038 4982234 24682187 .2018
# ND039 11949699 18634108 .6412
# ND040 8279266 13586684 .6093
