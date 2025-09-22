# Variant Calling and Filtering

## step 1

```
samples=$(ls /cfs/klemming/projects/supr/naiss2024-6-170/analyses/Bams_Herring/*.bam | rev | cut -d'/' -f1 | rev | cut -d'.' -f1 | sort -u | grep "MHER")

for sample in ${samples[*]}; do
    sbatch \
        --job-name=${sample}.step1 \
        --output=logs/calling_and_filtering/step1.${sample}.out \
        --error=logs/calling_and_filtering/step1.${sample}.err \
        src/calling_and_filtering/step1_mpileup_call.sh \
        ${sample}
done
```

## step 2

```
samples=$(ls /cfs/klemming/projects/supr/naiss2024-6-170/analyses/Bcfs_Herring/*.step1.bcf | rev | cut -d'/' -f1 | rev | cut -d'.' -f1 | sort -u | grep "MHER")

for sample in ${samples[*]}; do
    sbatch \
        --job-name=${sample}.step2 \
        --output=logs/calling_and_filtering/step2.${sample}.out \
        --error=logs/calling_and_filtering/step2.${sample}.err \
        src/calling_and_filtering/step2_intersect_repma.sh \
        ${sample}
done
```

## step 3

```
samples=$(ls /cfs/klemming/projects/supr/naiss2024-6-170/analyses/Bcfs_Herring/*.step2.bcf | rev | cut -d'/' -f1 | rev | cut -d'.' -f1 | sort -u)

for sample in ${samples[*]}; do
    sbatch \
        --job-name=${sample}.step3 \
        --output=logs/calling_and_filtering/step3.${sample}.out \
        --error=logs/calling_and_filtering/step3.${sample}.err \
        src/calling_and_filtering/step3_filter_indels_depth_qual_imbal.sh \
        ${sample}
done
```

## step 4 - create datasets

### modern samples: baltic - transition - NE-atlantic

```
rm data/datasets/modern_samples_bcfs.txt
BCFS=/cfs/klemming/projects/supr/naiss2024-6-170/analyses/Bcfs_Herring
# modern samples:
samples=($(awk -F, '$6 > 1940' data/samples_table.csv | grep -v "ND" | grep "harengus" | grep -E "NE-ATL|BALTIC|TRANS" | cut -d',' -f1 | grep -v "sample_id"))
for sample in ${samples[*]}; do
    echo "${BCFS}/${sample}.step3.bcf" >> data/datasets/modern_samples_bcfs.txt
done
```

### merge!
```
input_dataset=data/datasets/modern_samples_bcfs.txt
output_bcf=data/modern_samples.step3.bcf

sbatch \
    --job-name=merge_modern \
    --output=logs/calling_and_filtering/merge_modern.out \
    --error=logs/calling_and_filtering/merge_modern.err \
    src/calling_and_filtering/merge_bcfs.sh \
    ${input_dataset} ${output_bcf}
```