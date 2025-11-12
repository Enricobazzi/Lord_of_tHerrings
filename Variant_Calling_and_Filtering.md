# Variant Calling and Filtering

copy bams from generode folder in scratch to common folder
```
samples=$(cat config/historical_samples_paths.txt | cut -d'_' -f1 | grep -v samplename | sort -u)
for sample in ${samples[*]}; do
    echo ${sample}
    cp results/historical/mapping/GCA_040183275.1_Ch_v3.0_genomic/${sample}.merged.rmdup.merged.realn.rescaled.bam /cfs/klemming/projects/supr/naiss2024-6-170/analyses/Bams_Herring/
    cp results/historical/mapping/GCA_040183275.1_Ch_v3.0_genomic/${sample}.merged.rmdup.merged.realn.rescaled.bam.bai /cfs/klemming/projects/supr/naiss2024-6-170/analyses/Bams_Herring/
done
```

## step 1

```
samples=$(ls /cfs/klemming/projects/supr/naiss2024-6-170/analyses/Bams_Herring/*.bam | rev | cut -d'/' -f1 | rev | cut -d'.' -f1 | sort -u | grep "MHER")
samples=$(cat config/historical_samples_paths.txt | cut -d'_' -f1 | grep -v samplename | sort -u)

# for sample in ND157 ND158 ND164 ND165 ND169 ND184 ND319 ND325 ND463; do
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

## step 2-3

```
samples=$(ls /cfs/klemming/projects/supr/naiss2024-6-170/analyses/Bcfs_Herring/ND*.step1.bcf | rev | cut -d'/' -f1 | rev | cut -d'.' -f1 | sort -u | grep -v "ND0")

for sample in ${samples[*]}; do
    sbatch \
        --job-name=${sample}.step3 \
        --output=logs/calling_and_filtering/step3.${sample}.out \
        --error=logs/calling_and_filtering/step3.${sample}.err \
        src/calling_and_filtering/step3_filter_repma_indels_depth_qual_imbal.sh \
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
### vasa ship
```
BCFS=/cfs/klemming/projects/supr/naiss2024-6-170/analyses/Bcfs_Herring
for sample in ND364 ND365 ND374 ND375 ND387 AAL1 AAL2 AAL3 AK1 AK2 Fehmarn3 Fehmarn44 Fehmarn6 Gavle100 Gavle54 Gavle98 MHER016 MHER017 MHER018 NSSH33 NSSH34 NSSH36 MHER034 MHER035 MHER036; do
    echo "${BCFS}/${sample}.step3.bcf"
done > data/datasets/vasaship_quick_bcfs.txt
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

input_dataset=data/datasets/vasaship_quick_bcfs.txt
output_bcf=data/vasaship_quick.step3.bcf
sbatch \
    --job-name=merge_vasa \
    --output=logs/calling_and_filtering/merge_vasa.out \
    --error=logs/calling_and_filtering/merge_vasa.err \
    src/calling_and_filtering/merge_bcfs.sh \
    ${input_dataset} ${output_bcf}
```

### kk

```
plink --bcf <(bcftools view -Ou -i 'MAF > 0.02' data/vasaship_quick.step3.bcf) \
    --double-id --allow-extra-chr --set-missing-var-ids @:# \
    --indep-pairwise 50 10 0.1 \
    --out data/pca/vasaship_quick

plink --bcf <(bcftools view -Ou -i 'MAF > 0.02' data/vasaship_quick.step3.bcf) \
    --double-id --allow-extra-chr --set-missing-var-ids @:# \
    --extract data/pca/vasaship_quick.prune.in --pca 'header' \
    --out data/pca/vasaship_quick
```

```
BCFS=/cfs/klemming/projects/supr/naiss2024-6-170/analyses/Bcfs_Herring
samples=(ND364 ND365 ND374 ND375 ND387)

sample=ND364
remove=$(for s in ${samples[*]}; do if [ $s != $sample ]; then echo $s; fi; done | tr '\n' '|' | sed 's/.$//')

bcftools view -H ${BCFS}/${sample}.step3.bcf |
    awk '{print $1, $2-1, $2}' | tr ' ' '\t' | bedtools merge -i stdin \
> data/beds/${sample}.step3.bed

bedtools intersect -header \
    -a <(bcftools view data/vasaship_quick.step3.bcf) \
    -b data/beds/${sample}.step3.bed \
> data/vasaship_quick.step3.${sample}.vcf

plink --vcf data/vasaship_quick.step3.${sample}.vcf \
    --double-id --allow-extra-chr --set-missing-var-ids @:# \
    --indep-pairwise 50 10 0.1 --maf \
    --out data/pca/vasaship_quick.${sample}

plink --vcf data/vasaship_quick.step3.${sample}.vcf \
    --double-id --allow-extra-chr --set-missing-var-ids @:# \
    --remove <(grep -E $remove data/pca/vasaship_quick.${sample}.nosex) \
    --extract data/pca/vasaship_quick.${sample}.prune.in --pca 'header' \
    --out data/pca/vasaship_quick.${sample}

plink --vcf data/vasaship_quick.step3.${sample}.vcf \
    --double-id --allow-extra-chr --set-missing-var-ids @:# \
    --remove <(grep -E $remove data/pca/vasaship_quick.${sample}.nosex) \
    --extract data/pca/vasaship_quick.${sample}.prune.in --recode A \
    --out data/pca/vasaship_quick.${sample}


plink --vcf data/vasaship_quick.step3.${sample}.vcf \
    --double-id --allow-extra-chr --set-missing-var-ids @:# \
    --remove <(grep -E "ND365|ND374|ND375|ND387|ND364" data/pca/vasaship_quick.${sample}.nosex) \
    --extract data/pca/vasaship_quick.${sample}.prune.in --pca 'header' \
    --out data/pca/vasaship_quick.${sample}
```
