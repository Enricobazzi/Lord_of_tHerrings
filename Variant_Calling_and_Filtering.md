# Variant Calling

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

# Filtering variants