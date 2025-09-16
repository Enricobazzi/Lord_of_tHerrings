## Get mapstats for screening data

I ran alignments on data that needed screening. To obtain mapping statistics from them I run Nic's bash script [collect_ancient_seqstats.sh](src/screening/collect_ancient_seqstats.sh):
```
for sample in $(cat config/historical_samples_paths.txt | cut -f1 | cut -d'_' -f1); do
    sbatch collect_ancient_seqstats.sh $sample
done

echo "sample total_seqs total_merged total_mapped total_uniq total_MQ25 proportion_merged endogenous merged_endogenous complexity grr total_cov read_length" \
    > seqstats/summary_screening.txt

cat seqstats/*_ancient_screening_stats.txt >> seqstats/summary_screening.txt
```

## Analyze screening data together with poolseq data

After lifting over the SNPs datasets to translate coordinates between version 2 and 3 of the assembly (see [Liftover.md](./Liftover.md))
```
sed 's/\t/./' data/published_data/60.Neff.freq > tests/all_snps.matrix

for snp in $(cut -f4 Baltic_v_Atlantic_snps.v3.bed); do
    grep -wm1 $snp all_snps.matrix >> baltic_v_atlantic_snps.matrix
done

for snp in $(cut -f4 Spring_v_Autumn_snps.v3.bed); do
    grep -wm1 $snp all_snps.matrix >> spring_v_autumn_snps.matrix
done
```

```
for sample in ND375 ND374 ND364 ND365; do
    echo $sample
    bedtools intersect -loj \
        -a ../../liftover_herring/Baltic_v_Atlantic_snps.v3.bed \
        -b <(bcftools view results/historical/vcf/GCA_040183275.1_Ch_v3.0_genomic/${sample}.merged.rmdup.merged.realn.rescaled.Q30.sorted.bcf) \
        | cut -f4,14 | cut -d':' -f1 | sed 's/0\/0/0.0/' | sed 's/0\/1/0.5/' | sed 's/1\/1/1.0/' \
        > ../../liftover_herring/${sample}.baltic_v_atlantic_snps.afgt
done

for sample in ND375 ND374 ND364 ND365; do
    echo $sample
    bedtools intersect -loj \
        -a ../../liftover_herring/Spring_v_Autumn_snps.v3.bed \
        -b <(bcftools view results/historical/vcf/GCA_040183275.1_Ch_v3.0_genomic/${sample}.merged.rmdup.merged.realn.rescaled.Q30.sorted.bcf) \
        | cut -f4,14 | cut -d':' -f1 | sed 's/0\/0/0.0/' | sed 's/0\/1/0.5/' | sed 's/1\/1/1.0/' \
        > ../../liftover_herring/${sample}.spring_v_autumn_snps.afgt
done

for sample in ND375 ND374 ND364 ND365; do
    echo $sample
    bedtools intersect -loj \
        -a ../../liftover_herring/All_snps.v3.bed \
        -b <(bcftools view results/historical/vcf/GCA_040183275.1_Ch_v3.0_genomic/${sample}.merged.rmdup.merged.realn.rescaled.Q30.sorted.bcf) \
        | cut -f4,14 | cut -d':' -f1 | sed 's/0\/0/0.0/' | sed 's/0\/1/0.5/' | sed 's/1\/1/1.0/' \
        > ../../liftover_herring/${sample}.all_snps.afgt
done
```

## Get stats for Herring2 sequencing

```
samples=($(basename -a results/historical/mapping/GCA_040183275.1_Ch_v3.0_genomic/*.rescaled.bam | cut -d'.' -f1))
for sample in ${samples[*]}; do
    echo $sample
    sbatch \
        --job-name=${sample}_stats \
        --output=logs/seqstats/${sample}_stats.out \
        --error=logs/seqstats/${sample}_stats.err \
        get_stats.herring_2.sh ${sample}
done

echo sample total_seqs total_merged total_mapped total_uniq total_MQ25 proportion_merged endogenous merged_endogenous complexity grr total_cov > seqstats/seqstats_all.txt
cat seqstats/*_seqstats.txt >> seqstats/seqstats_all.txt

```