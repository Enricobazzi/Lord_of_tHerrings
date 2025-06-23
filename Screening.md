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

