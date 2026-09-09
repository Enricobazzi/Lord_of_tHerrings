# testing how to run MlRho

run script on samples:
```
for sample in $(cat data/angsd_matrix/bamlists/full_herr.sample_list.txt); do
    echo "${sample}"
    sbatch \
        --job-name=${sample}.mlrho \
        --output=logs/mlrho/${sample}.out \
        --error=logs/mlrho/${sample}.err \
        src/mlrho/run_mlrho_from_bam.sh ${sample}
done
```
