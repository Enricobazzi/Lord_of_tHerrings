```
for sample in $(cat data/bamlists/full_herr.sample_list.txt); do
    echo "${sample}"
    sbatch \
        --job-name=${sample}.mapdamage \
        --output=logs/mapdamage/mapdamage.${sample}.out \
        --error=logs/mapdamage/mapdamage.${sample}.err \
        src/mapdamage/run_mapdamage_plots.sh ${sample}
done
```