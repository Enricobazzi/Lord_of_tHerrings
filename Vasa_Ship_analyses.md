# Analysis of Vasa ship herring

## Generate covariance and distance matrices with angsd

Create a sample and bam list files of the samples to include:
```
dataset=vasa_ship

# sample list
grep "harengus" data/samples_table.csv | \
    grep -v NA | \
    grep -vE "NW-ATL|Lamichh|MartinezBarrio|HER135" | \
    awk -F',' '$9 == "current" || $10 >= 0.1' | \
    cut -d',' -f1 \
    > data/angsd_matrix/bamlists/${dataset}.sample_list.txt

# bam list
for sample in $(cat data/angsd_matrix/bamlists/${dataset}.sample_list.txt); do
    if [[ -f data/bams/${sample}.merged.rmdup.merged.realn.rescaled.bam ]]; then
      input_bam=data/bams/${sample}.merged.rmdup.merged.realn.rescaled.bam
    elif [[ -f data/bams/${sample}.merged.rmdup.merged.realn.bam ]]; then
      input_bam=data/bams/${sample}.merged.rmdup.merged.realn.bam
    else
      echo "error: bam file for sample ${sample} not found!"
    fi
    echo ${input_bam}
done > data/angsd_matrix/bamlists/${dataset}.bamlist
```

run `angsd -doIBS` using [angsd_doibs_makematrix.sh](src/angsd_matrix/angsd_doibs_makematrix.sh)
```
dataset=vasa_ship
sites=data/angsd_matrix/sites/supplementary_file_7.v2.sites

# sbatch
sbatch \
    --job-name=${dataset}.doibs \
    --output=logs/angsd_matrix/doibs.${dataset}.out \
    --error=logs/angsd_matrix/doibs.${dataset}.err \
    src/angsd_matrix/angsd_doibs_makematrix.sh ${dataset} ${sites}
```