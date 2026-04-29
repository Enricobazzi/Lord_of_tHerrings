## Use ANGSD to calculate Diversity metrix

pipeline from Maria:

https://github.com/mlucenaperez/ll_phylogeography/blob/master/3.per_pop_analisis/1.diversity.Rmd

### CREATE FOLDERS

```
mkdir data/diversity
mkdir data/diversity/bams
mkdir data/diversity/bamlists
mkdir data/diversity/beds
mkdir data/diversity/output
mkdir logs/diversity
mkdir src/diversity
```

###  BED OF GENOME WITHOUT INVERSIONS AND REPEATS

```
ml bedtools

awk '{print $1, 0, $2}' Reference/GCF_900700415.2_Ch_v2.0.2_genomic.fna.fai | tr ' ' '\t' \
    > Reference/GCF_900700415.2_Ch_v2.0.2_genomic.whole_genome.bed

bedtools merge -i <(cat Reference/GCF_900700415.2_Ch_v2.0.2_genomic.repeats.sorted.bed \
    data/angsd_matrix/sites/ns_inversions.bed | \
    sort -k1,1 -k2,2n -k3,3n) \
    > data/diversity/beds/reps_invs.bed

bedtools subtract \
    -a Reference/GCF_900700415.2_Ch_v2.0.2_genomic.whole_genome.bed \
    -b data/diversity/beds/reps_invs.bed \
    > data/diversity/beds/noreps_noinvs.bed
```

###  FILTER BAMS TO REMOVE INVERSIONS AND REPEATS

```
for sample in $(cat data/angsd_matrix/bamlists/full_herr.sample_list.txt); do
    echo "${sample}"
    sbatch \
        --job-name=${sample}.norepnoinv \
        --output=logs/diversity/norepnoinv.${sample}.out \
        --error=logs/diversity/norepnoinv.${sample}.err \
        src/diversity/bam_filter_reps_invs.sh ${sample}
done
```

### POPULATIONS SAMPLE LISTS

```
#### SKAGERRAK & KATTEGAT (SK)

# 1700s sill period (sk-17sp)
grep -i "sillperioder" data/samples_table.csv | grep -iE "masthugget|gullholmen" | \
    cut -d',' -f1 > data/diversity/bamlists/sk-17sp.sample_list.txt
# 1800s regular herring (sk-18rh)
grep -i "historical" data/samples_table.csv | grep -iE "dynekilen" | \
    cut -d',' -f1 > data/diversity/bamlists/sk-18rh.sample_list.txt
# 1800s sill period (sk-18sp)
grep -i "sillperioder" data/samples_table.csv | grep -iE "koster|kalvsund" | \
    cut -d',' -f1 > data/diversity/bamlists/sk-18sp.sample_list.txt
# modern herring (sk-mh)
grep -i "modern" data/samples_table.csv | grep -iE "idefjord|maseskar|risor" | \
    cut -d',' -f1 > data/diversity/bamlists/sk-mh.sample_list.txt

#### NORTH SEA (NS)

# ancient (pre-sillperiods) ns-ah
grep -i "historical" data/samples_table.csv | grep -iE "york|lyminge|netherlands" | \
    awk -F',' '$9 == "current" || $10 >= 0.1' | grep -v "HER135" | \
    cut -d',' -f1 > data/diversity/bamlists/ns-ah.sample_list.txt
# 1800s regular herring (ns-18rh)
grep -i "historical" data/samples_table.csv | grep -iE "stavanger|haugesund|foldfjorden|scotland" | \
    awk -F',' '$9 == "current" || $10 >= 0.1' | grep -v "HER135" | \
    cut -d',' -f1 > data/diversity/bamlists/ns-18rh.sample_list.txt
# 1800s sill period (ns-18sp)
grep -i "sillperioder" data/samples_table.csv | grep -iE "stavanger|haugesund|rover|scotland|norway_unknown" | \
    awk -F',' '$9 == "current" || $10 >= 0.1' | grep -v "HER135" | \
    cut -d',' -f1 > data/diversity/bamlists/ns-18sp.sample_list.txt
# modern herring (ns-mh)
grep -i "modern" data/samples_table.csv | grep -iE "northsea|celtic|downs|isleofman|karmoy" | \
    grep -v "Lamich" | \
    cut -d',' -f1 > data/diversity/bamlists/ns-mh.sample_list.txt

```

### GENERATE BAMLISTS

```
for dataset in sk-17sp sk-18rh sk-18sp sk-mh ns-ah ns-18rh ns-18sp ns-mh; do
    for sample in $(cat data/diversity/bamlists/${dataset}.sample_list.txt); do
        input_bam=data/diversity/bams/${sample}.subsampled_3X.noreps_noinvs.bam
        echo ${input_bam}
    done > data/diversity/bamlists/${dataset}.bamlist
done
```

### ANGSD DO SAF (LIKELIHOOD)

```
for dataset in sk-17sp sk-18rh sk-18sp sk-mh ns-ah ns-18rh ns-18sp ns-mh; do
    sbatch \
        --job-name=${dataset}.dosaflike \
        --output=logs/diversity/dosaflike.${dataset}.out \
        --error=logs/diversity/dosaflike.${dataset}.err \
        src/diversity/angsd_dosaf_likelihood.sh ${dataset}
done
```

### REALSFS

```
for dataset in sk-17sp sk-18rh sk-18sp sk-mh ns-ah ns-18rh ns-18sp ns-mh; do
    sbatch \
        --job-name=${dataset}.realsfs \
        --output=logs/diversity/realsfs.${dataset}.out \
        --error=logs/diversity/realsfs.${dataset}.err \
        src/diversity/run_realsfs.sh ${dataset}
done
```

### THETASTAT

```
for dataset in sk-17sp sk-18rh sk-18sp sk-mh ns-ah ns-18rh ns-18sp ns-mh; do
    sbatch \
        --job-name=${dataset}.dosafpest \
        --output=logs/diversity/dosafpest.${dataset}.out \
        --error=logs/diversity/dosafpest.${dataset}.err \
        src/diversity/run_thetastat.sh ${dataset}
done
```

## ALTERNATIVE - CALCULATE INDIVIDUAL HET

### RUN SAF AND REALSFS:

```
for sample in $(cat data/angsd_matrix/bamlists/full_herr.sample_list.txt); do
    echo "${sample}"
    sbatch \
        --job-name=${sample}.indhet \
        --output=logs/diversity/indhet.${sample}.out \
        --error=logs/diversity/indhet.${sample}.err \
        src/diversity/angsd_saf_het.sh ${sample}
done
```

### MANUALLY FROM BCF:

```
for sample in $(cat data/angsd_matrix/bamlists/full_herr.sample_list.txt); do
    echo "${sample}"
    sbatch \
        --job-name=${sample}.bcfhet \
        --output=logs/diversity/bcfhet.${sample}.out \
        --error=logs/diversity/bcfhet.${sample}.err \
        src/diversity/het_from_bcf.sh ${sample}
done

# This is shit because heterozygosity and depth(number of sites) is highly correlated
```

## CALCULATE DIVERSITY ONLY OF SUBSET OF SNPS

### CREATE BED FILES OF THETA VALUES

```
for dataset in sk-17sp sk-18rh sk-18sp sk-mh ns-ah ns-18rh ns-18sp ns-mh; do
    echo "${dataset}"
    sbatch \
        --job-name=${dataset}.thetabed \
        --output=logs/diversity/thetabed.${dataset}.out \
        --error=logs/diversity/thetabed.${dataset}.err \
        src/diversity/get_thetas_bed.sh ${dataset}
done
```

### INTERSECT WITH BED OF SNPS

```
ml bedtools

for dataset in sk-17sp sk-18rh sk-18sp sk-mh ns-ah ns-18rh ns-18sp ns-mh; do
    echo "${dataset}"
    bedtools intersect -wa \
        -a data/diversity/output/${dataset}.folded.thetas.bed \
        -b data/angsd_matrix/sites/supplementary_file_7.v2.bed \
        > data/diversity/output/${dataset}.supplementary_file_7.v2.bed
done
```