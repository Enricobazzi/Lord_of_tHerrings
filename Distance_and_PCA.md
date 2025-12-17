## Angsd makematrix

alternative run Angsd command to obtain covariance/ibs matrix (see [here](https://www.popgen.dk/angsd/index.php/PCA_MDS)):
```
# example code from Morgan
inds=68
angsd -bam ${bams} -rf ${autosomes} -minMapQ 30 -minQ 30 -GL 1 -doMajorMinor 1 -doMaf 1 -SNP_pval 1e-6 -doIBS 1 -doCounts 1 -doCov 1 -minInd ${inds} -makeMatrix 1 -minMaf 0.05 -P 20 -remove_bads 1 -only_proper_pairs 1 -uniqueOnly 1 -skipTriallelic 1 -out ${output}
```

### bamlists

#### initial test

test it out with the following individuals:
```
dataset=test

# idefjord 4
grep -i idefjord-outer data/samples_table.csv | \
    sort -t ',' -nrk 10 | cut -d',' -f1 | head -4 >> data/angsd_matrix/bamlists/${dataset}.sample_list.txt
# maseskar 4
grep -i maseskar data/samples_table.csv | \
    sort -t ',' -nrk 10 | cut -d',' -f1 | head -4 >> data/angsd_matrix/bamlists/${dataset}.sample_list.txt
# north sea 4
grep -i northsea data/samples_table.csv | grep -i kong | \
    sort -t ',' -nrk 10 | cut -d',' -f1 | head -4 >> data/angsd_matrix/bamlists/${dataset}.sample_list.txt
# norwegian sea 4
grep -i norwegian data/samples_table.csv | grep -i kong | \
    sort -t ',' -nrk 10 | cut -d',' -f1 | head -4 >> data/angsd_matrix/bamlists/${dataset}.sample_list.txt
# faroe 4
grep -i faroe data/samples_table.csv | grep -i kong | \
    sort -t ',' -nrk 10 | cut -d',' -f1 | head -4 >> data/angsd_matrix/bamlists/${dataset}.sample_list.txt
# iceland 4 
grep -i iceland data/samples_table.csv | grep -i kong | \
    sort -t ',' -nrk 10 | cut -d',' -f1 | head -4 >> data/angsd_matrix/bamlists/${dataset}.sample_list.txt
# uk 4
grep -iE "downs|celtic" data/samples_table.csv | \
    sort -t ',' -nrk 10 | cut -d',' -f1 | head -4 >> data/angsd_matrix/bamlists/${dataset}.sample_list.txt

# could add MODERN: risor & norway fjords

# historical transition (no kampinge_1300 no knastorp_600)
grep ND data/samples_table.csv | grep -v NA | grep TRANS | \
    grep -vE ",600,|,1300," | cut -d',' -f1 >> data/angsd_matrix/bamlists/${dataset}.sample_list.txt
```
#### all wp1 samples

select samples to include:
```
dataset=wp1_all

grep "harengus" data/samples_table.csv | grep -v NA | grep -vE "BALTIC|NW-ATL|Lamichh|MartinezBarrio" | cut -d',' -f1 > data/angsd_matrix/bamlists/${dataset}.sample_list.txt
```

#### subset to test out snps

```
dataset=wp1_subset

for pop in faroe norwegian northsea iceland celtic downs karmoy risor maseskar idefjord isleofman; do
    grep "harengus" data/samples_table.csv | grep -v NA | grep -vE "BALTIC|NW-ATL|Lamichh|MartinezBarrio" | \
        grep $pop | head -10 | cut -d',' -f1 >> data/angsd_matrix/bamlists/${dataset}.sample_list.txt
done

grep "harengus" data/samples_table.csv | grep -v NA | grep -vE "BALTIC|NW-ATL|Lamichh|MartinezBarrio" | \
    grep -E "more_1|more_3|more_4|netherlands_7|york_1|york_5|selso_9|york_9|lyminge_3|lyminge_4|ND" | \
    cut -d',' -f1 >> data/angsd_matrix/bamlists/${dataset}.sample_list.txt
```

#### create bamlist

create bamlist from sample file:
```
dataset=test
dataset=wp1_all
dataset=wp1_subset

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

### SNP sites

Get 794 SNPs under selection to use as `-sites` from [Supplementary file 7](https://elifesciences.org/articles/61076/figures#files), of supplementary material of [Han et al. 2020](https://elifesciences.org/articles/61076)

```
cut -d',' -f1,2 data/angsd_matrix/sites/supplementary_file_7.csv | \
    tr ',' ' ' | grep -v "CHR" | awk '{print $1, $2-1, $2}' | tr ' ' '\t' \
    > data/angsd_matrix/sites/supplementary_file_7.chrn.bed

python src/liftover/change_chrn_to_v2.py \
    --map_file data/liftover/v3_v2_chrn_chromosome_names.txt \
    --input_file data/angsd_matrix/sites/supplementary_file_7.chrn.bed \
    --output_file data/angsd_matrix/sites/supplementary_file_7.v2.bed

cut -f1,3 data/angsd_matrix/sites/supplementary_file_7.v2.bed \
    > data/angsd_matrix/sites/supplementary_file_7.v2.sites

angsd sites index data/angsd_matrix/sites/supplementary_file_7.v2.sites
```

### run angsd -doIBS (individual base sampling)

Individual base sampling avoids snp/gt-likelihood calling.

```
# angsd -bam ${bams} -rf ${autosomes} -minMapQ 30 -minQ 30 -GL 1 -doMajorMinor 1 -doMaf 1 -SNP_pval 1e-6 -doIBS 1 -doCounts 1 -doCov 1 -minInd ${inds} -makeMatrix 1 -minMaf 0.05 -P 20 -remove_bads 1 -only_proper_pairs 1 -uniqueOnly 1 -skipTriallelic 1 -out ${output}

dataset=test
dataset=wp1_subset

sites=data/angsd_matrix/sites/supplementary_file_7.v2.sites

# interactively
angsd \
    -bam data/angsd_matrix/bamlists/${dataset}.bamlist \
    -sites ${sites} \
    -out data/angsd_matrix/output/${dataset} \
    -minMapQ 30 -minQ 20 -GL 1 -doMajorMinor 1 -doMaf 1 -SNP_pval 1e-6 \
    -doIBS 1 -doCounts 1 -doCov 1 -makeMatrix 1 -minMaf 0.001 \
    -remove_bads 1 -only_proper_pairs 1 -uniqueOnly 1 -skipTriallelic 1

# sbatch
sbatch \
    --job-name=${dataset}.doibs \
    --output=logs/angsd_matrix/doibs.${dataset}.out \
    --error=logs/angsd_matrix/doibs.${dataset}.err \
    src/angsd_matrix/angsd_doibs_makematrix.sh ${dataset} ${sites}
```

## PCAngsd

see [PCangsd](https://www.popgen.dk/software/index.php/PCAngsd)

### calculate gt-likelihoods with angsd

```
dataset=wp1_all

# sbatch
sbatch \
    --job-name=${dataset}.gtlike \
    --output=logs/angsd_matrix/gtlike.${dataset}.out \
    --error=logs/angsd_matrix/gtlike.${dataset}.err \
    src/angsd_matrix/calculate_gtlike_angsd.sh ${dataset}

# or interactively:
angsd -GL 2 -nThreads 1 -doGlf 2 -doMajorMinor 1 -doMaf 2 -SNP_pval 1e-6 \
    -minMapQ 30 -minQ 20 -remove_bads 1 -only_proper_pairs 1 -uniqueOnly 1 -skipTriallelic 1 \
    -bam data/angsd_matrix/bamlists/${dataset}.bamlist \
    -sites data/angsd_matrix/sites/supplementary_file_7.v2.sites \
    -out data/angsd_matrix/gtlike/${dataset}.supplementary_file_7_sites
```

### run pcangsd

filter the beagle file to only include sites from sites file:
```
zcat data/angsd_matrix/gtlike/wp1_all.beagle.tmp.gz \
| awk '
NR==FNR { sites[$1 "_" $2]=1; next }
NR==1 || ($1 in sites)
' data/angsd_matrix/sites/supplementary_file_7.v2.sites - \
| gzip > data/angsd_matrix/gtlike/wp1_all.filtered.beagle.tmp.gz
```

generate a file for filtering desired samples:
```
dataset=wp1_subset

awk 'NR==FNR{a[$0];next}{print ($0 in a)?1:0}' \
    data/angsd_matrix/bamlists/${dataset}.sample_list.txt \
    data/angsd_matrix/bamlists/wp1_all.sample_list.txt \
    > data/angsd_matrix/bamlists/${dataset}.samplemask.txt
```

run pcangsd:
```
ml pcangsd

dataset=wp1_subset

pcangsd \
    -b data/angsd_matrix/gtlike/wp1_all.filtered.beagle.tmp.gz \
    --iter 10000 \
    --filter data/angsd_matrix/bamlists/${dataset}.samplemask.txt \
    -o ${dataset}.filtered.pcangsd
```