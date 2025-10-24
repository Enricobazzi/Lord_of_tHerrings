## Translate coordinates from Herring reference v2 to v3

### Download data:

Version 2 of herring reference genome:
[GCA_900700415](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_900700415.2/)
Version 3 of herring reference genome:
[GCA_040183275](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_040183275/)
```
ml bioinfo-tools NCBI-datasets/15.29.0
datasets download genome accession GCF_900700415.2 --include genome
```

### Run minimap2

```
ml minimap2/2.28
D_REF=/cfs/klemming/projects/supr/naiss2024-6-170/analyses/Reference
minimap2 -cx asm5 --cs \
    ${D_REF}/GCF_900700415.2_Ch_v2.0.2_genomic.fna \
    ${D_REF}/GCA_040183275.1_Ch_v3.0_genomic.fna \
    > v2_to_v3.paf
```

### Run paf2chain

/cfs/klemming/home/e/ebazzica/.cargo/bin/paf2chain -i v2_to_v3.paf > v2_to_v3.chain

### Convert coordinates of all and outlier snps

#### All SNPs

The allele frequencies matrix of all SNPs in all populations was downloaded from [Han et al. dryad database](https://doi.org/10.5061/dryad.pnvx0k6kr) in a file called `60.Neff.freq`.

To obtain a bed file of all the SNPs with the chromosome name as "chr#" (save original position as fourth column "chromosome.position"):
```
grep -v "CHROM" data/published_data/60.Neff.freq | cut -f1,2 | awk '{print $1, $2-1, $2, $1 "." $2}' | tr ' ' '\t' \
    > data/liftover/All_snps.chrn.bed
```
To change the chromosome name from chr# to v2 assembly version of the name:
```
python src/liftover/change_chrn_to_v2.py \
    --map_file data/liftover/v3_v2_chrn_chromosome_names.txt \
    --input_file data/liftover/All_snps.chrn.bed \
    --output_file data/liftover/All_snps.v2.bed
```

#### Outlier SNPs

Atlantic vs Baltic and Autumn vs Spring spawning significant SNPs were downloaded from Han et al. 2020 supplementary material (4 and 5) and saved as CSV (with a # added to the header for easy filtering).

To extract bed files with the chromosome name as "chr#" (save original position as fourth column "chromosome.position"):
```
grep -v "#" data/liftover/Baltic_v_Atlantic_snps.csv | cut -d',' -f1,2 | tr ',' ' ' \
    | awk '{print "chr" $1, $2-1, $2, "chr" $1 "." $2}' | tr ' ' '\t' \
    > data/liftover/Baltic_v_Atlantic_snps.chrn.bed
grep -v "#" data/liftover/Spring_v_Autumn_snps.csv | cut -d',' -f1,2 | tr ',' ' ' \
    | awk '{print "chr" $1, $2-1, $2, "chr" $1 "." $2}' | tr ' ' '\t' \
    > data/liftover/Spring_v_Autumn_snps.chrn.bed
```
To change the chromosome name from chr# to v2 assembly version of the name
```
python src/liftover/change_chrn_to_v2.py \
    --map_file data/liftover/v3_v2_chrn_chromosome_names.txt \
    --input_file data/liftover/Baltic_v_Atlantic_snps.chrn.bed \
    --output_file data/liftover/Baltic_v_Atlantic_snps.v2.bed

python src/liftover/change_chrn_to_v2.py \
    --map_file data/liftover/v3_v2_chrn_chromosome_names.txt \
    --input_file data/liftover/Spring_v_Autumn_snps.chrn.bed \
    --output_file data/liftover/Spring_v_Autumn_snps.v2.bed
```

### Run liftover to convert v2 bed files to v3

```
scp data/liftover/beds/Spring_v_Autumn_snps.v2.bed ebazzica@dardel.pdc.kth.se:/cfs/klemming/scratch/e/ebazzica/liftover/beds
scp data/liftover/beds/Baltic_v_Atlantic_snps.v2.bed ebazzica@dardel.pdc.kth.se:/cfs/klemming/scratch/e/ebazzica/liftover/beds
scp data/liftover/beds/All_snps.v2.bed ebazzica@dardel.pdc.kth.se:/cfs/klemming/scratch/e/ebazzica/liftover/beds

ml bioinfo-tools liftOver/2017-03-14

liftOver beds/All_snps.v2.bed v2_to_v3.chain beds/All_snps.v3.bed beds/All_snps.unmapped.bed
liftOver beds/Baltic_v_Atlantic_snps.v2.bed v2_to_v3.chain beds/Baltic_v_Atlantic_snps.v3.bed beds/Baltic_v_Atlantic_snps.unmapped.bed
liftOver beds/Spring_v_Autumn_snps.v2.bed v2_to_v3.chain beds/Spring_v_Autumn_snps.v3.bed beds/Spring_v_Autumn_snps.unmapped.bed

wc -l beds/All_snps.v2.bed beds/All_snps.v3.bed
# 15885982 All_snps.v2.bed
# 11635531 All_snps.v3.bed
wc -l beds/Baltic_v_Atlantic_snps.v2.bed beds/Baltic_v_Atlantic_snps.v3.bed
# 2292 Baltic_v_Atlantic_snps.v2.bed
# 2063 Baltic_v_Atlantic_snps.v3.bed
wc -l beds/Spring_v_Autumn_snps.v2.bed beds/Spring_v_Autumn_snps.v3.bed
# 1005 Spring_v_Autumn_snps.v2.bed
#  982 Spring_v_Autumn_snps.v3.bed
```

### bcfs of atlantic_v_autumn

```
samples=$(ls /cfs/klemming/projects/supr/naiss2024-6-170/analyses/Bcfs_Herring/*.step3.bcf | rev | cut -d'/' -f1 | rev | cut -d'.' -f1 | sort -u)
dataset=Baltic_v_Atlantic

for sample in ${samples[*]}; do
    sbatch \
        --job-name=${sample}.${dataset} \
        --output=logs/${sample}.${dataset}.out \
        --error=logs/${sample}.${dataset}.err \
        src/get_sample_dataset_bcf.sh \
        ${sample} ${dataset}
done

for sample in ND409 ND410 ND411 ND414 ND415 ND417 ND418 HER088 HER110 HER111 HER112 HER113 HER114 HER115 ND391 ND392 ND393 ND394 ND395 ND396 ND397 ND398 ND399 ND400 ND401 ND402 ND403 ND404 ND405 ND407 ND408 HER116 HER117 HER118 ND440 ND443 ND448 ND450 ND454 ND455 ND456 ND464 ND467 ND195 ND196 ND197 ND198 ND199 ND200 ND201 ND202 ND203 ND204 ND205 ND206 ND022 ND023 ND024 ND026 Fehmarn3 Fehmarn44 Fehmarn6 MHER009 MHER010 MHER011 MHER012 MHER013 MHER014 MHER015 MHER034 MHER035 MHER036 MHER037 MHER038 MHER039 MHER065 MHER066 MHER022 MHER023 MHER024 MHER028 MHER029 MHER030 MHER044 MHER045 MHER046 MHER061 MHER062 MHER063 MHER064 BF1 BF2 BF3 ND001 ND006 ND009 ND010 ND011; do
    echo bcfs/${sample}.${dataset}.bcf >> trans.bcflist
done
```