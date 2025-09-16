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

minimap2 -cx asm5 --cs \
    GCF_900700415.2_Ch_v2.0.2_genomic.fna \
    GCA_040183275.1_Ch_v3.0_genomic.fna \
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
scp data/liftover/Spring_v_Autumn_snps.v2.bed ebazzica@dardel.pdc.kth.se:/cfs/klemming/scratch/e/ebazzica/liftover_herring
scp data/liftover/Baltic_v_Atlantic_snps.v2.bed ebazzica@dardel.pdc.kth.se:/cfs/klemming/scratch/e/ebazzica/liftover_herring
scp data/liftover/All_snps.v2.bed ebazzica@dardel.pdc.kth.se:/cfs/klemming/scratch/e/ebazzica/liftover_herring

ml bioinfo-tools liftOver/2017-03-14

liftOver All_snps.v2.bed v2_to_v3.chain All_snps.v3.bed All_snps.unmapped.bed
liftOver Baltic_v_Atlantic_snps.v2.bed v2_to_v3.chain Baltic_v_Atlantic_snps.v3.bed Baltic_v_Atlantic_snps.unmapped.bed
liftOver Spring_v_Autumn_snps.v2.bed v2_to_v3.chain Spring_v_Autumn_snps.v3.bed Spring_v_Autumn_snps.unmapped.bed

wc -l All_snps.v2.bed All_snps.v3.bed
# 15885982 All_snps.v2.bed
# 11635531 All_snps.v3.bed
wc -l Baltic_v_Atlantic_snps.v2.bed Baltic_v_Atlantic_snps.v3.bed
# 2292 Baltic_v_Atlantic_snps.v2.bed
# 2063 Baltic_v_Atlantic_snps.v3.bed
wc -l Spring_v_Autumn_snps.v2.bed Spring_v_Autumn_snps.v3.bed
# 1005 Spring_v_Autumn_snps.v2.bed
#  982 Spring_v_Autumn_snps.v3.bed
```