# BayPass

## Install Baypass
```
git clone git@forge.inrae.fr:mathieu.gautier/baypass_public.git
unzip baypass_public-v3.1.zip
cd baypass_public-v3.1/sources
make clean all FC=gfortran
```
## Input preparation

### Allele Count file

generate 500 SNP subsets
```
awk 'NR>1 {
  n = ((NR-2) % 500) + 1         # which sequence (1..500)
  nn = sprintf("%03d", n)        # zero-padded filename part
  out = "data/selection_scans_poolseq/subsets_freqs/60.Neff." nn ".freq"
  print >> out
  close(out)
}' data/published_data/60.Neff.freq
```

convert allele frequency (AF) to allele count (AC):
```
python src/selection_scans_poolseq/generate_ac_files.py
```

### Environmental Data file

```
Rscript src/selection_scans_poolseq/get_env_data.R
```

### Population Structure Correction (omega-matrix) file

```
for n in {1..500}; do
    nn=$(printf "%03d\n" "$n")
    sbatch \
        --job-name=${nn}_core \
        --output=logs/core/${nn}.out \
        --error=logs/core/${nn}.err \
        baypass_core.sh ${nn}     
done
```

## Run!

### on indiviudal datasets

```
for n in {1..500}; do
    nn=$(printf "%03d\n" "$n")
    sbatch \
        --job-name=${nn}_aux \
        --output=logs/aux/${nn}.out \
        --error=logs/aux/${nn}.err \
        baypass_aux.sh ${nn}     
done
```

### join results

```
```

# RDA forest

## get environment associated snps

see [rda.R](src/selection_scans_poolseq/rda.R) to see how I got the candidate SNPs for subset 001

## get snps of sample *example with 001 subset*

### get bed of subset

```
awk '{print $1, $2 - 1, $2}' data/selection_scans_poolseq/subsets_acs/60.Neff.001.snps | tr ' ' '\t' > data/selection_scans_poolseq/subsets_beds/60.Neff.001.bed

python src/liftover/change_chrn_to_v2.py \
    --map_file data/liftover/v3_v2_chrn_chromosome_names.txt \
    --input_file data/selection_scans_poolseq/subsets_beds/60.Neff.001.bed \
    --output_file data/selection_scans_poolseq/subsets_beds/60.Neff.001.v2.bed
```

### liftover to v3

```
# on dardel:
cd /cfs/klemming/scratch/e/ebazzica/liftover

ml bioinfo-tools liftOver/2017-03-14
liftOver beds/60.Neff.001.v2.bed v2_to_v3.chain beds/60.Neff.001.v3.bed beds/60.Neff.001.unmapped.bed

wc -l beds/60.Neff.001.v2.bed beds/60.Neff.001.v3.bed
#  9584 beds/60.Neff.001.v2.bed
#  7653 beds/60.Neff.001.v3.bed

# to know which snp numbers are left
paste beds/60.Neff.001.v3.bed <(grep -vnf <(grep -v "#" beds/60.Neff.001.unmapped.bed) beds/60.Neff.001.v2.bed | cut -d':' -f1) > beds/60.Neff.001.v3.snpnumber.bed
```

### extract bcf of snps

```
# on dardel:
cd /cfs/klemming/scratch/e/ebazzica/rdaforest

ml bedtools
ml bcftools
sample=Z12
BCFS=/cfs/klemming/projects/supr/naiss2024-6-170/analyses/Bcfs_Herring
bed=/cfs/klemming/scratch/e/ebazzica/liftover/beds/60.Neff.001.v3.bed
inbcf=${BCFS}/${sample}.step3.bcf
outbcf=${sample}.subset_001.v3.bcf

bcftools view -Ob ${inbcf} -T ${bed} -o ${outbcf}
```

### transform to frequency + snp number table

```
# on dardel:
cd /cfs/klemming/scratch/e/ebazzica/rdaforest

ml bcftools
sample=Z12
bcf=${sample}.subset_001.v3.bcf
numfile=/cfs/klemming/scratch/e/ebazzica/liftover/beds/60.Neff.001.v3.snpnumber.bed

# frequency part
bcftools view $bcf | grep -v "#" | cut -f10 | cut -d':' -f1 | sed 's/0\/0/0.0/' | sed 's/0\/1/0.5/' | sed 's/1\/1/1.0/' > ${sample}.subset_001.fq.tmp

# snp number part
bedtools intersect -a ${numfile} \
    -b <(bcftools view $bcf | grep -v "#" | cut -f1,2 | awk '{print $1, $2 - 1, $2}' | tr ' ' '\t') | \
    cut -f4 > ${sample}.subset_001.snpnumber.tmp

# paste together
paste ${sample}.subset_001.snpnumber.tmp ${sample}.subset_001.fq.tmp > ${sample}.subset_001.freq
rm *.tmp
```
**DOES NOT WORK!**