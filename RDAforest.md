# Run RDAforest!

I start with modern dataset of NE-ATLANTIC, TRANSITION and BALTIC (see [this script](src/maps/modern_samples_ne-bal-trans.R) for sample selection and map). List of samples is in `data/rdaforest/modern_samples.txt`

## obtain dataset - hard calls

### get bcf of modern samples

```
# on daredel
cd /cfs/klemming/scratch/e/ebazzica/rdaforest

BCFS=/cfs/klemming/projects/supr/naiss2024-6-170/analyses/Bcfs_Herring
for sample in $(cat data/modern_samples.txt); do
    echo ${BCFS}/${sample}.step3.bcf >> data/modern_bcfs.txt
done

ml bcftools

bcftools merge -m snps -Ob \
    -o data/modern_samples.step3.bcf \
    --file-list data/modern_bcfs.txt

bcftools view -Ou -m2 -M2 -v snps data/modern_samples.step3.bcf | bcftools +fill-tags -Ou -- -t AC,AN,AF | bcftools annotate -Ov -x ^INFO/AC,^INFO/AN,^INFO/AF > data/modern_samples.step3.biall_acanaf.vcf

vcftools --vcf data/modern_samples.step3.biall_acanaf.vcf --missing-indv

bcftools view -S <(awk '$5 < 0.2' out.imiss | cut -f1 | grep -v "MB") data/modern_samples.step3.biall_acanaf.vcf | bcftools view -i 'F_MISSING = 0 & MAF > 0.05' > data/modern_samples.step3.biall_acanaf.02missmax_0fmiss_005maf.vcf

ml plink
plink --vcf data/modern_samples.step3.biall_acanaf.02missmax_0fmiss_005maf.vcf \
 --double-id --allow-extra-chr --set-missing-var-ids @:#\$1 \
 --recode A --out data/modern_samples.step3.biall_acanaf.02missmax_0fmiss_005maf
```

### get historical raw matrix entry

test with some samples
```
# on daredel
cd /cfs/klemming/scratch/e/ebazzica/rdaforest

# get bed of snps
grep -v "#" data/modern_samples.step3.biall_acanaf.02missmax_0fmiss_005maf.vcf | \
    awk '{print $1, $2-1, $2}' | tr ' ' '\t' \
    > data/02missmax_0fmiss_005maf.snps.bed

ml plink
ml bcftools

# loop through samples:
bed=data/02missmax_0fmiss_005maf.snps.bed
BCFS=/cfs/klemming/projects/supr/naiss2024-6-170/analyses/Bcfs_Herring
samples=$(ls ${BCFS}/*.step3.bcf | rev | cut -d'/' -f1 | rev | grep -E "HER" | grep -vE "MHER" | cut -d'.' -f1)
for sample in ${samples[@]}; do
    echo $sample
    bcf=${BCFS}/${sample}.step3.bcf
    bcftools view -Ov -T ${bed} ${bcf} > data/${sample}.02missmax_0fmiss_005maf.snps.vcf
    plink --vcf data/${sample}.02missmax_0fmiss_005maf.snps.vcf \
     --double-id --allow-extra-chr --set-missing-var-ids @:# \
     --recode A --out data/${sample}.02missmax_0fmiss_005maf.snps
done
```

## obtain dataset - distance from bam files:

see:
https://academic.oup.com/gigascience/article/doi/10.1093/gigascience/giaf032/8106438
https://github.com/gkanogiannis/fastreeR

### distAngsd
from the 209 modern samples + 124 our historical + 73 atmore historical = total of 406 individuals:
that makes 82215 pairs! see if feasible to calculate pairwise like this or there is more efficient way?

see [distAngsd.md](./distAngsd.md)

### Angsd makematrix
alternative run Angsd command to obtain covariance/ibs matrix (see [here](https://www.popgen.dk/angsd/index.php/PCA_MDS)):
```
# example code from Morgan
inds=68
angsd -bam ${bams} -rf ${autosomes} -minMapQ 30 -minQ 30 -GL 1 -doMajorMinor 1 -doMaf 1 -SNP_pval 1e-6 -doIBS 1 -doCounts 1 -doCov 1 -minInd ${inds} -makeMatrix 1 -minMaf 0.05 -P 20 -remove_bads 1 -only_proper_pairs 1 -uniqueOnly 1 -skipTriallelic 1 -out ${output}
```