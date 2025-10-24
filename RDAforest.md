# Run RDAforest!

I start with modern dataset of NE-ATLANTIC, TRANSITION and BALTIC (see [this script](src/maps/modern_samples_ne-bal-trans.R) for sample selection and map). List of samples is in `data/rdaforest/modern_samples.txt`

## obtain dataset

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

bcftools view data/modern_samples.step3.bcf > aaa.vcf

sed -i '$d' aaa.vcf

bcftools view -m2 -M2 -v snps aaa.vcf | bcftools +fill-tags -- -t AC,AN,AF | bcftools annotate -x ^INFO/AC,^INFO/AN,^INFO/AF > bbb.vcf

vcftools --vcf bbb.vcf --missing-indv
bcftools view -S <(awk '$5 < 0.5' out.imiss | cut -f1) -i 'F_MISS = 0' bbb.vcf

bcftools view -S <(awk '$5 < 0.2' out.imiss | cut -f1) bbb.vcf | bcftools view -i 'F_MISSING = 0 & MAF > 0.05' > ccc.vcf

ml plink
plink --vcf ccc.vcf \
 --double-id --allow-extra-chr --set-missing-var-ids @:# \
 --recode A --out ccc

```