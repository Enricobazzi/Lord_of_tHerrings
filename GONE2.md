# Run demographic reconstructions with GONE2

To run [GONE2](https://github.com/esrud/GONE2) ([Santiago et al., 2025](https://doi.org/10.1038/s41467-025-61378-w)) I need the following:

- gone2 installed
- input ped map files for each population
- input parameters

## GONE2 installation

Downloaded [latest release](https://github.com/esrud/GONE2/releases/tag/v1.0.2) which contains the `gone2` executable.

## Input PED/MAP files

To create a PED and MAP file for each population to be analyzed using plink, I first need to merge and filter the individual BCFs of the individuals of that population.

The following populations will be analyzed:
- risor_2003 (n=8): MHER009 MHER010 MHER011 MHER012 MHER038 MHER039 MHER065 MHER066
- ...

Files with list of individuals for each populations are in `data/GONE2/<population>.sample_list.txt`

### Merge individual BCFs into population BCF

```
pop=risor_2003

for sample in $(cat data/GONE2/${pop}.sample_list.txt); do
    echo data/bcfs/${sample}.step3.bcf
done > data/GONE2/${pop}.bcf_list.txt

sbatch \
    --job-name=${pop}_merge \
    --output=logs/GONE2/${pop}_merge.out \
    --error=logs/GONE2/${pop}_merge.err \
    -t 0-08:00:00 \
    src/calling_and_filtering/merge_bcfs.sh \
        data/GONE2/${pop}.bcf_list.txt \
        data/GONE2/${pop}.bcf
```

### Get SNPs in chromosomes, remove high missing data and recode to PED/MAP

```
bcftools view -Ou \
    -m2 -M2 -v snps data/GONE2/${pop}.bcf | \
bcftools view -Ou \
    -t $(echo $(bcftools index -s data/GONE2/${pop}.bcf | cut -f1 | grep CM0) | tr ' ' ',') | \
bcftools view -Ou \
    -i 'F_MISSING <= 0.15' \
    -o data/GONE2/${pop}.lowmiss_snps.bcf

plink \
 --bcf data/GONE2/${pop}.lowmiss_snps.bcf \
 --double-id --allow-extra-chr --set-missing-var-ids @:# \
 --indep-pairwise 50 10 0.5 \
 --out data/GONE2/${pop}.lowmiss_snps

plink --recode \
 --bcf data/GONE2/${pop}.lowmiss_snps.bcf \
 --double-id --allow-extra-chr --set-missing-var-ids @:# \
 --extract data/GONE2/${pop}.lowmiss_snps.prune.in \
 --out data/GONE2/${pop}.lowmiss_snps
```

## Run GONE2

```
i=1
src/GONE2/gone2 \
    -g 0 \
    -r 2.54 \
    -s 200000 \
    -o data/GONE2/${pop}.${i} \
    data/GONE2/${pop}.lowmiss_snps.ped
```