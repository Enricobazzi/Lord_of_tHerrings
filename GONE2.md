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
- risor_2008 (n=8): MHER009, MHER010, MHER011, MHER012, MHER038, MHER039, MHER065, MHER066
- maseskar_2003 (n=7): MHER013, MHER014, MHER015, MHER034, MHER035, MHER036, MHER037
- idefjord_2010 (n=13): MHER022, MHER023, MHER024, MHER028, MHER029, MHER030, MHER044, MHER045, MHER046, MHER061, MHER062, MHER063, MHER064
- outeridefjord_2010 (n=8): MHER028, MHER029, MHER030, MHER044, MHER045, MHER046, MHER063, MHER064
- northsea_2016 (n=17): 201650792, 201650793, 201650794, 201650795, 201650796, 201650797, 201650798, 201650799, 2016507910, 2016507911, 2016507912, 2016507913, 2016507915, 2016507916, 2016507918, 2016507919, 2016507920
- norwegiansea_2015 (n=29): 201550562, 201550563, 201550564, 201550565, 201550566, 201550567, 201550568, 201550569, 1552004523, 1552004524, 1552004526, 1552004722, 1552004723, 1552004724, 1552005123, 1552005124, 1552005125, 1552005524, 1552005525, 1552005923, 1552005924, 1552006326, 1552006327, 2015505610, 2015505612, 2015505613, 2015505622, 2015505623, 2015505624
- faroe_2016 (n=27): 201450544, 201450545, 201550141, 201550142, 201550143, 201550145, 201550147, 201550149, 201604591, 201750143, 201750144, 201750145, 201750146, 201750147, 201750149, 201750153, 201750154, 201750155, 201750157, 201750158, 201750159, 2015501410, 2015501411, 2017501410, 2017501411, 2017501414, 2017501510
- dynekilen_1874 (n=12): ND195, ND196, ND197, ND198, ND199, ND200, ND201, ND202, ND203, ND204, ND205, ND206
- masthugget_1740 (n=8): ND443, ND448, ND450, ND454, ND455, ND456, ND464, ND467
- kampinge_1300 (n=17): ND391, ND392, ND393, ND394, ND395, ND396, ND397, ND398, ND399, ND400, ND401, ND402, ND403, ND404, ND405, ND407, ND408
- knastorp_600 (n=7): ND409, ND410, ND411, ND414, ND415, ND417, ND418
- foldfjorden_1875 (n=15): ND324, ND325, ND326, ND328, ND329, ND330, ND331, ND332, ND333, ND334, ND335, ND336, ND337, ND342, ND343
- stavanger_1880 (n=10): ND178, ND180, ND182, ND183, ND184, ND187, ND188, ND319, ND320, ND321

Files with list of individuals for each populations are in `data/GONE2/${pop}.sample_list.txt`

### Merge individual BCFs into population BCF

```
pop=risor_2008
pop=maseskar_2003
pop=idefjord_2010
pop=northsea_2016
pop=norwegiansea_2015
pop=faroe_2016
pop=dynekilen_1874
pop=masthugget_1740
pop=kampinge_1300
pop=knastorp_600
pop=foldfjorden_1875
pop=stavanger_1880

for sample in $(cat data/GONE2/${pop}.sample_list.txt); do
    echo data/bcfs/${sample}.step3.bcf
done > data/GONE2/${pop}.bcf_list.txt

sbatch \
    --job-name=${pop}_merge \
    --output=logs/GONE2/${pop}_merge.out \
    --error=logs/GONE2/${pop}_merge.err \
    -t 0-18:00:00 \
    src/calling_and_filtering/merge_bcfs.sh \
        data/GONE2/${pop}.bcf_list.txt \
        data/GONE2/${pop}.bcf
```

### Get SNPs in chromosomes, remove high missing data and recode to PED/MAP

```
for pop in idefjord_2010 faroe_2016; do
    echo ${pop}
done

bcftools view -Ou \
    -m2 -M2 -v snps data/GONE2/${pop}.bcf | \
bcftools view -Ou \
    -t $(echo $(bcftools index -s data/GONE2/${pop}.bcf | cut -f1 | grep CM0) | tr ' ' ',') | \
bcftools view -Ou \
    -i 'F_MISSING <= 0.15 & MAF > 0.05' \
    -o data/GONE2/${pop}.lowmiss_snps.bcf


# vcftools --bcf data/GONE2/${pop}.bcf --missing-indv
bcftools view -Ou \
    -m2 -M2 -v snps data/GONE2/${pop}.bcf | \
bcftools view -Ou \
    -S <(awk '$5 < 0.6' out.imiss | cut -f1) | \
bcftools view -Ou \
    -t $(echo $(bcftools index -s data/GONE2/${pop}.bcf | cut -f1 | grep CM0) | tr ' ' ',') | \
bcftools view -Ou \
    -i 'F_MISSING <= 0.15 & MAF > 0.05' \
    -o data/GONE2/${pop}.lowmiss_snps.bcf

# subsetting extract outerfoldjord
# pop=outeridefjord_2010
bcftools view -Ou \
    -S data/GONE2/${pop}.sample_list.txt \
    data/GONE2/idefjord_2010.bcf | \
bcftools view -Ou \
    -m2 -M2 -v snps | \
bcftools view -Ou \
    -t $(echo $(bcftools index -s data/GONE2/idefjord_2010.bcf | cut -f1 | grep CM0) | tr ' ' ',') | \
bcftools view -Ou \
    -i 'F_MISSING <= 0.15 & MAF > 0.05' \
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
pop=risor_2008
pop=maseskar_2003
pop=idefjord_2010
pop=outeridefjord_2010
pop=northsea_2016
pop=norwegiansea_2015
pop=faroe_2016
pop=dynekilen_1874
pop=masthugget_1740
pop=kampinge_1300
pop=knastorp_600
pop=foldfjorden_1875
pop=stavanger_1880

for pop in idefjord_2010 outeridefjord_2010; do
    echo ${pop}
    for i in {1..50}; do
        sbatch \
            --job-name=${i}.${pop}.gone2 \
            --output=logs/GONE2/${pop}.${i}_GONE2.out \
            --error=logs/GONE2/${pop}.${i}_GONE2.err \
            src/GONE2/run_gone.sh \
            ${pop} ${i} 0
    done
done
```