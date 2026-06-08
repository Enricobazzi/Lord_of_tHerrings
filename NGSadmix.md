# run ngsadmix

## test with chr1

```
ml bioinfo-tools NGSadmix/32
dataset=wp1_final_bal
k=4
i=1
THR=8

beaglefile=data/gtlike/${dataset}.beagle.gz
out=data/ngsadmix/${dataset}.admix.k_${k}.seed_${i}

time NGSadmix \
  -likes ${beaglefile} \
  -K ${k} \
  -P ${THR} \
  -seed ${i} \
  -o ${out}
```

## test with chr1 and sf7 sites only

```
# filter the beagle file to only include sites from sites file:
ml pigz

dataset=wp1_final_bal
sites=supplementary_file_7.v2

pigz -dc data/gtlike/${dataset}.beagle.gz \
| awk '
NR==FNR { sites[$1 "_" $2]=1; next }
NR==1 || ($1 in sites)
' data/sites/${sites}.sites - \
| pigz -p 8 -1 > data/gtlike/${dataset}.${sites}.beagle.gz

####

# then sbatch/run this:

ml bioinfo-tools NGSadmix/32

dataset=wp1_final_bal.NC_045152.1
sites=supplementary_file_7.v2

beaglefile=data/gtlike/${dataset}.${sites}.beagle.gz
k=4
THR=8
for i in {1..20}; do
    out=data/ngsadmix/${dataset}.${sites}.admix.k_${k}.seed_${i}
    NGSadmix \
      -likes ${beaglefile} \
      -K ${k} \
      -P ${THR} \
      -seed ${i} \
      -o ${out}
done

evalAdmix=src/ngsadmix/evalAdmix/evalAdmix

${evalAdmix} -beagle ${beaglefile} \
    -fname data/ngsadmix/${dataset}.${sites}.admix.k_${k}.seed_${i}.fopt.gz \
    -qname data/ngsadmix/${dataset}.${sites}.admix.k_${k}.seed_${i}.qopt \
    -o evaladmixOut.${dataset}.${sites}.admix.k_${k}.seed_${i}.corres -P ${THR}

```

## From Morgan:

```
#!/bin/bash
#SBATCH --job-name=ngsadmix # Job name
#SBATCH -e ngsadmix.error
#SBATCH --mail-type=ALL # Mail events (BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=morgan@bio.ku.dk # Where to send mail
#SBATCH --ntasks=1 # Run a single task
#SBATCH --cpus-per-task=20 # Number of CPU cores per task
#SBATCH --mem=20gb
#SBATCH --nodes=1
#SBATCH --time=72:00:00 # Time limit hrs:min:sec
#SBATCH --output=ngsadmix.log # Standard output/error
############################################

module load angsd

for k in `seq 18`
do
mkdir -p ${k}
cd ${k}
for i in `seq 1 100`
do
beagle=/projects/mjolnir1/people/nrb613/README_gray_seal_2022/google_drive/cleaned_reads/genotype_likelihoods/depth_of_3_wo_harbour_seals_outliers_50/gl.map.30_2e-6_postQC_loci.beagle.gz
beagle=/projects/mjolnir1/people/nrb613/README_gray_seal_2022/google_drive/cleaned_reads/genotype_likelihoods/depth_of_3_wo_harbour_seals_outliers_50/2024/postqc.gl.map.30.2024.beagle.gz
nfile=ngsadmix_result
NGSadmix -likes ${beagle} -K ${k} -o ${nfile}.${k}.${i} -P 20
done
cd ..
done
```

```
vim evalAdmix.sh
#!/bin/bash
#SBATCH --job-name=evalAdmix # Job name
#SBATCH -e evalAdmix.error
#SBATCH --mail-type=ALL # Mail events (BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=morgan@bio.ku.dk # Where to send mail
#SBATCH --ntasks=1 # Run a single task
#SBATCH --cpus-per-task=20 # Number of CPU cores per task
#SBATCH --mem=25gb
#SBATCH --nodes=1
#SBATCH --time=24:00:00 # Time limit hrs:min:sec
#SBATCH --output=evalAdmix.log # Standard output/error
############################################

for file in *fopt.gz
do
shortID=`echo ${file}| cut -f2,3 -d "."`

evalAdmix=/projects/mjolnir1/people/nrb613/X204SC22021885-Z01-F004/combined.nuclear/redo/ngsadmix/evalAdmix/evalAdmix
inputBeagleFile=/projects/mjolnir1/people/nrb613/X204SC22021885-Z01-F004/combined.nuclear/redo/snakemake.runs.final.global/postQC_genotype_likelihoods/EB.mapQ30.autosomes.postQC-p1.no.relatives.beagle.gz

${evalAdmix} -beagle ${inputBeagleFile} -fname ngsadmix_result.${shortID}.fopt.gz -qname ngsadmix_result.${shortID}.qopt -o evaladmixOut.${shortID}.corres -P 20

done
```

```
# https://raw.githubusercontent.com/GenisGE/evalAdmix/refs/heads/master/visFuns.R

source("/Users/morganmccarthy/Downloads/visFuns.R") 

cols_admix <- c("red","orange","yellow", "green", "blue", "purple","green","grey", "gold")

#Getting order info
ids <- read.table("/Users/morganmccarthy/Documents/Manuscripts/2025/WGS_blue_whales/ngsadmix/ids.txt")
ids <- ids %>%
  rename(Run=V1)

pop <- left_join(ids, B_musculus_meta, by = "Run") %>%
  pull(Abb2)

popord <- level_loc_order

ord <- orderInds(pop=pop,popord=ord)

ord <- orderInds(pop=pop, q=q, popord=popord)

#K1
r <- as.matrix(read.table("/Users/morganmccarthy/Documents/Manuscripts/2025/WGS_blue_whales/ngsadmix/evalAdmix/evaladmixOut.K_1.15.corres"), header=FALSE)
#pop <- bam.list.meta$AbbrevMap
q <- as.matrix(read.table("/Users/morganmccarthy/Documents/Manuscripts/2025/WGS_blue_whales/ngsadmix/evalAdmix/ngsadmix_result.K_1.15.qopt"), header=FALSE)
plotAdmix(q=q, pop=pop, ord=ord, rotatelab = 90, padj = 0.05,cex.main =  1.5,cex.lab = 1,cex.inds = 1)
plotCorRes(cor_mat = r, pop = pop, ord=ord, title = "K=1", max_z=0.1, min_z=-0.1,
           rotatelabpop=90, cex.lab=0.8, cex.lab.2 = 0.8,cex.legend = 0.8,cex.main = 0.8,)
```

## use script:

[run_ngsadmix.sh](src/ngsadmix/run_ngsadmix.sh)

```
dataset=wp1_final_bal

for k in {2..10}; do
    for i in {1..10}; do
        sbatch \
            --job-name=${dataset}.ngsadmix \
            --output=logs/ngsadmix/ngsadmix.${dataset}.out \
            --error=logs/ngsadmix/ngsadmix.${dataset}.err \
            src/ngsadmix/run_ngsadmix.sh ${dataset} ${k} ${i}
    done
done
```