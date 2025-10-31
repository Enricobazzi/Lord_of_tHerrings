## install

https://academic.oup.com/mbe/article/39/6/msac119/6596627

```
git clone https://github.com/lz398/distAngsd.git

module load eigen/3.4.0
module load htslib/1.20
module load gsl/2.7

# modify all files in folder so that:
# #include <eigen3/Eigen/Core> -> #include <Eigen/Core>
# include <eigen3/Eigen/Eigenvalues> -> #include <Eigen/Eigenvalues>
# HTS_LIBDIR=$(realpath $(HTSSRC))/libhts.a -> HTS_LIBDIR=$(realpath $(HTSSRC))/../lib/libhts.a
# LDFLAGS=-lz -llzma -lbz2 -lpthread -lcurl -lgsl -> LDFLAGS=-lcrypto -ldeflate -lz -llzma -lbz2 -lpthread -lcurl -lgsl -lgslcblas

make \
    EIGEN=/pdc/software/eb/software/eigen/3.4.0/include \
    HTSSRC=/pdc/software/eb/software/htslib/1.20/include \
    distAngsd
```

## shortcut to run the program

```
module load eigen/3.4.0
module load htslib/1.20
module load gsl/2.7
distAngsd=~/scripts/distAngsd/distAngsd
${distAngsd}
```

## input files

from the [distAngsd README](https://github.com/lz398/distAngsd?tab=readme-ov-file#preparation-and-command-examples)
```
ml bcftools

ref=/cfs/klemming/projects/supr/naiss2024-6-170/analyses/Generode_Herring/Reference/GCA_040183275.1_Ch_v3.0_genomic.fna
BAMS=/cfs/klemming/projects/supr/naiss2024-6-170/analyses/Bams_Herring
sample1=Z12
sample2=MHER008

if [[ -f ${BAMS}/${sample1}.merged.rmdup.merged.realn.rescaled.bam ]]; then
  bam1=${BAMS}/${sample1}.merged.rmdup.merged.realn.rescaled.bam
elif [[ -f ${BAMS}/${sample1}.merged.rmdup.merged.realn.bam ]]; then
  bam1=${BAMS}/${sample1}.merged.rmdup.merged.realn.bam
else
  echo "BAM file for sample 1 ${sample1} not found!"
  exit 1
fi

if [[ -f ${BAMS}/${sample2}.merged.rmdup.merged.realn.rescaled.bam ]]; then
  bam2=${BAMS}/${sample2}.merged.rmdup.merged.realn.rescaled.bam
elif [[ -f ${BAMS}/${sample2}.merged.rmdup.merged.realn.bam ]]; then
  bam2=${BAMS}/${sample2}.merged.rmdup.merged.realn.bam
else
  echo "BAM file for sample 2 ${sample2} not found!"
  exit 1
fi

outbcf=pair_bcfs/${sample1}_${sample2}.mpileup_noindel.bcf

# for distAngsd-geno:
bcftools mpileup -Ou -f ${ref} \
    ${bam1} ${bam2} | \
    bcftools filter -Ou -e INFO/INDEL!=0 \
    > ${outbcf}
# bgzip ${outbcf}
# bcftools index ${outbcf}.gz
bcftools index ${outbcf}
```

## run
```
module load eigen/3.4.0
module load htslib/1.20
module load gsl/2.7
distAngsd=~/scripts/distAngsd/distAngsd

sample1=Z12
sample2=MHER008
bcf=pair_bcfs/${sample1}_${sample2}.mpileup_noindel.bcf
ofile=ofiles/${sample1}_${sample2}

${distAngsd} -o ${ofile} -vcf ${bcf} 
```