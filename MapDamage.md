```
ml PDCOLD/23.12 R/4.1.1 mapdamage/2.2.3

for sample in $(cat data/bamlists/full_herr.sample_list.txt); do
    echo "${sample}"
done



ibam=data/bams/${sample}.subsampled_3X.bam
ref=Reference/GCF_900700415.2_Ch_v2.0.2_genomic.fna
odir=data/mapdamage/${sample}

mapDamage -i ${ibam} -r ${ref} -d ${odir} -n 50000 --merge-reference-sequences

```