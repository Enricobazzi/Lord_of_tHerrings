# Keep track of information regarding the samples

These are notes on where stuff is and where information is found.

## Published data

In order to download and unpack some of the fastq files I've created a conda environment where I installed the [SRA tools](https://github.com/ncbi/sra-tools) software:
```
conda create --name sra-tools
conda activate sra-tools
conda install -c bioconda sra-tools
```

Then I can use `fastq-dump` to extract the reads from a SRR file like this:
```
prefetch -X 104857600 <SRR_code> # 104857600 is 100GB for maximum file size
fastq-dump --split-3 <SRR_code> --outdir <output_directory>
```

The published herring data from [Han et al. 2020](https://elifesciences.org/articles/61076), [Atmore et al. 2022](https://www.pnas.org/doi/10.1073/pnas.2208703119), [Atmore et al. 2024](https://onlinelibrary.wiley.com/doi/10.1111/gcb.70010), [Lamichhaney et al. 2012](https://www.pnas.org/doi/full/10.1073/pnas.1216128109), [Lamichhaney et al. 2017](https://doi.org/10.1073/pnas.1617728114),
was downloaded by me and Nic from either the SRA ftp server or using sra-tools. The metadata for their respective projects can be found at:
- Han et al. 2020: [PRJNA642736](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA642736&o=acc_s%3Aa)
- Atmore et al. 2022 [PRJEB52723](https://www.ncbi.nlm.nih.gov/Traces/study/?page=2&acc=PRJEB52723&o=acc_s%3Aa)
- Atmore et al. 2024 [PRJEB77597](https://www.ncbi.nlm.nih.gov/Traces/study/?page=2&acc=PRJEB77597&o=acc_s%3Aa)
- Lamichhaney et al. 2017 [PRJNA338612](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA338612&o=acc_s%3Aa)
- Lamichhaney et al. 2012 [SRA057909](https://www.ncbi.nlm.nih.gov/sra/?term=sra057909)
- Martinez et al. 2016 [SRP056617](https://trace.ncbi.nlm.nih.gov/Traces/study/?acc=SRP056617&o=acc_s%3Aa), [SRP017094](https://trace.ncbi.nlm.nih.gov/Traces/study/?acc=SRP017094&o=acc_s%3Aa), [SRP017095](https://trace.ncbi.nlm.nih.gov/Traces/study/?acc=SRP017095&o=acc_s%3Aa)
- Hill et al. 2019 [PRJEB32358](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJEB32358&o=acc_s%3Aa)
- Fuentes-Pardo et al. 2024 [PRJNA930418](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA930418&o=acc_s%3Aa)
- Kongsstovu et al. 2022 [PRJEB25669](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJEB25669&o=acc_s%3Aa)

As SRR accession codes for each sample in each published study have been saved in files called `SRR_Acc_List.txt` within each folder, I can download and unpack them:
```
cd </path/to/study/folder/>

conda activate sra-tools
for acc in $(cat SRR_Acc_List.txt); do 
    echo $acc
    prefetch -X 104857600 $acc
    fastq-dump --split-3 $acc
    gzip ${acc}_1.fastq
    gzip ${acc}_2.fastq
done
```

## New sequences

I save the sample information tables that come with the sequences in `data/sampleinfo_files`. These are used to write metadata files for the GenErode pipeline.

The ones I have now are:
```
P21365
```
