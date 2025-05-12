# Keep track of information regarding the samples

These are notes on where stuff is and where information is found.

## Published data

The published herring data from [Han et al. 2020](https://elifesciences.org/articles/61076), [Atmore et al. 2022](https://www.pnas.org/doi/10.1073/pnas.2208703119) and [Atmore et al. 2024](https://onlinelibrary.wiley.com/doi/10.1111/gcb.70010) was downloaded by Nic from the SRA ftp server. The metadata for their respective projects can be found at:
- Han et al. 2020: [PRJNA642736](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA642736&o=acc_s%3Aa)
- Atmore et al. 2022 [PRJEB52723](https://www.ncbi.nlm.nih.gov/Traces/study/?page=2&acc=PRJEB52723&o=acc_s%3Aa)
- Atmore et al. 2024 [PRJEB77597](https://www.ncbi.nlm.nih.gov/Traces/study/?page=2&acc=PRJEB77597&o=acc_s%3Aa)

The Han et al. data come in a packed form that needs to be unpacked before processing. In order to unpack the fastq files I've created a conda environment where I installed the [SRA tools](https://github.com/ncbi/sra-tools) software:
```
conda create --name sra-tools
conda activate sra-tools
conda install -c bioconda sra-tools
```

Then I can use `fastq-dump` to extract the reads from a SRR file like this:
```
fastq-dump --split-3 <SRR_file> --outdir <output_directory>

for acc in $(cat individual_accessions.lst); do
    echo "dumping $acc"
    fastq-dump --split-3 /cfs/klemming/projects/supr/naiss2024-6-170/raw_data/published_data/Han_et_al2020/${acc}
done
```

I will momentarily locate the fastq files in my scratch GenErode folder (in the subfolder `data/Han_2020`) waiting for more space in the storage folder.

## New sequences

I save the sample information tables that come with the sequences in `data/sampleinfo_files`. These are used to write metadata files for the GenErode pipeline.

The ones I have now are:
```
P21365
```


