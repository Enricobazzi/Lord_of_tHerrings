# Notes on learning and running GenErode

## Installation

To install the pipeline I clone the github folder:
```
git clone https://github.com/NBISweden/GenErode.git
```
and create a conda environment (called `generode`) using the environment.yml file in the folder:
```
conda env create -n generode -f environment.yml
```

## Configuration files

To run the pipeline two main configuration files need to be modified in the GenErode folder:
- `slurm/config.yaml` manages SLURM stuff (e.g. computation project, memory/cpu specifications)
- `config/config.yaml` manages what steps of the pipeline will be run

Both of these were copied and modified from the [template slurm file](https://github.com/NBISweden/GenErode/blob/main/config/slurm/profile/config_plugin_dardel.yaml) and the [template config file](https://github.com/NBISweden/GenErode/blob/main/config/config.yaml) available from the GenErode github folder.

I have replicated these here locally to easily modify them, but they should be located in the GenErode folder where the analysis is run on the cluster. To copy them:
```
scp data/generode/config.yaml ebazzica@dardel.pdc.kth.se:/cfs/klemming/scratch/e/ebazzica/Herring1/GenErode/config/config.yaml
scp data/generode/slurm_config.yaml ebazzica@dardel.pdc.kth.se:/cfs/klemming/scratch/e/ebazzica/Herring1/GenErode/slurm/config.yaml
scp data/generode/*_samples_paths.txt ebazzica@dardel.pdc.kth.se:/cfs/klemming/scratch/e/ebazzica/Herring1/GenErode/config/
```

## Samples MetaData files

To run alignments on a group of samples I need to create metadata files for both historical and modern samples separately (since they follow distinct pipelines). These are located in `config/historical_samples_paths.txt` and `config/modern_samples_paths.txt`.

Description on how to fill the information in these can be found at:
https://github.com/NBISweden/GenErode/wiki/2.-Requirements-&-pipeline-configuration#1-prepare-metadata-files-of-samples

*MY SCRIPT TO OBTAIN METADATA FOR NEWLY SEQUENCED*
```
# header first
echo "samplename_index_lane readgroup_id readgroup_platform path_to_R1_fastq_file path_to_R2_fastq_file" \
    > historical_samples_paths.txt

# you need pandas - my base env on dardel has it installed
conda activate
pmff=/cfs/klemming/home/e/ebazzica/scripts/print_metadata_from_folder.py

python ${pmff} --idir /cfs/klemming/projects/supr/naiss2024-6-170/raw_data/Herring_DeepSeq_2/files/P36109 \
    >> historical_samples_paths.txt

python ${pmff} --idir /cfs/klemming/projects/supr/naiss2024-6-170/raw_data/Herring_DeepSeq_3/files/P37012 \
    >> historical_samples_paths.txt

# I did:
# grep -vE "ND157|ND158|ND164|ND165|ND169|ND184|ND319|ND325|ND463" historical_samples_paths.txt > tmp && mv tmp historical_samples_paths.txt
# because those samples were mapped already when Herring_DeepSeq_2 arrived and didn't receive any
# additional sequencing in Herring_DeepSeq_3
# If starting from scratch it's not needed.
```

## Run alignments

To align raw reads to the reference genome and make standard filtering of BAM files I set:
```
bam_rmdup_realign_indels: True
```
in the `config/config.yaml` file.

Check the [Snakemake file](https://github.com/NBISweden/GenErode/blob/main/Snakefile) to see which steps will be included.

By default `bam_rmdup_realign_indels: True` will include the following steps of the pipeline automatically if not run already:
```
workflow/rules/0.1_reference_genome_preps.smk
workflow/rules/0.2_repeat_identification.smk
workflow/rules/1.1_fastq_processing.smk
workflow/rules/2_mapping.smk
workflow/rules/3.1_bam_rmdup_realign_indels.smk
```

Since I know that `0.1_reference_genome_preps` and `0.2_repeat_identification` have already been run, I will comment them out in the Snakemake file to avoid unwanted repetitions of these analyses (they shouldn't happen but could by error).

## Run the pipeline

Once all files are ready (regardless of what steps I want to run) I run the pipeline like this:
```
# create a screen for the run:
screen -S <run>

# in the screen prepare the environment:
module load PDC bioinfo-tools apptainer tmux
conda activate generode

# dry run + main run:
snakemake --profile slurm -n &> <YYMMDD>_dry.out
snakemake --profile slurm &> <YYMMDD>_main.out
```

## Alignment QC

### Depth

After obtaining bam files run the [depth_from_bam.sh](src/mapping/depth_from_bam.sh) to extract alignment depth from the bam:
```
for bam in $(ls /cfs/klemming/projects/supr/naiss2024-6-170/analyses/Herring/data/bams/*.bam); do
  sbatch src/mapping/depth_from_bam.sh ${bam}
done

# Use this dirty script to get each sample's mean:
grep "#" /cfs/klemming/projects/supr/naiss2024-6-170/analyses/Herring/data/bams/*.depth | \
    rev | cut -d'/' -f1 | rev | sed 's/.depth:# mean=/ /' | cut -d' ' -f1,2
```

After adding depth to the wg_depth column of [samples_table.csv](data/samples_table.csv), plot the mean depth of each sample in a list:
```
# for ND samples:
python src/mapping/plot_depth.py \
    --sfile <(grep "ND" data/samples_table.csv | cut -d',' -f1) \
    --ofile plots/mapping/mean_depth_ND_samples.png
```