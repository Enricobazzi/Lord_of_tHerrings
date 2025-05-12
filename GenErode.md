# Notes on learning and running GenErode

## Configuration files

To run the pipeline two main configuration files need to be modified in the GenErode folder:
- `slurm/config.yaml` manages SLURM stuff (e.g. computation project, memory/cpu specifications)
- `config/config.yaml` manages what steps of the pipeline will be run

Both of these were copied and modified from the [template slurm file](https://github.com/NBISweden/GenErode/blob/main/config/slurm/profile/config_plugin_dardel.yaml) and the [template config file](https://github.com/NBISweden/GenErode/blob/main/config/config.yaml) available from the GenErode github folder.

I have replicated these here locally to easily modify them, but they should be located in the GenErode folder where the analysis is run on the cluster. To copy them:
```
scp config/config.yaml ebazzica@dardel.pdc.kth.se:/cfs/klemming/scratch/e/ebazzica/GenErode/config/config.yaml
scp slurm/config.yaml ebazzica@dardel.pdc.kth.se:/cfs/klemming/scratch/e/ebazzica/GenErode/slurm/config.yaml
```

## Samples MetaData files

To run alignments on a group of samples I need to create metadata files for both historical and modern samples separately (since they follow distinct pipelines). These are located in `config/historical_samples_paths.txt` and `config/modern_samples_paths.txt`.

Description on how to fill the information in these can be found at:
https://github.com/NBISweden/GenErode/wiki/2.-Requirements-&-pipeline-configuration#1-prepare-metadata-files-of-samples

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

