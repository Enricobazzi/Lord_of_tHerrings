# Purpose: Create metadata file from CSV files I generated for each dataset from 
# the NCBI SRA database.

# for the Han et al. 2020 data:
paste \
    <(paste <(tail -n +2 data/published_data/Han_2020_individuals.csv | cut -d',' -f14 | cut -d'_' -f1) \
        <(tail -n +2 data/published_data/Han_2020_individuals.csv | cut -d',' -f6 | rev | cut -d'_' -f1 | cut -c 1) \
        <(yes PRJNA642736 | head -113) | tr '\t' '_') \
    <(paste <(tail -n +2 data/published_data/Han_2020_individuals.csv | cut -d',' -f14 | cut -d'_' -f1) \
        <(yes PRJNA642736 | head -113) \
        <(tail -n +2 data/published_data/Han_2020_individuals.csv | cut -d',' -f6 | rev | cut -d'_' -f1 | cut -c 1) \
        | tr '\t' '.') \
    <(yes illumina | head -113) \
    <(paste <(yes '/cfs/klemming/scratch/e/ebazzica/GenErode/data/Han_2020/' | head -113) \
        <(tail -n +2 data/published_data/Han_2020_individuals.csv | cut -d',' -f1) \
        <(yes '_1.fastq' | head -113) | tr -d '\t' ) \
    <(paste <(yes '/cfs/klemming/scratch/e/ebazzica/GenErode/data/Han_2020/' | head -113) \
        <(tail -n +2 data/published_data/Han_2020_individuals.csv | cut -d',' -f1) \
        <(yes '_2.fastq' | head -113) | tr -d '\t') \
    | tr '\t' ' ' >> config/modern_samples_paths.txt

# for the Atmore et al. 2022 data:
paste \
    <(paste <(grep 'modern' data/published_data/Atmore_2022.csv | cut -d',' -f14 | cut -d'_' -f3 | tr -d '-') \
        <(yes "1" | head -52) \
        <(yes PRJEB52723 | head -52) | tr '\t' '_') \
    <(paste <(grep 'modern' data/published_data/Atmore_2022.csv | cut -d',' -f14 | cut -d'_' -f3 | tr -d '-') \
        <(yes PRJEB52723 | head -52) \
        <(yes "1" | head -52) \
        | tr '\t' '.') \
    <(yes illumina | head -52) \
    <(paste <(yes '/cfs/klemming/projects/supr/naiss2024-6-170/raw_data/published_data/Atmore_et_al2022_2024/' | head -52) \
        <(
        for sample in $(grep 'modern' data/published_data/Atmore_2022.csv | cut -d',' -f14 | cut -d'_' -f3); do
            grep $sample data/published_data/atamore_et_al2022_2024.lst | grep R1_001
        done 
        ) | tr -d '\t') \
    <(paste <(yes '/cfs/klemming/projects/supr/naiss2024-6-170/raw_data/published_data/Atmore_et_al2022_2024/' | head -52) \
        <(
        for sample in $(grep 'modern' data/published_data/Atmore_2022.csv | cut -d',' -f14 | cut -d'_' -f3); do
            grep $sample data/published_data/atamore_et_al2022_2024.lst | grep R2_001
        done 
        ) | tr -d '\t') \
    | tr '\t' ' ' >> config/modern_samples_paths.txt
