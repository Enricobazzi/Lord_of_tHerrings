# Purpose: Create metadata file from CSV files I generated for each dataset from 
# the NCBI SRA database. This is hard-coded for the datasets we downloaded.

pdata=/cfs/klemming/projects/supr/naiss2024-6-170/raw_data/published_data

# add header historic
# DON'T RUN THIS IF OUR HISTORIC DATA IS ALREADY IN THE TABLE 
# echo "samplename_index_lane readgroup_id readgroup_platform path_to_R1_fastq_file path_to_R2_fastq_file" \
#     > config/historical_samples_paths.txt

# add header modern
echo "samplename_index_lane readgroup_id readgroup_platform path_to_R1_fastq_file path_to_R2_fastq_file" \
    > config/modern_samples_paths.txt

# for the Atmore et al. 2022 historical individual data:
project=PRJEB52723
paste -d' ' \
    <(
        for sample in $(grep 'historical' data/published_data/Atmore_2022.csv | cut -d',' -f14 | sort -u); do
            n=$(grep -w $sample <(grep 'historical' data/published_data/Atmore_2022.csv | cut -d',' -f14) | wc -l)
            if [[ $sample != $previous_sample ]]; then
                for i in $(seq 1 $n); do
                    echo "${sample}_${i}_${project} ${sample}.${project}.${i} illumina"
                done
                previous_sample=$sample
            fi
        done
    ) \
    <(
        for sample in $(grep 'historical' data/published_data/Atmore_2022.csv | cut -d',' -f14 | sort -u); do
            if [[ $sample != $previous_sample ]]; then
                fqs=($(grep $sample data/published_data/atamore_et_al2022_2024.lst | grep -v "M-$sample" | grep R1_001))
                for fq in ${fqs[*]}; do echo "${pdata}/Atmore_et_al2022_2024/${fq}"; done
                previous_sample=$sample
            fi
        done
    ) \
    <(
        for sample in $(grep 'historical' data/published_data/Atmore_2022.csv | cut -d',' -f14 | sort -u); do
            if [[ $sample != $previous_sample ]]; then
                fqs=($(grep $sample data/published_data/atamore_et_al2022_2024.lst | grep -v "M-$sample" | grep R2_001))
                for fq in ${fqs[*]}; do echo "${pdata}/Atmore_et_al2022_2024/${fq}"; done
                previous_sample=$sample
            fi
        done
    ) >> config/historical_samples_paths.txt

# for the Atmore et al. 2022 modern individual data:
project=PRJEB52723
paste -d' ' \
    <(
        for sample in $(grep 'modern' data/published_data/Atmore_2022.csv | cut -d',' -f14 | sort -u); do
            n=$(grep -w $sample <(grep 'modern' data/published_data/Atmore_2022.csv | cut -d',' -f14) | wc -l)
            if [[ $sample != $previous_sample ]]; then
                for i in $(seq 1 $n); do
                    echo "${sample}_${i}_${project} ${sample}.${project}.${i} illumina"
                done
                previous_sample=$sample
            fi
        done
    ) \
    <(
        for sample in $(grep 'modern' data/published_data/Atmore_2022.csv | cut -d',' -f14 | sort -u); do
            if [[ $sample != $previous_sample ]]; then
                msample=$(echo $sample | sed 's/M/M-/')
                fqs=($(grep $msample data/published_data/atamore_et_al2022_2024.lst | grep R1_001))
                for fq in ${fqs[*]}; do echo "${pdata}/Atmore_et_al2022_2024/${fq}"; done
                previous_sample=$sample
            fi
        done
    ) \
    <(
        for sample in $(grep 'modern' data/published_data/Atmore_2022.csv | cut -d',' -f14 | sort -u); do
            if [[ $sample != $previous_sample ]]; then
                msample=$(echo $sample | sed 's/M/M-/')
                fqs=($(grep $msample data/published_data/atamore_et_al2022_2024.lst | grep R2_001))
                for fq in ${fqs[*]}; do echo "${pdata}/Atmore_et_al2022_2024/${fq}"; done
                previous_sample=$sample
            fi
        done
    ) >> config/modern_samples_paths.txt

# for the Atmore et al. 2024 historical individual data:
project=PRJEB77597
paste -d ' ' \
    <(
        for sample in $(grep 'historical' data/published_data/Atmore_2024.csv | cut -d',' -f14); do
            n=$(grep -w $sample <(grep 'historical' data/published_data/Atmore_2024.csv | cut -d',' -f14) | wc -l)
            if [[ $sample != $previous_sample ]]; then
                for i in $(seq 1 $n); do
                    echo "${sample}_${i}_${project} ${sample}.${project}.${i} illumina"
                done
                previous_sample=$sample
            fi
        done
    ) \
    <(
        for run in $(tail -n +2 data/published_data/Atmore_2024.csv | cut -d',' -f1); do
            echo "${pdata}/Atmore_et_al2022_2024/${run}_1.fastq.gz ${pdata}/Atmore_et_al2022_2024/${run}_2.fastq.gz"
        done
    ) >> config/historical_samples_paths.txt

# for Fuentes-Pardo et al. 2024 modern individual data
project=PRJNA930418
paste -d' ' \
    <(
        for sample in $(grep 'individual' data/published_data/FuentesPardo_2024.csv | cut -d',' -f14); do
            n=$(grep -w $sample <(grep 'individual' data/published_data/FuentesPardo_2024.csv | cut -d',' -f14) | wc -l)
            if [[ $sample != $previous_sample ]]; then
                for i in $(seq 1 $n); do
                    echo "${sample}_${i}_${project} ${sample}.${project}.${i} illumina"
                done
                previous_sample=$sample
            fi
        done
    ) \
    <(
        for run in $(grep 'individual' data/published_data/FuentesPardo_2024.csv | cut -d',' -f1); do
            echo "${pdata}/FuentesPardo_et_al2024/${run}_1.fastq.gz ${pdata}/FuentesPardo_et_al2024/${run}_2.fastq.gz"
        done
    ) >> config/modern_samples_paths.txt

# for Han et al. 2020 modern individual data
project=PRJNA642736
paste -d' ' \
    <(
        for sample in $(grep 'individual' data/published_data/Han_2020.csv | cut -d',' -f14); do
            n=$(grep -w $sample <(grep 'individual' data/published_data/Han_2020.csv | cut -d',' -f14) | wc -l)
            if [[ $sample != $previous_sample ]]; then
                for i in $(seq 1 $n); do
                    echo "${sample}_${i}_${project} ${sample}.${project}.${i} illumina"
                done
                previous_sample=$sample
            fi
        done
    ) \
    <(
        for run in $(grep 'individual' data/published_data/Han_2020.csv | cut -d',' -f1); do
            echo "${pdata}/Han_et_al2020/${run}_1.fastq.gz ${pdata}/Han_et_al2020/${run}_2.fastq.gz"
        done
    ) >> config/modern_samples_paths.txt

# for Kongsstovu et al. 2022 modern individual data 
project=PRJEB25669
paste -d' ' \
    <(
        for sample in $(grep 'individual' data/published_data/Kongsstovu_2022.csv | cut -d',' -f14); do
            n=$(grep -w $sample <(grep 'individual' data/published_data/Kongsstovu_2022.csv | cut -d',' -f14) | wc -l)
            if [[ $sample != $previous_sample ]]; then
                for i in $(seq 1 $n); do
                    echo "${sample}_${i}_${project} ${sample}.${project}.${i} illumina"
                done
                previous_sample=$sample
            fi
        done
    ) \
    <(
        for run in $(grep 'individual' data/published_data/Kongsstovu_2022.csv | cut -d',' -f1); do
            echo "${pdata}/Kongsstovu_et_al2022/${run}_1.fastq.gz ${pdata}/Kongsstovu_et_al2022/${run}_2.fastq.gz"
        done
    ) >> config/modern_samples_paths.txt

# for Lamichhaney et al. 2017 modern individual data
project=PRJNA338612
paste -d' ' \
    <(
        for sample in $(grep 'individual' data/published_data/Lamichhaney_2017.csv | cut -d',' -f14); do
            n=$(grep -w $sample <(grep 'individual' data/published_data/Lamichhaney_2017.csv | cut -d',' -f14) | wc -l)
            if [[ $sample != $previous_sample ]]; then
                for i in $(seq 1 $n); do
                    echo "${sample}_${i}_${project} ${sample}.${project}.${i} illumina"
                done
                previous_sample=$sample
            fi
        done
    ) \
    <(
        for run in $(grep 'individual' data/published_data/Lamichhaney_2017.csv | cut -d',' -f1); do
            echo "${pdata}/Lamichhaney_et_al2017/${run}_1.fastq.gz ${pdata}/Lamichhaney_et_al2017/${run}_2.fastq.gz"
        done
    ) >> config/modern_samples_paths.txt

# for Martinez-Barrio et al. 2016 modern individual data
project=SRP056617
paste -d' ' \
    <(
        for sample in $(grep 'individual' data/published_data/MartinezBarrio_2016.csv | cut -d',' -f14); do
            n=$(grep -w $sample <(grep 'individual' data/published_data/MartinezBarrio_2016.csv | cut -d',' -f14) | wc -l)
            if [[ $sample != $previous_sample ]]; then
                for i in $(seq 1 $n); do
                    echo "${sample}_${i}_${project} ${sample}.${project}.${i} illumina"
                done
                previous_sample=$sample
            fi
        done
    ) \
    <(
        for run in $(grep 'individual' data/published_data/MartinezBarrio_2016.csv | cut -d',' -f1); do
            echo "${pdata}/MartinezBarrio_et_al2016/${run}_1.fastq.gz ${pdata}/MartinezBarrio_et_al2016/${run}_2.fastq.gz"
        done
    ) >> config/modern_samples_paths.txt
