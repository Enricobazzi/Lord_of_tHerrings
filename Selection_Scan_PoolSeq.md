# Install Baypass
```
git clone git@forge.inrae.fr:mathieu.gautier/baypass_public.git
unzip baypass_public-v3.1.zip
cd baypass_public-v3.1/sources
make clean all FC=gfortran
```
# Input preparation

## Allele Count file

generate 500 SNP subsets
```
awk -v n=500 '{
    file = sprintf("data/selection_scans_poolseq/subsets_freqs/60.Neff.%03d.freq", ((NR-1) % n) + 1);
    print >> file
    close(file)
}' data/published_data/60.Neff.freq
```

convert allele frequency (AF) to allele count (AC):
```
python src/selection_scans_poolseq/generate_ac_files.py
```

## Environmental Data file

```
Rscript src/selection_scans_poolseq/get_env_data.R
```

## Population Structure Correction (omega-matrix) file

```
for n in {1..500}; do
    nn=$(printf "%03d\n" "$n")
    sbatch \
        --job-name=${nn}_core \
        --output=logs/core/${nn}.out \
        --error=logs/core/${nn}.err \
        baypass_core.sh ${nn}     
done
```

# Run!

## on indiviudal datasets

```
for n in {1..500}; do
    nn=$(printf "%03d\n" "$n")
    sbatch \
        --job-name=${nn}_aux \
        --output=logs/aux/${nn}.out \
        --error=logs/aux/${nn}.err \
        baypass_aux.sh ${nn}     
done
```

## join results

```

```
