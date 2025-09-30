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
awk 'NR>1 {
  n = ((NR-2) % 500) + 1         # which sequence (1..500)
  nn = sprintf("%03d", n)        # zero-padded filename part
  out = "data/selection_scans_poolseq/subsets_freqs/60.Neff." nn ".freq"
  print >> out
  close(out)
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
