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

## Population Structure Correction (omega-matrix) file
