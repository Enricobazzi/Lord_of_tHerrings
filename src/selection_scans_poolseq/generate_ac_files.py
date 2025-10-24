"""
This script generates allele count (AC) files from allele frequency (AF) data from the published Pool-Seq dataset.

The allele count file is required for BayPass analysis.

Two columns are created for each population (one for reference allele counts and one for alternative allele counts).

NA values are replaced with 0 counts, and the data is converted to integers.

Additionally, a file is generated for each subset that contains the CHROM and POS columns of the SNPs. 
"""
import pandas as pd

# info table with sample sizes for each population
info_table = "data/selection_scans_poolseq/poolseqdata_info.csv"
info_table = pd.read_csv(info_table, sep=",")

# column names for the af files
colnames=pd.read_csv("data/published_data/60.Neff.freq", nrows=0, sep='\t').columns.tolist()
# list of populations to exclude (Pacific populations)
pacific_pops = ["HWS1_Japan_SeaOfJapan", "HWS2_PechoraSea_BarentsSea", "HWS3_WhiteSea_WhiteSea", "HWS4_KandalakshaBay_WhiteSea", "HWS5_KandalakshaBay_WhiteSea", "HWS6_Balsfjord_Atlantic", "PB8_Pacific_Pacific_Spring"]

# read the 500 af files:
for n in range(1, 501):
    nn = str(n).zfill(3)
    af_file = f"data/selection_scans_poolseq/subsets_freqs/60.Neff.{nn}.freq"
    af = pd.read_csv(af_file, sep="\t", header=None)
    af.columns = colnames
    af = af.drop(columns=pacific_pops)
    # number of snps before filtering maf
    initial_snps = af.shape[0]
    # filter rows where total AF is not between 0.05 and 0.95
    af["Total_AF"] = af.iloc[:, 2:].mean(axis=1, skipna=True)
    af = af[(af["Total_AF"] >= 0.05) & (af["Total_AF"] <= 0.95)]
    af.reset_index(drop=True, inplace=True)
    # filter rows with more than 20% NA values
    af = af[af.iloc[:, 2:-1].isna().mean(axis=1) <= 0.2]  # exclude CHROM, POS and Total_AF columns
    af.reset_index(drop=True, inplace=True)
    # number of snps after filtering maf
    filtered_snps = af.shape[0]
    # create ac dataframe
    new_colnames = [f'{name}_{suffix}' for name in af.columns[2:-1] for suffix in ("ref", "alt" )]
    nrows = af.shape[0]
    ac = pd.DataFrame(index=range(nrows), columns=new_colnames)
    # fill ac dataframe with af values multiplied by sample size from info_table
    for col in af.columns[2:-1]:  # exclude CHROM, POS and Total_AF columns
        # Find the row where 'Sample name' matches the column name
        sample_size = info_table.loc[info_table['Sample name'] == col, 'Sample size']
        sample_size_value = sum(sample_size.values)
        ac[f'{col}_alt'] = round(af[col] * sample_size_value)
        ac[f'{col}_ref'] = sample_size_value - ac[f'{col}_alt']
    
    # OPTION 1 - set missing data to 0
    #nona_snps = ac.dropna().shape[0]
    ac = ac.fillna(0)
    
    # OPTION 2 - remove missing data rows
    # nona_snps = ac.dropna().shape[0]
    #ac = ac.dropna()
    
    # convert to int
    ac = ac.astype(int)
    # save ac file - without index and header (only data)
    out_file = f"data/selection_scans_poolseq/subsets_acs/60.Neff.{nn}.ac"
    ac.to_csv(out_file, sep="\t", index=False, header=False)
    print(f'Saved {out_file} with shape {ac.shape} from {initial_snps} SNPs before MAF and missing data filtering')
    # save snp positions file
    snp_pos = af[['CHROM', 'POS']].iloc[ac.index]
    snp_pos_file = f"data/selection_scans_poolseq/subsets_acs/60.Neff.{nn}.snps"
    snp_pos.to_csv(snp_pos_file, sep="\t", index=False, header=False)
    print(f'Saved {snp_pos_file} with shape {snp_pos.shape}')
