"""
This script generates allele count (AC) files from allele frequency (AF) data from the published Pool-Seq dataset.

The allele count file is required for BayPass analysis.

Two columns are created for each population (one for reference allele counts and one for alternative allele counts).

NA values are replaced with 0 counts, and the data is converted to integers.
"""
import pandas as pd

# info table with sample sizes for each population
info_table = "data/selection_scans_poolseq/poolseqdata_info.csv"
info_table = pd.read_csv(info_table, sep=",")

# column names for the af files
colnames=["CHROM", "POS", "A_Kalix_Baltic_Spring", "B_Vaxholm_Baltic_Spring", "DalBoB_Atlantic_Autumn", "DalFB_Atlantic_Spring", "DalGeB_Atlantic_Autumn", "DalInB_Atlantic_Spring", "DalNsF_Atlantic_Autumn", "DalNsS_Atlantic_Spring", "G_Gamleby_Baltic_Spring", "HGS10_Downs_EnglishChannel_Winter", "HGS11_RingkobingFjord_NorthSea_Spring", "HGS12_BornholmBasin_Baltic_Autumn", "HGS15_NSSH_Atlantic_Spring", "HGS16_Orkney_NorthSea_Autumn", "HGS17_IsleOfMan_IrishSea_Autumn", "HGS18_CelticSea_Atlantic_AutumnWinter", "HGS19_TeelinBay_Atlantic_Winter", "HGS1_Riga_Baltic_Spring", "HGS20_CapeWrath_Atlantic_Spring", "HGS21_Hebrides_Atlantic_Mixed", "HGS22_CapeWrath_Atlantic_Autumn", "HGS23_Clyde_Atlantic_Spring", "HGS24_Landvik_Atlantic_Spring", "HGS25_Lindas_Atlantic_Spring", "HGS26_Lusterfjorden_Atlantic_Spring", "HGS27_Gloppen_Atlantic_Spring", "HGS2_Riga_Baltic_Spring", "HGS3_Riga_Baltic_Autumn", "HGS4_Riga_Baltic_Autumn", "HGS5_Schlei_Baltic_Autumn", "HGS6_Schlei_Baltic_Spring", "HGS71_Rugen_Baltic_Spring", "HGS72_Rugen_Baltic_Spring", "HGS8_KattegatNorth_Atlantic_Spring", "HGS9_Greenland_Atlantic_Spring", "HWS1_Japan_SeaOfJapan", "HWS2_PechoraSea_BarentsSea", "HWS3_WhiteSea_WhiteSea", "HWS4_KandalakshaBay_WhiteSea", "HWS5_KandalakshaBay_WhiteSea", "HWS6_Balsfjord_Atlantic", "H_Fehmarn_Baltic_Autumn", "J_Traslovslage_Baltic_Spring", "LandvikS17_Norway_Baltic_Spring", "N_NorthSea_Atlantic_Autumn", "O_Hamburgsund_Atlantic_Spring", "PB10_Skagerrak_Atlantic_Spring", "PB11_Kalmar_Baltic_Spring", "PB12_Karlskrona_Baltic_Spring", "PB1_HastKar_Baltic_Spring", "PB2_Iceland_Atlantic_Spring", "PB4_Hudiksvall_Baltic_Spring", "PB5_Galve_Baltic_Spring", "PB6_Galve_Baltic_Summer", "PB7_Galve_Baltic_Autumn", "PB8_Pacific_Pacific_Spring", "PB9_Kattegat_Atlantic_Spring", "PN3_CentralBaltic_Baltic_Spring", "Q_Norway_Atlantic_Atlantic_Spring", "TysklandS18_Germany_Baltic"]
# list of populations to exclude (Pacific populations)
pacific_pops = ["HWS1_Japan_SeaOfJapan", "HWS2_PechoraSea_BarentsSea", "HWS3_WhiteSea_WhiteSea", "HWS4_KandalakshaBay_WhiteSea", "HWS5_KandalakshaBay_WhiteSea", "HWS6_Balsfjord_Atlantic", "PB8_Pacific_Pacific_Spring"]

# read the 500 af files:
for n in range(1, 501):
    nn = str(n).zfill(3)
    af_file = f"data/selection_scans_poolseq/subsets_freqs/60.Neff.{nn}.freq"
    af = pd.read_csv(af_file, sep="\t", header=None)
    af.columns = colnames
    af = af.drop(columns=pacific_pops)
    # create ac dataframe
    new_colnames = [f'{name}_{suffix}' for name in af.columns[2:] for suffix in ("ref", "alt" )]
    nrows = af.shape[0]
    ac = pd.DataFrame(index=range(nrows), columns=new_colnames)
    # fill ac dataframe with af values multiplied by sample size from info_table
    for col in af.columns[2:]:
        # Find the row where 'Sample name' matches the column name
        sample_size = info_table.loc[info_table['Sample name'] == col, 'Sample size']
        if not sample_size.empty:
            sample_size_value = sum(sample_size.values)
        else:
            print(f'No sample size found for {col}')
        ac[f'{col}_alt'] = round(af[col] * sample_size_value)
        ac[f'{col}_ref'] = sample_size_value - ac[f'{col}_alt']
    # set missing data to 0
    nona_snps = ac.dropna().shape[0]
    ac = ac.fillna(0)
    # convert to int
    ac = ac.astype(int)
    # save ac file - without index and header (only data)
    out_file = f"data/selection_scans_poolseq/subsets_acs/60.Neff.{nn}.ac"
    ac.to_csv(out_file, sep="\t", index=False, header=False)
    print(f'Saved {out_file} with shape {ac.shape} and {nona_snps} SNPs with no NA values')
