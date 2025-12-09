"""
Plot effective population size (Ne) over time for herring populations.

This script reads GONE2 output files containing Ne estimates across iterations,
computes median values per generation, and visualizes Ne trends over historical years.
Sillperioder are highlighted in the plot.
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Configuration
POPULATIONS = ['faroe_2016', 'outeridefjord_2010', 'idefjord_2010', 'northsea_2016', 'norwegiansea_2015', 'maseskar_2003', 'dynekilen_1874', 'masthugget_1740', 'foldfjorden_1875', 'stavanger_1880']
POPULATIONS = ['outeridefjord_2010', 'maseskar_2003', 'risor_2008', 'dynekilen_1874', 'masthugget_1740', 'foldfjorden_1875', 'stavanger_1880']
POPULATIONS = ['faroe_2016', 'norwegiansea_2015', 'northsea_2016', 'maseskar_2003', 'dynekilen_1874', 'masthugget_1740', 'foldfjorden_1875', 'stavanger_1880']
GENERATION_TIME: float = 4  # Generation time in years
SILL_PERIODS = [(1556, 1589), (1650, 1680), (1747, 1809), (1877, 1906)]
OUTPUT_FILE = f"plots/GONE2/{'.'.join(POPULATIONS)}.ne.png"
# Assign colors to populations
COLOR_MAP = {
    'maseskar_2003': 'C0',
    'dynekilen_1874': 'C1',
    'masthugget_1740': 'C2',
    'kampinge_1300': 'C3',
    'knastorp_600': 'C4',
    'foldfjorden_1875': 'C3',
    'stavanger_1880': 'C4',
    'risor_2008': 'C5',
    'northsea_2016': 'C6',
    'norwegiansea_2015': 'C7',
    'faroe_2016': 'C8',
    'outeridefjord_2010': 'C9',
    'idefjord_2010': 'C10'
}   

# Initialize plot
plt.figure(figsize=(10, 6))

# Process each population and compute median Ne values across iterations
for population in POPULATIONS:
    pop_df = {}
    
    # Read Ne data from all 50 iterations
    for iteration in range(1, 51):
        ne_file = f"data/GONE2/output/{population}.{iteration}_GONE2_Ne"
        ne_data = pd.read_csv(ne_file, sep="\t")
        # Remove first generations to avoid artifacts
        ne_data = ne_data[ne_data['Generation'] > 4]
        pop_df[iteration] = ne_data
    
    # Combine all iterations and compute median and CI Ne per generation
    plot_df = pd.concat(pop_df.values(), ignore_index=True)
    ne_data = plot_df.groupby('Generation', as_index=False).median()
    #CI_lower = plot_df.groupby('Generation')['Ne_diploids'].quantile(0.025).reset_index(name='CI_lower')
    #CI_upper = plot_df.groupby('Generation')['Ne_diploids'].quantile(0.975).reset_index(name='CI_upper')
    #ne_data = ne_data.merge(CI_lower, on='Generation').merge(CI_upper, on='Generation')
    
    # Convert generations to historical years
    year_of_sampling: int = int(population.split('_')[-1])
    generations = ne_data['Generation']
    years = year_of_sampling - (generations * GENERATION_TIME)
    ne_values = ne_data['Ne_diploids']
    
    # Plot Ne trajectory for this population
    plt.plot(years, ne_values, marker=None, linestyle='-', 
             label=population, color=COLOR_MAP[population])
    # Plot confidence intervals
    #plt.fill_between(years, ne_data['CI_lower'], ne_data['CI_upper'], 
    #                 color=COLOR_MAP[population], alpha=0.2)

# Highlight sill periods (historical fishing restrictions)
for start, end in SILL_PERIODS: 
    plt.axvspan(start, end, color='gray', alpha=0.3)

# Configure plot appearance
plt.xlabel(f'Year (Generation time = {GENERATION_TIME} years)')
plt.xticks(np.arange(1300, 2025, 100))
plt.ylabel('Median Ne from 50 Iterations')
plt.title('Effective Population Size (Ne)\n in the last 120 generations')

# Create legend with sill periods and populations
import matplotlib.patches as mpatches
sill_patch = mpatches.Patch(color='gray', alpha=0.3, label='Sillperioder')
population_patches = []
for population in POPULATIONS:
    if population == "maseskar_2003":
        population_patches.append(mpatches.Patch(color='C0', label='Maseskär (2003)'))
    elif population == "dynekilen_1874":
        population_patches.append(mpatches.Patch(color='C1', label='Dynekilen (1874)'))
    elif population == "masthugget_1740":
        population_patches.append(mpatches.Patch(color='C2', label='Masthugget (1740)'))
    elif population == "foldfjorden_1875":
        population_patches.append(mpatches.Patch(color='C3', label='Foldfjorden (1875)'))
    elif population == "stavanger_1880":
        population_patches.append(mpatches.Patch(color='C4', label='Stavanger (1880)'))
    elif population == "risor_2008":
        population_patches.append(mpatches.Patch(color='C5', label='Risor (2008)'))
    elif population == "northsea_2016":
        population_patches.append(mpatches.Patch(color='C6', label='North Sea (2016)'))
    elif population == "norwegiansea_2015":
        population_patches.append(mpatches.Patch(color='C7', label='Norwegian Sea (2015)'))
    elif population == "faroe_2016":
        population_patches.append(mpatches.Patch(color='C8', label='Faroe (2016)'))
    elif population == "outeridefjord_2010":
        population_patches.append(mpatches.Patch(color='C9', label='Outeridefjord (2010)'))
    elif population == "idefjord_2010":
        population_patches.append(mpatches.Patch(color='C10', label='Idefjord (2010)'))
plt.legend(handles=[sill_patch] + population_patches, loc='lower left')

# Save the plot
plt.savefig(OUTPUT_FILE)