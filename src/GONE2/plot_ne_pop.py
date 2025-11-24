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
POPULATIONS = ['maseskar_2003', 'dynekilen_1874', 'masthugget_1740', 'kampinge_1300', 'knastorp_600']
GENERATION_TIME: float = 4  # Generation time in years
SILL_PERIODS = [(1556, 1589), (1650, 1680), (1747, 1809), (1877, 1906)]

# Assign colors to populations
COLOR_MAP = {
    'maseskar_2003': 'C0',
    'dynekilen_1874': 'C1',
    'masthugget_1740': 'C2',
    'kampinge_1300': 'C3',
    'knastorp_600': 'C4'
}

# Process each population and compute median Ne values across iterations
for population in POPULATIONS:
    # Initialize plot
    plt.figure(figsize=(10, 6))
    pop_df = {}
    # Extract the year of sampling from the population string
    year_of_sampling: int = int(population.split('_')[-1])

    # Read Ne data from all 50 iterations
    for iteration in range(1, 51):
        ne_file = f"data/GONE2/output/{population}.{iteration}_GONE2_Ne"
        ne_data = pd.read_csv(ne_file, sep="\t")
        # Remove first 10 generations to avoid artifacts
        ne_data = ne_data[ne_data['Generation'] > 10]
        pop_df[iteration] = ne_data
        years = year_of_sampling - (ne_data['Generation'] * GENERATION_TIME)
        ne_values = ne_data['Ne_diploids']
        # Plot Ne trajectory for this population
        plt.plot(years, ne_values, marker=None, linestyle='-', 
                 label=population, color=COLOR_MAP[population], alpha=0.3)
    
    # Combine all iterations and compute median Ne per generation
    plot_df = pd.concat(pop_df.values(), ignore_index=True)
    ne_data = plot_df.groupby('Generation', as_index=False).median()
    
    # Convert generations to historical years
    generations = ne_data['Generation']
    years = year_of_sampling - (generations * GENERATION_TIME)
    ne_values = ne_data['Ne_diploids']
    
    # Plot Ne trajectory for this population
    plt.plot(years, ne_values, marker=None, linestyle='-', 
             label=population, color=COLOR_MAP[population])

    # Highlight sill periods (historical fishing restrictions)
    for start, end in SILL_PERIODS:
        plt.axvspan(start, end, color='gray', alpha=0.3)
    
    # Configure plot appearance
    plt.xlabel(f'Year (Generation time = {GENERATION_TIME} years)')
    plt.xticks(np.arange(1300, 2025, 100))
    plt.ylabel('Median Ne from 50 Iterations')
    plt.title('Effective Population (Ne) Size over Years')
    
    # Create legend with sill periods and populations
    import matplotlib.patches as mpatches
    sill_patch = mpatches.Patch(color='gray', alpha=0.3, label='Sillperioder')
    if population == "maseskar_2003":
        population_patches = [mpatches.Patch(color='C0', label='Maseskar (2003)')]
    elif population == "dynekilen_1874":
        population_patches = [mpatches.Patch(color='C1', label='Dynekilen (1874)')]
    elif population == "masthugget_1740":
        population_patches = [mpatches.Patch(color='C2', label='Masthugget (1740)')]
    elif population == "kampinge_1300":
        population_patches = [mpatches.Patch(color='C3', label='Kampinge (1300)')]
    elif population == "knastorp_600":
        population_patches = [mpatches.Patch(color='C4', label='Knastorp (600)')]
    plt.legend(handles=[sill_patch] + population_patches, loc='lower left')
    
    # Save the plot
    output_file = f"plots/GONE2/{population}.ne.png"
    plt.savefig(output_file)
    plt.close()    