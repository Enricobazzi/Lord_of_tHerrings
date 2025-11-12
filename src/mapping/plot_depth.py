import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import argparse

# parse arguments - sfile and ofile
parser = argparse.ArgumentParser(description='Plot mean depth per sample with region colors.')
parser.add_argument('--sfile', type=str, default='data/samples.txt', help='Path to samples file')
parser.add_argument('--ofile', type=str, default='plots/mapping/mean_depth_ND_samples.png', help='Output plot file path')
args = parser.parse_args()

# read the samples table (one sample name per line)
sfile = args.sfile
samples = pd.read_csv(sfile, header=None)[0].tolist()
# output file
ofile = args.ofile
# color dictionary
region_colors = {
    'BALTIC': "#DE292C",
    'NE-ATL': "#567BD7FF",
    'NW-ATL': "#41B4F2",
    'TRANS': '#A035AF',
    'E-PACIFIC': '#FFFF00',
    'W-PACIFIC': '#00FF00',
    'ARCTIC': "#69F0E0",
    'UNKNOWN': '#808080'
}

# samples table:
table = "data/samples_table.csv"
df = pd.read_csv(table)
df = df[df['sample_id'].isin(samples) & (df['wg.depth'] > 0)]

# prepare data for plotting
samples = tuple(df['sample_id'])
y_pos = np.arange(len(samples))
mean_d = df['wg.depth'].to_numpy()
colors = df['region'].map(region_colors)

# Cap values at 6
max_value = 6
mean_d_capped = np.minimum(mean_d, max_value + 0.01)

# Create horizontal bar plot with color based on region
fig, ax = plt.subplots(figsize=(8, 20))
ax.barh(y_pos, mean_d_capped, align='center', color=colors)
ax.set_yticks(y_pos, labels=samples)
ax.invert_yaxis()  # labels read top-to-bottom

# Add vertical grid lines in light grey
ax.xaxis.grid(True, linestyle='--', which='both', color='lightgrey')
ax.set_xlabel('Mean Depth')
ax.set_title('Mean Depth per Sample')

# Add legend of colors
legend_handles = [mpatches.Patch(color=color, label=region) for region, color in region_colors.items()]
ax.legend(handles=legend_handles, title='Regions', bbox_to_anchor=(1.05, 1), loc='upper left')

# Add text for capped bars (actual value)
for i, (val, cap_val) in enumerate(zip(mean_d, mean_d_capped)):
    if val > max_value:
        ax.text(max_value + 0.05, i, f'...{val:.2f}X', va='center', ha='left', fontsize=7, color='darkred')

ax.margins(y=0.01)
plt.savefig(ofile, bbox_inches='tight', dpi=300)