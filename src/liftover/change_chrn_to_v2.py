# this script changes chromosome names from chrn to v2 format in a BED file.

# load arguments:
import argparse
argparser = argparse.ArgumentParser(description="Change chromosome names from v3 to v2 in a BED file.")
argparser.add_argument("--map_file", type=str, help="Path to the mapping file")
argparser.add_argument("--input_file", type=str, help="Path to the input BED file")
argparser.add_argument("--output_file", type=str, help="Path to the output BED file")
args = argparser.parse_args()

map_file = args.map_file
input_file = args.input_file
output_file = args.output_file

# load the mapping into a dictionary
mapping = {}
with open(map_file) as mf:
    for line in mf:
        v3, v2, chrn = line.strip().split()
        mapping[chrn] = v2

# process the input file
with open(input_file) as infile, open(output_file, "w") as outfile:
    for line in infile:
        for chrn, v2 in mapping.items():
            # Replace whole words only (chrX not inside other words)
            line = line.replace(f"{chrn}\t", f"{v2}\t")
            # Add other variants as needed depending on input format
        outfile.write(line)
