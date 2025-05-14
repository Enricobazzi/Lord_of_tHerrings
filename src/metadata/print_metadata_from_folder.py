import os
import pandas as pd # type: ignore
import argparse

def parse_args():
    """
    Parse the command line arguments.
    """
    parser = argparse.ArgumentParser(
        description='Print metadata file (no header) from folder with sequencing data'
    )
    parser.add_argument(
        'idir',
        type=str,
        help='Path to the sequencing directory - Needs to be absolute if you want absolute paths to R1 and R2 files'
    )
    return parser.parse_args()

def get_lst_files(seq_dir):
    """
    Get the list of .lst files in the given directory.
    """
    lst_files = [
        os.path.join(seq_dir, f) for f in os.listdir(seq_dir) if f.endswith('.lst')
    ]
    return lst_files

def get_si_df(seq_dir):
    """
    Get the sample_info file from the given directory.
    """
    # get the sample_info file from the subdirectory 00-Reports (the file ending with *_sample_info.txt)
    sample_info_files = [
        os.path.join(seq_dir, '00-Reports', f)
        for f in os.listdir(os.path.join(seq_dir, '00-Reports'))
        if f.endswith('_sample_info.txt')
    ]
    if len(sample_info_files) == 0:
        raise ValueError('No sample_info file found in the 00-Reports subdirectory')
    if len(sample_info_files) > 1:
        raise ValueError('More than one sample_info file found in the 00-Reports subdirectory')
    sample_info_file = sample_info_files[0]
    return pd.read_csv(sample_info_file, sep='\t')

def get_samples(si_df):
    """
    Get the samples from the sample_info file.
    """
    # get the samples
    samples = [str(s) for s in si_df['User ID']]
    return samples

def get_idx(sample):
    """
    Get the index from the sample name.
    """
    if len(sample.split('_')) > 1:
        idx = sample.split('_')[1]
    else:
        idx = 1
    return idx

def get_ngi_id(sample, si_df):
    """
    Get the NGI ID from the sample name.
    """
    ngi_id = si_df.loc[si_df['User ID'] == sample, 'NGI ID'].values[0]
    return ngi_id

def get_r1_r2(lst_files, ngi_id, seq_dir):
    """
    Get the R1 and R2 files from the lst file.
    """
    # get the lst file for the ngi_id
    lst = [
        lst for lst in lst_files if ngi_id in lst
    ]
    if len(lst) == 0:
        raise ValueError(f'No lst file found for the sample {ngi_id}')
    if len(lst) > 1:
        raise ValueError(f'More than one lst file found for the sample {ngi_id}')
    lst = lst[0]
    # read the lst file
    with open(lst, 'r') as f:
        lines = f.readlines()
    if len(lines) == 0:
        raise ValueError(f'No lines found in the lst file: {lst}')
    # get the R1 and R2 lines
    r1 = [l for l in lines if 'R1' in l]
    r2 = [l for l in lines if 'R2' in l]
    if len(r1) == 0:
        raise ValueError(f'No R1 line found in the lst file: {lst}')
    if len(r2) == 0:
        raise ValueError(f'No R2 line found in the lst file: {lst}')
    if len(r1) > 1:
        raise ValueError(f'More than one R1 line found in the lst file: {lst}')
    if len(r2) > 1:
        raise ValueError(f'More than one R2 line found in the lst file: {lst}')
    r1 = f'{seq_dir}/{r1[0]}'
    r2 = f'{seq_dir}/{r2[0]}'
    return r1.strip(), r2.strip()

def get_lane(r1, r2):
    """
    Get the lane number from the R1 and R2 files.
    """
    # get the lane number from the R1 and R2 files
    lane = r1.split('/')[-1].split('_')[3]
    if lane != r2.split('/')[-1].split('_')[3]:
        raise ValueError(f'Lane number is different in R1 and R2 files: {lane} and {r2.split("/")[-1].split("_")[3]}')
    return lane

def get_sample_line(sample, idx, lane, r1, r2):
    """
    Get the sample line from the sample name.
    """
    if '_' in sample:
        sample_name = sample.split('_')[0]
    else:
        sample_name = sample
    # get the sample line from the sample name
    sample_line = f'{sample_name}_{idx}_{lane} {sample_name}.{lane}.{idx} illumina {r1} {r2}'
    return sample_line

def main():
    """
    Main function to create the metadata file.
    """
    args = parse_args()
    seq_dir = args.idir
    if not os.path.exists(seq_dir):
        raise ValueError(f'Sequence directory {seq_dir} does not exist')

    # get the list of .lst files
    lst_files = get_lst_files(seq_dir)

    # get the sample_info file
    si_df = get_si_df(seq_dir)

    # get the samples
    samples = get_samples(si_df)

    # print the metadata file
    for sample in samples:
        idx = get_idx(sample)
        ngi_id = get_ngi_id(sample, si_df)
        r1, r2 = get_r1_r2(lst_files, ngi_id, seq_dir)
        lane = get_lane(r1, r2)
        sample_line = get_sample_line(sample, idx, lane, r1, r2)
        print(sample_line)
    
if __name__ == '__main__':
    main()
