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
        '--idir',
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

def get_ngi_ids(si_df):
    """
    Get the samples from the sample_info file.
    """
    # get the samples
    ngi_ids = [str(s) for s in si_df['NGI ID']]
    return ngi_ids

def get_sample_and_idx_from_ngi_id(ngi_id, si_df):
    """
    Get the sample from the sample_info file based on the NGI ID.
    """
    # get the sample
    sample_idx = si_df[si_df['NGI ID'] == ngi_id]["User ID"]
    samples = [str(s).split('_')[0] for s in sample_idx]
    idxs = [int(s.split('_')[1]) for s in sample_idx]
    if len(samples) == 0:
        raise ValueError(f'No sample found for NGI ID {ngi_id}')
    if len(samples) > 1:
        raise ValueError(f'More than one sample found for NGI ID {ngi_id}')
    sample = samples[0]
    if len(idxs) == 0:
        raise ValueError(f'No index found for NGI ID {ngi_id}')
    if len(idxs) > 1:
        raise ValueError(f'More than one index found for NGI ID {ngi_id}')
    idx = idxs[0]
    return sample, idx

def get_rgids_lanes_fastqs_pairs_from_lst(lst_file, seq_dir):
    """
    Get the fastq files from the .lst file.
    """
    with open(lst_file, 'r') as f:
        lines = f.readlines()
    fastqs = [line.strip() for line in lines if line.strip().endswith('.fastq.gz')]
    if len(fastqs) == 0:
        raise ValueError(f'No fastq files found in {lst_file}')
    if len(fastqs) % 2 != 0:
        raise ValueError(f'Uneven number of fastq files found in {lst_file}')
    # create pairs based on same name where the ending is:
    # R1_001.fastq.gz in one and R2_001.fastq.gz in the other
    fastq_bases = set([f.replace('_R1_001.fastq.gz', '').replace('_R2_001.fastq.gz', '') for f in fastqs])
    fastq_pairs = [f'{seq_dir}/{fb}_R1_001.fastq.gz {seq_dir}/{fb}_R2_001.fastq.gz' for fb in fastq_bases]
    # lane is assumed to be the third element in the fastq base name
    lanes = [fb.split('/')[-1].split('_')[3] for fb in fastq_bases]
    # Rigid ID is assumed to be the last element in the name of the folder with the fastq files
    rgids = [fb.split('/')[-2].split('_')[-1] for fb in fastq_bases]
    if len(rgids) != len(lanes):
        raise ValueError('Number of RGIDs and lanes do not match')
    if len(lanes) != len(fastq_pairs):
        raise ValueError('Number of lanes and fastq pairs do not match')
    if len(set(lanes)) != len(lanes):
        raise ValueError('Lanes are not unique')
    return rgids, lanes, fastq_pairs

def main():
    args = parse_args()
    seq_dir = args.idir
    if not os.path.exists(seq_dir):
        raise ValueError(f'Sequence directory {seq_dir} does not exist')
    si_df = get_si_df(seq_dir)
    ngi_ids = get_ngi_ids(si_df)
    for ngi_id in ngi_ids:
        sample, idx = get_sample_and_idx_from_ngi_id(ngi_id, si_df)
        lst_file = f'{seq_dir}/{ngi_id}.lst'
        if not os.path.exists(lst_file):
            raise ValueError(f'No .lst file found for NGI ID {ngi_id}')
        rgids, lanes, fastqs = get_rgids_lanes_fastqs_pairs_from_lst(lst_file, seq_dir)
        for i in range(len(lanes)):
            print(f'{sample}_{idx}_{lanes[i]} {rgids[i]}.{lanes[i]}.{idx} illumina {fastqs[i]}')

if __name__ == '__main__':
    main()

