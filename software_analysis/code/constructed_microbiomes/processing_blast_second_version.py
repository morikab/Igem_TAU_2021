###create csv only with needed info

# -*- coding: utf-8 -*-
import pandas as pd
import json
from Bio import SeqIO
from os import listdir
from os.path import isfile, join
import time
###create json with all information

def refseq_to_blast_name(refseq_name):
    blast_name = refseq_name.split('.')[0]
    return blast_name

def tls_file_to_entry(file):
    return file.split('.')[1]

def blast_to_refseq_name(blast_name, full_refseq_list):
    refseq_name = [i for i in full_refseq_list if refseq_to_blast_name(i) == blast_name][0]
    return refseq_name


def single_blast_run(blast_df, genomes_df):
    'dict structure is that the key is the refseq id and the value is the highest percent of identity'
    'between the specific resfeq id (genome) and the tls scores (meaning the highest identity score between all tls '
    '(that correspond to a specific 16s denoted by its refseq id'
    full_refseq_list = genomes_df.index.to_list()

    #filter
    tic1 = time.time()
    blast_df.sort_values(by=['pident'], inplace=True)
    blast_df['pident'] = blast_df[blast_df['pident'] > 85]['pident']
    blast_df.dropna(inplace=True)
    tic2 = time.time()
    print(1, tic2-tic1)
    blast_df['qseqid'] = [blast_to_refseq_name(q, full_refseq_list) for q in blast_df['qseqid'].to_list()]
    tic3 = time.time()
    blast_df.set_index('qseqid', inplace=True)
    print(2, tic3-tic2)
    df_counts = blast_df['qseqid'].value_counts()
    tic4 = time.time()
    print(3, tic4-tic3)
    blast_df.join(df_counts)
    tic5 = time.time()
    print(4, tic5-tic4)
    blast_df.drop_duplicates(subset=['qseqid'], inplace=True, keep='last')
    blast_df['match_len'] = blast_df['qend']-blast_df['qstart']
    print(blast_df)
    blast_df = blast_df.join(genomes_df)
    print(blast_df)
    return blast_df


def check_all_blast_res(genomes_df, tls_dir, out_fid:str):
    blast_col = ['sseqid','qseqid',  'pident','length',
                        'mismatch', 'gapopen', 'qstart', 'qend',
                        'sstart', 'send', 'evalue', 'bitscore'] #fmt 10
    tls_files = [f for f in listdir(tls_dir) if isfile(join(tls_dir, f))]
    entries  = [tls_file_to_entry(f) for f in tls_files if 'mstr' in f]
    blast_files = [tls_dir + f for f in tls_files if '_5_hits.csv' in f ]
    for entry in entries:
        blast_df = pd.DataFrame(columns= blast_col)
        for blast_csv in [f for f in blast_files if entry in f]:
            blast_df_new = pd.read_csv(blast_csv)
            blast_df_new.columns = blast_col
            blast_df = pd.concat([blast_df, blast_df_new], ignore_index=True)
        entry_df = single_blast_run(blast_df, genomes_df)
        # save_data(entry, entry_df, out_fid)

def save_data(entry:str, entry_df:pd.DataFrame, out_fid:str):
    entry_df.to_csv(out_fid + entry+ '_5hits_85idth_with_genomes.csv')

if __name__ == "__main__":
    print('Start')
    tic = time.time()
    genomes_df = pd.read_csv('../../data/processed_genomes/filtered/cai_and_16s_for_genomes_filtered.csv', index_col=0)
    tls_metadata_df = pd.read_csv('../../data/processed_tls/tls_assembly_metadata.csv', index_col=0)
    tls_dir = '../../data/genbank_tls/'
    out_fid = '../../data/tls_genome_match/'
    check_all_blast_res(genomes_df, tls_dir, out_fid)
    print(time.time() - tic)
