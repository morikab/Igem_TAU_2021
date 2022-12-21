# -*- coding: utf-8 -*-
import pandas as pd
import json
from Bio import SeqIO
from os import listdir
from os.path import isfile, join
import time

### filtering the 5_hit blast results by 95% pident, and keeping the best pident score for every genome
nt_to_aa = {
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
    'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
    'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
    'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_',
    'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W',
}

###create json with all information

def refseq_to_blast_name(refseq_name):
    blast_name = refseq_name.split('.')[0]
    return blast_name

def blast_to_refseq_name(blast_name, full_refseq_list):
    refseq_name = [i for i in full_refseq_list if refseq_to_blast_name(i) == blast_name][0]
    return refseq_name


def entry_to_data_files(entry:str, tls_files:list, n_hits = '200'):
    ###find the data files for each master study
    fasta_files = [tls_dir + f for f in tls_files if 'fsa_nt' in f and entry in f]
    data_files = [fasta[:-7]+'_'+n_hits+'_hits.csv' for fasta in fasta_files]
    return data_files

def tls_metadata(tls_metadata:pd.DataFrame(), tls_files:list, n_hits):
    blast_results_dict = {}
    for entry_name, row in tls_metadata.iterrows():
        tls_dict = row.to_dict()
        del tls_dict['fasta']
        tls_dict['files'] = entry_to_data_files(entry_name, tls_files, n_hits)
        blast_results_dict[entry_name] = tls_dict
    return blast_results_dict

def blast_df_to_dict(blast_df):
    'dict structure is that the key is the refseq id and the value is the highest percent of identity'
    'between the specific resfeq id (genome) and the tls scores (meaning the highest identity score between all tls '
    '(that correspond to a specific 16s denoted by its refseq id'

    n_seq = len(set(blast_df['sseqid'].to_list()))
    blast_df['pident'] = blast_df[blast_df['pident'] > 95]['pident']
    blast_df.sort_values(by=['evalue'], inplace=True)
    blast_df['match_len'] = blast_df['qend']-blast_df['qstart']
    avg_match_len = blast_df['match_len'].mean()

    blast_df.drop(['length', 'match_len', 'mismatch', 'gapopen', 'qstart', 'qend',
                        'sstart', 'send', 'bitscore'], axis= 1, inplace=True)
    grouped_blast = blast_df.groupby("sseqid")
    qseq_list = []
    evalue_list = []
    for read_name, read_hits in grouped_blast:
        read_hits_filtered = read_hits.loc[read_hits['evalue'] == read_hits['evalue'].min()]
        qseq_list += read_hits_filtered['qseqid'].to_list()
        evalue_list += read_hits_filtered['evalue'].to_list()
    scores_dict = dict(sorted(zip(qseq_list, evalue_list)))
    return scores_dict, n_seq, avg_match_len


def blastn_run(scores_dict, genomes_df, full_refseq_list):

    cai_columns = list(nt_to_aa.keys())
    cai_df = genomes_df.drop([i for i in genomes_df.columns if i not in cai_columns], inplace=False, axis = 1)
    non_cai_df = genomes_df.drop(cai_columns, inplace=False, axis=1)
    match_data = {}
    for genome_name, scores_dict in scores_dict.items():
        refseq = blast_to_refseq_name(genome_name, full_refseq_list)
        cai_values = cai_df.loc[refseq, :].transpose().to_dict()
        other_values_dict = non_cai_df.loc[refseq, :].transpose().to_dict()
        match_data[refseq] = {
            'scores': scores_dict,
            'cai': cai_values
        }
        match_data[refseq].update(other_values_dict)
    return match_data



def check_all_blast_res(genomes_df, tls_metadata:dict, out_fid:str, id_th, n_hits):
    full_refseq_list = genomes_df.index.to_list()
    blast_col = ['sseqid','qseqid',  'pident','length',
                        'mismatch', 'gapopen', 'qstart', 'qend',
                        'sstart', 'send', 'evalue', 'bitscore']

    for entry, entry_dict in tls_metadata.items():
        if isfile(out_fid + entry+ '_' +id_th+ '_' +n_hits+'.json'):
            continue
        if len(entry_dict['files']) == 0:
            continue
        blast_df = pd.DataFrame(columns= blast_col)
        for blast in entry_dict['files']:
            print(blast)
            blast_df_new = pd.read_csv(blast)
            blast_df_new.columns = blast_col
            blast_df = pd.concat([blast_df, blast_df_new], ignore_index=True)
            print(entry)
        scores_dict, n_seq, avg_match_len = blast_df_to_dict(blast_df)
        match_data = blastn_run(scores_dict, genomes_df, full_refseq_list)
        entry_dict['n_seq'] = n_seq
        entry_dict['avg_match_len'] = avg_match_len
        entry_dict['match_data'] = match_data
        save_data(entry, entry_dict, out_fid, id_th, n_hits)


def save_data(entry, blast_results_dict, out_fid, id_th, n_hits):
    with open(out_fid + entry+ '_' +id_th+ '_' +n_hits+'.json', "w") as outfile:
        json.dump(blast_results_dict, outfile)

if __name__ == "__main__":
    print('Start')
    tic = time.time()
    id_th = '95'
    n_hits = '5'
    genomes_df = pd.read_csv('../../data/processed_genomes/filtered/cai_and_16s_for_genomes_filtered.csv', index_col=0)
    genomes_df.drop(columns=['5s', '23s'], inplace=True)
    tls_metadata_df = pd.read_csv('../../data/processed_tls/tls_assembly_metadata.csv', index_col=0)
    tls_dir = '../../data/genbank_tls/'
    tls_files = [f for f in listdir(tls_dir) if isfile(join(tls_dir, f))]
    out_fid = '../../data/tls_genome_match/'
    blast_results_dict= tls_metadata(tls_metadata_df, tls_files, n_hits)
    check_all_blast_res(genomes_df, blast_results_dict, out_fid, id_th, n_hits)
    print('time ', time.time() - tic)







