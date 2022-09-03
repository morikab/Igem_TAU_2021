# -*- coding: utf-8 -*-
import pandas as pd
import json
destination_dir = '../../data/refseq_genomes/'

def refseq_to_blast_name(refseq_name):
    blast_name = refseq_name.split('.')[0]
    return blast_name

def blast_to_refseq_name(blast_name, full_refseq_list):
    refseq_name = [i for i in full_refseq_list if refseq_to_blast_name(i) == blast_name][0]
    return refseq_name

def blast_df_to_dict(blast_df, full_refseq_list):
    'dict structure is that the key is the refseq id and the value is the highest percent of identity'
    'between the specific resfeq id (genome) and the tls scores (meaning the highest identity score between all tls '
    '(that correspond to a specific 16s denoted by its refseq id'
    blast_name_list = blast_df.iloc[:,1].to_list()
    alignment_score_list = blast_df.iloc[:,2].to_list()
    blast_dict = {}
    for blast_name in sorted(set(blast_name_list)):
        refseq_name = blast_to_refseq_name(blast_name, full_refseq_list)
        blast_scores = [alignment_score_list[i] for i, i_blast_name in enumerate(blast_name_list)
                        if i_blast_name == blast_name]
        blast_dict[refseq_name] = max(blast_scores)
    return blast_dict


def blastn_run(blast_csv, genomes_df):
    blast_df = pd.read_csv(blast_csv)
    full_refseq_list = genomes_df.index.to_list()
    blast_dict = blast_df_to_dict(blast_df, full_refseq_list)

    #adding the alignment scores to the genome df
    tls_scores_for_df = []
    for refseq_name in full_refseq_list:
        try:
            tls_scores_for_df.append(blast_dict[refseq_name])
        except:
            tls_scores_for_df.append(None)
    genomes_df[blast_csv] = tls_scores_for_df

    # extracting data for the dict
    genome_idexes = [i for i, refseq_name in enumerate(full_refseq_list) if refseq_name in list(blast_dict.keys())]
    cai_weights = genomes_df.iloc[genome_idexes, list(range(4,68))].transpose().to_dict()
    return blast_dict, cai_weights, genomes_df


def tls_sequencing_info(tls_fasta):
    with open(tls_fasta, 'r') as fp:
        sequencing_content = fp.readlines()
        n_seqs = round(len(sequencing_content)/2)
        amplicon_len = len(sequencing_content[1])

    return n_seqs, amplicon_len


def check_all_blast_res(genomes_df, tls_new_metadata_df):
    tls_new_metadata_df = tls_new_metadata_df.iloc[:, 1:]

    blast_results_dict = {}
    for idx, blast_csv in enumerate(tls_new_metadata_df['blast_csv'].to_list()):
        tls_dict = tls_new_metadata_df.iloc[idx, :].to_dict()

        n_seqs, amplicon_len = tls_sequencing_info(blast_csv[:-3] + 'fsa_nt')
        tls_dict['n_seqs'] = n_seqs
        tls_dict['amplicon_len'] = amplicon_len

        blast_dict, cai_weights, genomes_df = blastn_run(blast_csv, genomes_df)
        tls_dict['cai'] = cai_weights
        tls_dict['alignment_score'] = blast_dict

        entry_name = blast_csv.split('.')[-3]
        blast_results_dict[entry_name] = tls_dict
        print(idx, n_seqs, amplicon_len)
    return blast_results_dict, genomes_df



def save_data(blast_results_dict, genomes_df, out_fid):
    genomes_df.to_csv(out_fid + 'tls_genome_matches.csv')
    with open(out_fid + 'tls_genome_matches.json', "w") as outfile:
        json.dump(blast_results_dict, outfile)

if __name__ == "__main__":
    print('Start')

    genomes_df = pd.read_csv('../../data/processed_genomes/filtered/cai_and_16s_for_genomes_filtered.csv')
    genomes_df.set_index('Unnamed: 0', inplace=True, drop=True,)
    genomes_df.index.name = None

    tls_new_metadata_df = pd.read_csv('../../data/processed_tls/tls_assembly_metadata_with_blast.csv', index_col=None)
    out_fid = '../../data/tls_genome_match/'
    blast_results_dict, genomes_df = check_all_blast_res(genomes_df, tls_new_metadata_df)
    save_data(blast_results_dict, genomes_df, out_fid)







