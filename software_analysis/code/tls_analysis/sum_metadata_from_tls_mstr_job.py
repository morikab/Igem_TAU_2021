# -*- coding: utf-8 -*-

from Bio import SeqIO
import os
import json
import pandas as pd
#in progress! just started working on this


def intersection(a, b):
    c = set(a).intersection(b)
    return list(c)

def mstr_fid_to_microbiome_id(mstr_fid):
    return mstr_fid.split('.')[-3]

def filter_metadata_in_mstr(mstr_fid):
    'finds the relevent information about tls files from the mstr files and returns it in a form of a dictionary'

    mstr_file = SeqIO.read(mstr_fid, 'genbank')
    id =  mstr_fid_to_microbiome_id(mstr_fid)
    relevant_metadata = {}
    if '16S' in mstr_file.description or '16s' in mstr_file.description:
        try:
            relevant_metadata['date'] = mstr_file.annotations['date']
        except:
            relevant_metadata['date'] = None

        try:
            relevant_metadata['title'] = mstr_file.annotations['references'][0].title
        except:
            relevant_metadata['title'] = None

        try:
            relevant_metadata['Assembly Method'] = mstr_file.annotations['structured_comment']['Assembly-Data']['Assembly Method']
        except:
            relevant_metadata['Assembly Method'] = None

        try:
            relevant_metadata['Sequencing Technology'] = mstr_file.annotations['structured_comment']['Assembly-Data']['Sequencing Technology']
        except:
            relevant_metadata['Sequencing Technology'] = None

        try:
            relevant_metadata['isolation_source'] = '_'.join(mstr_file.features[0].qualifiers['isolation_source'])
        except:
            relevant_metadata['isolation_source'] = None

        try:
            relevant_metadata['mol_type'] = '_'.join(mstr_file.features[0].qualifiers['mol_type']) ###should be 'genomic DNA' for all
        except:
            relevant_metadata['mol_type'] = None

        try:
            relevant_metadata['host'] = '_'.join(mstr_file.features[0].qualifiers['host'])
        except:
            relevant_metadata['host'] = None

    return id, relevant_metadata


def fasta_mstr_files_to_use(base_fid):
    'creates a list of dictionaries, where each dictionary contains the metadata and fid of the fastafile'
    'if a mstr fiel has no accompanying fasta, it will be printed out to the user'
    fasta_files = [os.path.join(base_fid, file) for file in os.listdir(base_fid) if '.1.fsa_nt' in file]
    tls_assem = {}
    mstr_files = [os.path.join(base_fid, file_name)
                     for file_name in os.listdir(base_fid) if '.mstr.gbff' in file_name]
    print(mstr_files)
    for mstr_fid in mstr_files:
        print('\n', mstr_fid)
        id, metadata = filter_metadata_in_mstr(mstr_fid)
        print(id)
        if metadata:
            fasta_fid = [file for file in fasta_files if id in file]
            if fasta_fid:
                metadata['fasta'] = fasta_fid
                tls_assem[id] = metadata
            else:
                print(id, ' no fasta file was found')
        else:
            print(id, ' uses a non-16S tls')
    return tls_assem

def save_data(tls_assem, output_dir):
    'saving the files as json and csv files'

    with open(os.path.join(output_dir, 'tls_assembly_metadata.json'), 'w') as handle:
        json.dump(tls_assem, handle)

    assem_csv = pd.DataFrame(tls_assem).transpose()
    assem_csv.to_csv(os.path.join(output_dir, 'tls_assembly_metadata.csv'))






if __name__ == "__main__":
    print("Start")
    output_dir = '../../data/processed_tls'
    tls_assem = fasta_mstr_files_to_use('../../data/genbank_tls')
    save_data(tls_assem, output_dir)
    print("The end")

# def parse_tls(fid):
#     tls_seqs = {}
#
#     with open(fid, "rt") as handle:
#         for record in SeqIO.parse(handle, "fasta"):
#             tls_seqs[record.description] = record.seq
#             # print(record.description)
#             # print(record.seq)



# filter_metadata_in_mstr('../data/genbank_tls/id.vdb_wgsnc.0301.2019.KAAB.mstr.gbff')
# parse_tls('../data/genbank_tls/id.vdb_wgsnc.0301.2019.KAAB.1.fsa_nt')