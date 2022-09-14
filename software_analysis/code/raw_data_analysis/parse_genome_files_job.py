import pandas as pd
from Bio import SeqIO
import numpy as np
import os
import json
import sys
sys.path.insert(1, '../../../modules/user_IO')
sys.path.insert(1, '../../../modules/ORF')
sys.path.insert(1, '../../../modules/')
from input_functions import calculate_cai_weights_for_input
from calculating_cai import general_geomean
import csv


'''
calculate CAI and extract the rrna sequences from the refseq downloaded genomes
'''


def mean(data):
  n = len(data)
  mean = sum(data) / n
  return mean

def intersection(a, b):
    c = set(a).intersection(b)
    return c

def genomic_cds_to_cub(fid):
    cai_final_dict = {}
    cds_dict= {}
    for record in SeqIO.parse(fid, "fasta"):
        cds_dict[record.description] = record.seq
    cai_weights = calculate_cai_weights_for_input(cds_dict, estimated_expression={}, expression_csv_fid=None)
    cai_scores = general_geomean(sequence_lst=cds_dict.values(), weights=cai_weights)

    cai_final_dict.update(cai_weights)
    cai_final_dict['std'] = np.std(cai_scores)
    cai_final_dict['avg'] = np.mean(cai_scores)
    return cai_weights

def genomic_rna_to_rrna(fid):
    rrna_dict= {}
    for record in SeqIO.parse(fid, "fasta"):
        if '5S ribosomal RNA' in record.description:
            rrna_dict['5s'] = str(record.seq)
        elif '23S ribosomal RNA' in record.description:
            rrna_dict['23s'] = str(record.seq)
        elif '16S ribosomal RNA' in record.description:
            rrna_dict['16s'] = str(record.seq)
    return rrna_dict


def data_for_every_org(cds_dir, rna_dir, out_dir):
    cds_suffix = '_cds_from_genomic.fna'
    rna_suffix = '_rna_from_genomic.fna'
    final_dict = {}

    cds_dir_org = [file[:-len(cds_suffix)] for file in os.listdir(cds_dir)]
    rna_dir_org = [file[:-len(rna_suffix)] for file in os.listdir(rna_dir)]
    org_list = intersection(cds_dir_org, rna_dir_org)
    print('the number of cds files is: ', len(cds_dir_org),
          '\nthe number of rna files is: ', len(rna_dir_org),
          '\nthe number of genomes containing both a cds and a rna file is: ', len(org_list))
    print('files not in cds but in rna: ', [i for i in rna_dir_org if i not in cds_dir_org],
          '\nfiles not in rna but in cds: ', [i for i in cds_dir_org if i not in rna_dir_org])

    f = open(out_dir + '/cai_and_16s_for_genomes.csv', 'w')
    writer = csv.writer(f)
    writer.writerow(['5s', '16s', '23s', 'TTT', 'TTC', 'TTA', 'TTG', 'TCT', 'TCC', 'TCA', 'TCG', 'TAT', 'TAC', 'TGT', 'TGC', 'TGG', 'CTT', 'CTC', 'CTA', 'CTG', 'CCT', 'CCC', 'CCA', 'CCG', 'CAT', 'CAC', 'CAA', 'CAG', 'CGT', 'CGC', 'CGA', 'CGG', 'ATT', 'ATC', 'ATA', 'ATG', 'ACT', 'ACC', 'ACA', 'ACG', 'AAT', 'AAC', 'AAA', 'AAG', 'AGT', 'AGC', 'AGA', 'AGG', 'GTT', 'GTC', 'GTA', 'GTG', 'GCT', 'GCC', 'GCA', 'GCG', 'GAT', 'GAC', 'GAA', 'GAG', 'GGT', 'GGC', 'GGA', 'GGG', 'TGA', 'TAA', 'TAG'])


    for org in org_list:
        print(org)
        org_dict = {}
        cds_path = os.path.join(cds_dir + org+cds_suffix)
        rna_path = os.path.join(rna_dir+ org+rna_suffix)

        cai_dict = genomic_cds_to_cub(cds_path)
        rrna_dict = genomic_rna_to_rrna(rna_path)

        org_dict.update(rrna_dict)
        org_dict.update(cai_dict)

        final_dict[org] = org_dict

        row_to_write = list(org_dict.values())
        row_to_write.insert(0, org)
        writer.writerow(row_to_write)

    return final_dict


def write_fasta(fid, list_seq, list_name):
    ofile = open(fid + '.fasta', "w+")
    for i in range(len(list_seq)):
        ofile.write(">" + str(list_name[i]) + "\n" + str(list_seq[i]) + "\n")
    ofile.close()

def save_data(final_dict, out_dir):
    print(final_dict)
    with open(out_dir + "cai_and_16s_for_genomes3.json", 'w') as handle:
        json.dump(final_dict, handle)

    csv_data = pd.DataFrame(final_dict).transpose()
    csv_data.to_csv(out_dir+ "cai_and_16s_for_genomes3.csv")

    fasta_dict = {key: value['16s'] for key, value in final_dict.items()}
    write_fasta(out_dir+ "cai_and_16s_for_genomes3.fasta",
                list(fasta_dict.values()),
                list(fasta_dict.keys()))




if __name__ == "__main__":
    print("Start")
    output_dir = '../../data/processed_genomes_new/'

    org_dict = data_for_every_org("../../data/refseq_genomes/ncbi_genome_cds/",
                                  "../../data/refseq_genomes/ncbi_genome_rna/",
                                  out_dir = output_dir)
    save_data(org_dict, output_dir)
    print("The end")


