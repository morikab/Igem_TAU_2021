import gzip
from Bio import SeqIO
from modules.ORF.calculating_cai import general_geomean
from modules.user_IO.input_functions import calculate_cai_weights_for_input
import numpy as np
import os
import json

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
    with gzip.open(fid, "rt") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            cds_dict[record.description] = record.seq
    cai_weights = calculate_cai_weights_for_input(cds_dict, estimated_expression={}, expression_csv_fid=None)
    cai_scores = general_geomean(sequence_lst=cds_dict.values(), weights=cai_weights)

    cai_final_dict.update(cai_weights)
    cai_final_dict['std'] = np.std(cai_scores)
    cai_final_dict['avg'] = np.mean(cai_scores)

    return cai_weights

def genomic_rna_to_rrna(fid):
    rrna_dict= {}
    with gzip.open(fid, "rt") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if '5S ribosomal RNA' in record.description:
                rrna_dict['5s'] = str(record.seq)
            elif '23S ribosomal RNA' in record.description:
                rrna_dict['23s'] = str(record.seq)
            elif '16S ribosomal RNA' in record.description:
                rrna_dict['16s'] = str(record.seq)
    return rrna_dict




def cds_rna_dirs_to_cub_rrna_csv(cds_dir, rna_dir ):
    cds_suffix = '_cds_from_genomic.fna.gz' #todo: fix the way i write and reference genes
    rna_suffix = '_rna_from_genomic.fna.gz'
    final_dict = {}

    output_csv = open('../data/refseq_genomes/refseq_analysed.csv', 'w')


    output_csv.write('assembly, 5s, 23s, 16s, TTT, TTC, TTA, TTG, TCT, TCC, TCA, TCG, TAT, TAC, TGT, TGC, TGG, CTT, CTC, CTA, CTG, CCT, CCC, CCA, CCG, CAT, CAC, CAA, CAG, CGT, CGC, CGA, CGG, ATT, ATC, ATA, ATG, ACT, ACC, ACA, ACG, AAT, AAC, AAA, AAG, AGT, AGC, AGA, AGG, GTT, GTC, GTA, GTG, GCT, GCC, GCA, GCG, GAT, GAC, GAA, GAG, GGT, GGC, GGA, GGG, TGA, TAA, TAG\n')

    cds_dir_org = [file[:-24] for file in os.listdir(cds_dir)]
    rna_dir_org = [file[:-24] for file in os.listdir(rna_dir)]

    org_list = intersection(cds_dir_org, rna_dir_org)


    for org in org_list:
        org_dict = {}
        cds_path = os.path.join(cds_dir, org+cds_suffix)
        rna_path = os.path.join(rna_dir, org+rna_suffix)

        cai_dict = genomic_cds_to_cub(cds_path)
        rrna_dict = genomic_rna_to_rrna(rna_path)

        org_dict.update(rrna_dict)
        org_dict.update(cai_dict)

        final_dict[org] = org_dict

        csv_list = [org] + [str(i) for i in org_dict.values()]
        output_csv.write(', '.join(csv_list)+'\n')
    with open("../data/refseq_genomes/analysed_refseq.json", 'w') as handle:
        json.dump(final_dict, handle)


cds_gz = "data/GCF_003434225.1_ASM343422v1_cds_from_genomic.fna.gz"
rna_gz = "data/GCF_003434225.1_ASM343422v1_rna_from_genomic.fna.gz"
genomic_rna_to_rrna(rna_gz)

cds_rna_dirs_to_cub_rrna_csv("../data/refseq_genomes/ncbi_genome_cds", 'ncbi_genome_rna')
