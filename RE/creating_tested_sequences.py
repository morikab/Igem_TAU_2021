from RE.functions import *
import time
from Bio import SeqIO
import os
org2 = 'Bacillus subtilis'

base_path = os.path.join(os.path.dirname(__file__), 'example_data')
cds_nt = SeqIO.read('..\example_data\\mCherry_original.fasta',
                    "fasta").seq
cds_aa = translate(cds_nt)

tic2 = time.time()

EC_enzyme_dict = REbase_org('Escherichia coli K-12', cds_aa)
print('number of relevant Escherichia coli K-12 REs: ', len(EC_enzyme_dict))

EC_full_dict = REbase_org('Escherichia coli K-12', cds_aa) #should be 'Escherichia coli'
print('number of relevant Escherichia coli REs: ', len(EC_full_dict))


tic3 = time.time()
print(2, tic3-tic2)


BS_enzyme_dict = REbase_org('Bacillus subtilis 168', cds_aa) ##REbase of Bacillus
print('number of relevant Bacillus subtilis 168 REs: ', len(BS_enzyme_dict))

BS_full_dict = REbase_org('Bacillus subtilis 168', cds_aa) # should be 'Bacillus subtilis'
print('number of relevant Bacillus subtilis REs: ', len(BS_full_dict))

tic4 = time.time()
print(3, tic4-tic3)

add_BS_seq = insert_site_CDS(cds_nt, BS_enzyme_dict)
add_EC_seq = insert_site_CDS(cds_nt, EC_enzyme_dict)
tic5 = time.time()
print(4, tic5-tic4)

remove_EC_seq = remove_site_from_plasmid(add_BS_seq, EC_full_dict)
remove_BS_SEQ = remove_site_from_plasmid(add_EC_seq, BS_full_dict)
tic6 = time.time()
print(5, tic6-tic5)

all_res_dict = {'cds_add_bs':add_BS_seq, 'cds_add_ec':add_EC_seq,
                'add_bs_remove_ec':remove_EC_seq, 'add_ec_remove_bs':remove_BS_SEQ,
                'original': cds_nt}

# testing differences between sequences:
for seq_name1, seq1 in all_res_dict.items():
    print('\n', seq_name1)
    if sites_in_cds(EC_enzyme_dict, seq1):
        print('EC enzymes:', len(sites_in_cds(EC_enzyme_dict, seq1)))
        print(sites_in_cds(EC_enzyme_dict, seq1), '\n')
    if sites_in_cds(BS_enzyme_dict, seq1):
        print('BS enzymes:', len(sites_in_cds(BS_enzyme_dict, seq1)))
        print(sites_in_cds(BS_enzyme_dict, seq1), '\n')

    for seq_name2, seq2 in all_res_dict.items():
        change_num = sum(1 for a, b in zip(seq1, seq2) if a != b)
        print(seq_name1, seq_name2, change_num)
        if translate(seq1) != translate(seq2):
            raise ValueError('problem with site editing!!')


