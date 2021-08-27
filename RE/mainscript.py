from RE.functions import REbase_org, translate
import time

org1 = 'Escherichia coli K-12'
org2 = 'Bacillus subtilis 168'


tic1 = time.time()
cds_nt = user_inp2.seq
tic2 = time.time()
cds_aa = translate(cds_nt)




EC_ntaa_dict = REbase_org('Escherichia coli', cds_aa)  ##REbase of E.coli
tic3 = time.time()
print(2, tic3-tic2)
BS_ntaa_dict = REbase_org('Bacillus subtilis', cds_aa) ##REbase of Bacillus
tic4 = time.time()
print(3, tic4-tic3)

#create final sequence variable

final_seq = 'str after ORF optimization, taking out deoptimized and adding optimized'







# EC_ntaa_dict,EC_NT_list = REseq_org('Bacillus stearothermophilus Z130', cds_seq)  ##REbase of E.coli
# print(EC_ntaa_dict,EC_NT_list)
#
# cds_add_bs, add_bs_start_end_idex_dict=insert_site_CDS(cds_seq,BS_ntaa_dict)
# cds_add_ec, add_ecoli_start_end_idex_dict=insert_site_CDS(cds_seq,EC_ntaa_dict)
#
# add_bs_remove_ec=remove_site_from_plasmid(cds_add_bs, BS_NT_list)
# add_ec_remove_bs=remove_site_from_plasmid(cds_add_ec, EC_NT_list)
#
# all_res_dict = {'cds_add_bs':cds_add_bs, 'cds_add_ec':cds_add_ec,
#                 'add_bs_remove_ec':add_bs_remove_ec, 'add_ec_remove_bs':add_ec_remove_bs}
#
# for seq_name1, seq1 in all_res_dict.items():
#     for seq_name2, seq2 in all_res_dict.items():
#         change_num = sum(1 for a, b in zip(seq1, seq2) if a != b)
#         print(seq_name1, seq_name2, change_num)
#
#
# # #TODO: test, turn into multi organisis
# # print(all_res_dict)
# # write_fasta('all_res_dict', all_res_dict.values(), all_res_dict.keys())