from RB_functions import *
from Bio.Seq import *
from Bio import SeqIO




cds_path = r'final sequenecs IDT.fasta'
cds_seq = Seq('')
for record in SeqIO.parse(cds_path, "fasta"):
        cds_seq += record.seq           

BS_ntaa_dict,BS_NT_list = REseq_org('Bacillus subtilis') ##REbase of Bacillus
##Bacillus there's 23 pure sites,17 repeated sites,and 6 unkown "?" ==442 sites

cds_add_bs, add_bs_start_end_idex_dict=insert_site_CDS(cds_seq,BS_ntaa_dict)
cds_remove_bs=remove_site_from_plasmid(str(cds_seq), BS_NT_list)
cds_add_remove_bs=remove_site_from_plasmid(cds_add_bs, BS_NT_list)


bsfile = open("Bacillus_RE.txt", "w")

bsfile.write(">add_Bacillus_RE\n" +cds_add_bs + "\n")
bsfile.write(">Remove_Bacillus_RE\n" +cds_remove_bs+ "\n")
bsfile.write(">Add_Remove_Bacillus_RE\n" +cds_add_remove_bs + "\n")

bsfile.close()

#EC_ntaa_dict,EC_NT_list = REseq_org('Escherichia coli')  ##REbase of E.coli
## in E.coli there's 93 pure sites,341 repeated sites,and 8 unkown "?" ==442 sites 

#cds_add_ecoli, add_ecoli_start_end_idex_dict=insert_site_CDS(cds_seq,EC_ntaa_dict)
#cds_remove_ecoli=remove_site_from_plasmid(cds_seq, EC_NT_list)
#cds_add_remove_ecoli=remove_site_from_plasmid(cds_add_ecoli, EC_NT_list)



"""for site_list in BC_options_ntlist:
    for site in site_list:
        
        [site,opptions]=seq2aaS(site)
        BC_aalist.append(opptions)
        BC_ntlist.append(site)
#### In the DB itself!!! ###         
'''BC_options_aalist=[seq2aaS(seq) for seq in BC_ntlist]'''


no_there=[]
for enzy_EC in range(len(EC_seq_list)):
    for option_enzy_EC in EC_options_ntlist[enzy_EC]:
        for enzy_BC in range(len(EC_seq_list)):
            BC_seq_opp=BC_options_ntlist[enzy_BC+1]
            if ~BC_seq_opp.count(option_enzy_EC):
                x=1
                #no_there[enzy_EC+1]=enzy_BC"""""





