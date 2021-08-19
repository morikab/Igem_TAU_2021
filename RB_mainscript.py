from RB_functions import *


EC_rest = REbase_org('Escherichia coli')  ##REbase of E.coli
## in E.coli there's 93 pure sites,341 repeated sites,and 8 unkown "?" ==442 sites 
BC_rest = REbase_org('Bacillus subtilis') ##REbase of Bacillus
##Bacillus there's 23 pure sites,17 repeated sites,and 6 unkown "?" ==442 sites 

EC_seq_list=[*EC_rest.keys()]
BC_seq_list=[*BC_rest.keys()]

EC_options_ntdict={seq:seq_opps(seq) for seq in EC_seq_list}
BC_options_ntdict={seq:seq_opps(seq) for seq in BC_seq_list}

print(BC_options_ntdict)
EC_options_aadict={}


for site_list in BC_options_ntlist:
    for site in site_list:
        
        [site,opptions]=seq2aaS(site)
        BC_aalist.append(opptions)
        BC_ntlist.append(site)
#### In the DB itself!!! ###         
'''BC_options_aalist=[seq2aaS(seq) for seq in BC_ntlist]'''

"""
no_there=[]
for enzy_EC in range(len(EC_seq_list)):
    for option_enzy_EC in EC_options_ntlist[enzy_EC]:
        for enzy_BC in range(len(EC_seq_list)):
            BC_seq_opp=BC_options_ntlist[enzy_BC+1]
            if ~BC_seq_opp.count(option_enzy_EC):
                x=1
                #no_there[enzy_EC+1]=enzy_BC"""



cds_path = r'data\cds_ecoli_new'
cds_seq = Seq('')
for record in SeqIO.parse(self.cds_path, "fasta"):
        cds_seq += record.seq            


def insert_site_CDS(cds,aalist1,ntlist1):
    AA_cds_seq=translate(cds)
    n=0
    for aasite in aalist1:
        n+=1
        if aasite in AA_cds_seq:
            stratindx=AA_cds_seq.find(aasite)
            endindx=stratindx+len(aasite)
            cds[3*startindx:3*endindx+1]=ntlist1[n]
            
    return cds

def remove_site_Plasmid(cds):
    AA_cds_seq=translate(cds)
    n=0
    for aasite in list1:
        n+=1
        if aasite in AA_cds_seq:
            stratindx=AA_cds_seq.find(aasite)
            endindx=stratindx+len(aasite)
            if cds[3*startindx:3*endindx+1]==ntlist1[n]:
                cds[3*startindx:3*endindx+1]=''
            
    return cds

