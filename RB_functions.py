import re
from collections import defaultdict
from dnachisel import DnaOptimizationProblem, AvoidPattern,EnforceTranslation

f = open(r'REbase_16.8.txt')
content = f.readlines()
REbase = [x.strip() for x in content]

##open REbase from website : http://rebase.neb.com/rebase/rebase.f12.html
## REbase varible contain every line in the REbase.txt
#"""            Every enzyme data called Block

#every block of enzyme is :
#1. Enzyme Name
#2. Recognition Site
#3. Microrganism                       """


### Standard abbreviations to represent ambiguity Dic. from the descripition of REbase DB.
sa={'R':['G','A'],'Y':['C','T'],'M':['A','C'],'K':['G','T'],
    'S':['G','C'],'W':['A','T'],'B':['G','C','T'],
    'D':['A','G','T'],'H':['A','C','T'],'V':['A','G','C'],
    'N':['A','G','C','T']}


table = {
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

def translate(seq,table=table): #DONE
    protein = ""
    if len(seq) % 3 == 0:
        for i in range(0, len(seq), 3):
            codon = seq[i:i + 3]
            protein += table[codon]
    return protein

def REbase_org(org): #DONE
    n_new = 0
    n_repeat = 0
    n_unknown = 0
    # function return REbase of microrganism (org):
    #### Key - Recogn. Site as string
    #### Value - {'Enzyme Name':----,'MicroOrganism':-----}

    RE_org_ind = [i for i, k in enumerate(REbase)
                  if org in k] ## index of the 3rd section of the "block"

    rest_dict_org = {}  ## REbase dict.
    for k in RE_org_ind:
        enzym=REbase[k + 2][3:]
        if  enzym =='?' or ',' in enzym:
            n_unknown +=1
        elif enzym in rest_dict_org.keys():
            n_repeat +=1
        #if enzym!='?':
        else:
            n_new +=1
            rest_dict_org[enzym] = {
                'Enzyme Name': REbase[k - 2][3:],
                'MicroOrganism': REbase[k][3:]}
    print(n_new, n_repeat, n_unknown)
    return rest_dict_org

#### returen the NT options of the seq. :
def relevant_seq(re_raw_seq):
    re_raw_seq=re_raw_seq.replace('^','')
    re_raw_seq=re.sub(r'\([^)]*\)', '', re_raw_seq)
    return re_raw_seq


def seq_opps(seq,sa=sa): #dpne
    seq=seq.replace('^','')
    seq=re.sub(r'\([^)]*\)', '', seq)
    sa=defaultdict(list,sa)
    curr_seq=[seq]
    for nt in range(len(seq)):
        if not sa[seq[nt]]:
            curseq=seq[nt]
        elif sa[seq[nt]]:
            new_curr_seqs=[]
            for val in sa[seq[nt]]:
                for s in curr_seq:
                    s = list(s)
                    s[nt]=val
                    s=''.join(s)
                    new_curr_seqs+=[s]
            curr_seq=new_curr_seqs
    return curr_seq

## Comparison:
def comp_seq(seq1,seq2):
# These are written from 5' to 3', only one strand being given.
    seq1=seq_opps(seq1)
    seq2=seq_opps(seq2)
    return re.search(seq2,seq1)
def unique(list1):
    return sorted(set(list1))

def flatten(t):
    return [item for sublist in t for item in sublist]

def unpack_opt(list_lists1):
    n=1
    for sublist in list_lists1:
        if len(sublist):
            n=n*len(sublist)
    list_frame=['']*n
    for elemnt in list_lists1:## 3 options []
        if isinstance(elemnt,list) and (len(elemnt)>1):
            m=0
            if "wind" in locals():
                wind=int(wind/len(elemnt))
            else:
                wind=int(n/len(elemnt))
            for repeat in range(0,int(n/wind)):
                for subelemnt in elemnt: ## aa of the option
                    list_frame[m:m+wind]=[opt+subelemnt for opt in list_frame[m:m+wind]]
                    m=m+wind
        elif isinstance(elemnt,str)& len(elemnt): 
            list_frame=[opt+elemnt for opt in list_frame]
        elif isinstance(elemnt,list)& len(elemnt)==1:
            list_frame=[opt+elemnt[0] for opt in list_frame]
    return list_frame          
                
            


def seq2aaS(recoSeq):
    ''''recoSeq is the recognition sequence 
    the function return the aa possibile sequences for the recoSeq'''
    n=len(recoSeq)
    recoSeq1=recoSeq[:n-n%3]
    recoSeq1rest=recoSeq[n-n%3:]
    endAA1=[v for k, v in table.items() if k.startswith(recoSeq1rest)and (recoSeq1rest!='')]
    aa_options1=unpack_opt(translate(recoSeq1).split()+[unique(endAA1)])
    endNT1=[k for k, v in table.items() if k.startswith(recoSeq1rest)and (recoSeq1rest!='')]
    NT_options1=unpack_opt(recoSeq1.split()+[unique(endNT1)])

    
    recoSeq2 = recoSeq[1:n - (n-1) % 3]
    recoSeq2rest = recoSeq[n - (n-1) % 3:]
    startAA2=[v for k, v in table.items() if k.endswith(recoSeq[0])]
    startNT2=[k for k, v in table.items() if k.endswith(recoSeq[0])]

    endAA2=[v for k, v in table.items() if k.startswith(recoSeq2rest)and (recoSeq2rest!='')]
    aa_options2 = unpack_opt([unique(startAA2)]+translate(recoSeq2).split() + [unique(endAA2)])
    endNT2=[k for k, v in table.items() if k.startswith(recoSeq2rest)and (recoSeq2rest!='')]
    NT_options2 = unpack_opt([unique(startNT2)]+recoSeq2.split() + [unique(endNT2)])

    recoSeq3 = recoSeq[2:n - (n - 2) % 3]
    recoSeq3rest = recoSeq[n - (n - 2) % 3:]
    startAA3 = [v for k, v in table.items() if k.endswith(recoSeq[0:2])]
    startNT3 = [k for k, v in table.items() if k.endswith(recoSeq[0:2])]
    
    endAA3 = [v for k, v in table.items() if k.startswith(recoSeq3rest) and (recoSeq3rest!='')]
    aa_options3 = unpack_opt([unique(startAA3)] + translate(recoSeq3).split() + [unique(endAA3)])
    endNT3 = [k for k, v in table.items() if k.startswith(recoSeq3rest) and (recoSeq3rest!='')]
    NT_options3 = unpack_opt([unique(startNT3)] + recoSeq3.split() + [unique(endNT3)])
    
    aa_options=[aa_options1,aa_options2,aa_options3]
    NT_options=[NT_options1,NT_options2,NT_options3]
    
    return aa_options,NT_options

def ntaa_options_2_ntaa_dict(aa_options,NT_options):
    aa_options=flatten(aa_options)
    NT_options=flatten(NT_options)
    ntaa_dict={}
    for indx in range(len(aa_options)):
        ntaa_dict[aa_options[indx]]=NT_options[indx]
    return ntaa_dict



"""def ntdict_to_aadict(options_ntdict):
    aa_dict = {}
    for seq_id, nt_opt in options_ntdict.items():
        for nt in nt_opt:
            [aa_options,NT_options]=seq2aaS(nt)
            aa_dict[nt] = unique(flatten(aa_options))
    return aa_dict"""


def insert_site_CDS(cds,ntaa_dict):
    """this function insert req.site to cds according to ntaa_dict but it can insert more than site in same index
and insert over what have been done before"""
    AA_cds_seq=translate(cds)
    start_end_idex_dict = {}
    
    for aasite,ntsite in ntaa_dict.items():
        if aasite in AA_cds_seq:
            stratindx=AA_cds_seq.find(aasite)
            endindx=stratindx+len(aasite)
            cds=list(cds)
            cds[3*stratindx:3*endindx]=list(ntsite)
            cds=''.join(cds)
            start_end_idex_dict[ntsite] = {'start': stratindx,
                                      'end':stratindx + len(aasite),
                                      'aa': aasite}
            break
                
                
            
    return cds, start_end_idex_dict

def remove_site_Plasmid(plasmid,ntaa_dict,NT_options):
    AA_plasmid_seq=translate(plasmid)
    plasmid_opt=[]
    rest_plasmid_opt=[]
    start_end_idex_dict_delete = {}
    NT_options=flatten(NT_options)
    for ntsite in NT_options:
        if ntsite in plasmid:
            stratindx=plasmid.find(ntsite)
            endindx=stratindx+len(ntsite)
            plasmid1=list(plasmid)
            plasmid1[stratindx:endindx]='-'*(endindx-stratindx)
            plasmid1=''.join(plasmid1)
            plasmid_opt.append(plasmid1)

            rest_plasmid=['-']*len(plasmid)
            rest_plasmid[stratindx:endindx]=ntsite
            rest_plasmid=''.join(rest_plasmid)
            rest_plasmid_opt.append(rest_plasmid)
    return plasmid_opt,rest_plasmid_opt
    

aa_options,NT_options=seq2aaS('CG')
ntaa_dict=ntaa_options_2_ntaa_dict(aa_options,NT_options)
plasmid='ACGTTGC'
remove_site_Plasmid(plasmid,ntaa_dict,NT_options)

"""def remove_site_Plasmid(cds):
    AA_cds_seq=translate(cds)
    n=0
    for aasite in list1:
        n+=1
        if aasite in AA_cds_seq:
            stratindx=AA_cds_seq.find(aasite)
            endindx=stratindx+len(aasite)
            if cds[3*startindx:3*endindx+1]==ntlist1[n]:
                cds[3*startindx:3*endindx+1]=''
            
    return cds    """




def remove_site_from_plasmid(plasmid_nt_cds, re_nt_list):
    #plasmid must be in the correct reading frame
    site_constraints = [AvoidPattern(nt) for nt in re_nt_list]+[EnforceTranslation()]
    problem = DnaOptimizationProblem(sequence = plasmid_nt_cds, constraints = site_constraints)
    problem.resolve_constraints()
    return problem.sequence

    






