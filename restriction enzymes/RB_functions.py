import re
from collections import defaultdict
from dnachisel import DnaOptimizationProblem, AvoidPattern,EnforceTranslation
import time

def write_fasta(fid, list_seq, list_name):
    ofile = open(fid + '.fasta', "w+")
    for i in range(len(list_seq)):
        ofile.write(">" + list_name[i] + "\n" + list_seq[i] + "\n")
    ofile.close()


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
    else:
        return ValueError('len(seq)%3 !=0')



def REbase_org(org): #DONE
    RE_org_ind = [i for i, k in enumerate(REbase)
                  if org in k] ## index of the 3rd section of the "block"
    rest_dict_org = {}  ## REbase dict.
    for k in RE_org_ind:
        site=REbase[k + 2][3:]
        if  '?' in site or ',' in site:
            continue
        rest_dict_org[site] = {
            'Enzyme Name': REbase[k - 2][3:],
            'MicroOrganism': REbase[k][3:]}
    return rest_dict_org

def relevant_seq(re_raw_seq):
    re_raw_seq=re_raw_seq.replace('^','')
    re_raw_seq=re.sub(r'\([^)]*\)', '', re_raw_seq)
    return re_raw_seq


def seq_opps(seq,sa=sa): #dpne
    seq=relevant_seq(seq)
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
                



def seq2aaS(recoSeq, plasmid_aa_cds):
    ''''recoSeq is the recognition sequence 
    the function return the aa possibile sequences for the recoSeq'''
    n=len(recoSeq)
    nuc_list = ['A', 'G', 'C', 'T' ]
    if n%3 ==0:
        NT_options = [recoSeq]
    elif n%3 ==1:
        double_nuc_list = []
        for i in nuc_list:
            for k in nuc_list:
                double_nuc_list.append(k + i)
        NT_options = [recoSeq+k for k in double_nuc_list ]+ [k+recoSeq for k in double_nuc_list]
    else:
        NT_options = [recoSeq+k for k in nuc_list ]+ [k+recoSeq for k in nuc_list]


    print(1, NT_options)
    NT_options_final = []
    ntaa_dict={}
    for nt_seq in NT_options:
        aa_seq=translate(nt_seq)
        print(2, aa_seq)
        if aa_seq in plasmid_aa_cds:
            print(1)
            if aa_seq not in ntaa_dict:
                ntaa_dict[aa_seq]=[nt_seq]
            else:
                ntaa_dict[aa_seq].append(nt_seq)
            NT_options_final.append(nt_seq)
    return ntaa_dict,NT_options_final


def REseq_org(org, plasmid_nt_cds):
    plasmid_aa_cds = translate(plasmid_nt_cds)
    org_dict=REbase_org(org) #{'site': ['name', microorganism']}
    org_seq_list=[*org_dict.keys()] #list of sites
    org_seq_opts={}
    org_NT_options=[]
    seq_opts_list = []

    for seq in org_seq_list:
        seq_opts_list =seq_opts_list + seq_opps(seq) #all options according to SA
    seq_opts_list = unique(seq_opts_list)

    for seq_opt in seq_opts_list:
        ntaa_dict,NT_options=seq2aaS(seq_opt, plasmid_aa_cds)
        if len (ntaa_dict) == 0 :
            continue
        org_NT_options=org_NT_options+NT_options
        print(org_NT_options)
        for aa_site,nt_site in ntaa_dict.items():
            if aa_site not in org_seq_opts:
                org_seq_opts[aa_site]=nt_site
            else:
                org_seq_opts[aa_site]=org_seq_opts[aa_site]+nt_site

    return org_seq_opts,org_NT_options

        
def insert_site_CDS(cds,ntaa_dict):
    """this function insert req.site to cds according to ntaa_dict but it can insert more than site in same index
and insert over what have been done before"""
    AA_cds_seq=translate(cds)
    start_end_idex_dict = {}
    for aasite,ntsites in ntaa_dict.items():
        if aasite in AA_cds_seq:
            stratindx=AA_cds_seq.find(aasite)
            endindx=stratindx+len(aasite)
            cds=list(cds)
            cds[3*stratindx:3*endindx]=list(ntsites[0])
            cds=''.join(cds)
            start_end_idex_dict[ntsites[0]] = {'start': stratindx,
                                      'end':stratindx + len(aasite),
                                      'aa': aasite}
            break
    return cds, start_end_idex_dict

def remove_site_from_plasmid(plasmid_nt_cds, re_nt_list):
    #plasmid must be in the correct reading frame
    plasmid_nt_cds = str(plasmid_nt_cds)
    tested_re_nt_list = [nt for nt in re_nt_list if nt in plasmid_nt_cds]
    site_constraints = [AvoidPattern(nt) for nt in tested_re_nt_list]+[EnforceTranslation()]
    problem = DnaOptimizationProblem(sequence = plasmid_nt_cds, constraints = site_constraints)
    problem.resolve_constraints()
    return problem.sequence
