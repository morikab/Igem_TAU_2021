import re
from collections import defaultdict
from dnachisel import DnaOptimizationProblem, AvoidPattern,EnforceTranslation, EnforcePatternOccurence
import time
import itertools
import string

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

table_extended = {'AAA': ['K'], 'AAG': ['K'], 'AAC': ['N'], 'AAT': ['N'], 'AAN': ['K', 'N'], 'AGA': ['R'], 'AGG': ['R'],
                  'AGC': ['S'], 'AGT': ['S'], 'AGN': ['R', 'S'], 'ACA': ['T'], 'ACG': ['T'], 'ACC': ['T'], 'ACT': ['T'],
                  'ACN': ['T'], 'ATA': ['I'], 'ATG': ['M'], 'ATC': ['I'], 'ATT': ['I'], 'ATN': ['I', 'M'],
                  'ANA': ['I', 'K', 'R', 'T'], 'ANG': ['K', 'M', 'R', 'T'], 'ANC': ['I', 'N', 'S', 'T'],
                  'ANT': ['I', 'N', 'S', 'T'], 'ANN': ['I', 'K', 'M', 'N', 'R', 'S', 'T'], 'GAA': ['E'], 'GAG': ['E'],
                  'GAC': ['D'], 'GAT': ['D'], 'GAN': ['D', 'E'], 'GGA': ['G'], 'GGG': ['G'], 'GGC': ['G'], 'GGT': ['G'],
                  'GGN': ['G'], 'GCA': ['A'], 'GCG': ['A'], 'GCC': ['A'], 'GCT': ['A'], 'GCN': ['A'], 'GTA': ['V'],
                  'GTG': ['V'], 'GTC': ['V'], 'GTT': ['V'], 'GTN': ['V'], 'GNA': ['A', 'E', 'G', 'V'],
                  'GNG': ['A', 'E', 'G', 'V'], 'GNC': ['A', 'D', 'G', 'V'], 'GNT': ['A', 'D', 'G', 'V'],
                  'GNN': ['A', 'D', 'E', 'G', 'V'], 'CAA': ['Q'], 'CAG': ['Q'], 'CAC': ['H'], 'CAT': ['H'],
                  'CAN': ['H', 'Q'], 'CGA': ['R'], 'CGG': ['R'], 'CGC': ['R'], 'CGT': ['R'], 'CGN': ['R'], 'CCA': ['P'],
                  'CCG': ['P'], 'CCC': ['P'], 'CCT': ['P'], 'CCN': ['P'], 'CTA': ['L'], 'CTG': ['L'], 'CTC': ['L'],
                  'CTT': ['L'], 'CTN': ['L'], 'CNA': ['L', 'P', 'Q', 'R'], 'CNG': ['L', 'P', 'Q', 'R'],
                  'CNC': ['H', 'L', 'P', 'R'], 'CNT': ['H', 'L', 'P', 'R'], 'CNN': ['H', 'L', 'P', 'Q', 'R'],
                  'TAA': ['_'], 'TAG': ['_'], 'TAC': ['Y'], 'TAT': ['Y'], 'TAN': ['Y', '_'], 'TGA': ['_'], 'TGG': ['W'],
                  'TGC': ['C'], 'TGT': ['C'], 'TGN': ['C', 'W', '_'], 'TCA': ['S'], 'TCG': ['S'], 'TCC': ['S'],
                  'TCT': ['S'], 'TCN': ['S'], 'TTA': ['L'], 'TTG': ['L'], 'TTC': ['F'], 'TTT': ['F'], 'TTN': ['F', 'L'],
                  'TNA': ['L', 'S', '_'], 'TNG': ['L', 'S', 'W', '_'], 'TNC': ['C', 'F', 'S', 'Y'],
                  'TNT': ['C', 'F', 'S', 'Y'], 'TNN': ['C', 'F', 'L', 'S', 'W', 'Y', '_'], 'NAA': ['E', 'K', 'Q', '_'],
                  'NAG': ['E', 'K', 'Q', '_'], 'NAC': ['D', 'H', 'N', 'Y'], 'NAT': ['D', 'H', 'N', 'Y'],
                  'NAN': ['D', 'E', 'H', 'K', 'N', 'Q', 'Y', '_'], 'NGA': ['G', 'R', '_'], 'NGG': ['G', 'R', 'W'],
                  'NGC': ['C', 'G', 'R', 'S'], 'NGT': ['C', 'G', 'R', 'S'], 'NGN': ['C', 'G', 'R', 'S', 'W', '_'],
                  'NCA': ['A', 'P', 'S', 'T'], 'NCG': ['A', 'P', 'S', 'T'], 'NCC': ['A', 'P', 'S', 'T'],
                  'NCT': ['A', 'P', 'S', 'T'], 'NCN': ['A', 'P', 'S', 'T'], 'NTA': ['I', 'L', 'V'], 'NTG': ['L', 'M', 'V'],
                  'NTC': ['F', 'I', 'L', 'V'], 'NTT': ['F', 'I', 'L', 'V'], 'NTN': ['F', 'I', 'L', 'M', 'V'],
                  'NNA': ['A', '(E', 'G', 'I', 'K', 'L', 'P', 'Q', 'R', 'S', 'T', 'V', '_'],
                  'NNG': ['A', 'E', 'G', 'K', 'L', 'M', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', '_'],
                  'NNC': ['A', 'C', 'D', 'F', 'G', 'H', 'I', 'L', 'N', 'P', 'R', 'S', 'T', 'V', 'Y'],
                  'NNT': ['A', 'C', 'D', 'F', 'G', 'H', 'I', 'L', 'N', 'P', 'R', 'S', 'T', 'V', 'Y']}


def translate(seq,table=table): #DONE
    protein = ""
    if len(seq) % 3 == 0:
        for i in range(0, len(seq), 3):
            codon = seq[i:i + 3]
            protein += table[codon]
        return protein
    else:
        return ValueError('len(seq)%3 !=0')


def itertools_prod_for_strs(list_of_lists):
    '''
    :param list_of_lists: ex:[['a', 'b'], ['c', 'd', 'g'], ['e', 'f']]
    :return: ex:['ace', 'acf', 'ade', 'adf', 'age', 'agf', 'bce', 'bcf', 'bde', 'bdf', 'bge', 'bgf']
    '''
    list_of_idxlist = [list(range(len(i))) for i in list_of_lists]
    itertools_prod = itertools.product(*list_of_idxlist)

    final_str_list = []
    for idx_tup in itertools_prod:
        string = ""
        for sublist_idx, str_idex in enumerate(idx_tup):
            string +=list_of_lists[sublist_idx][str_idex]
        final_str_list.append(string)
    return final_str_list





# def str_itertools_product():
#     for item in itertools.product(string.ascii_lowercase, repeat=2):
#         yield "".join(list(item))


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
                

def translate_withN(seq, cds_aa, used_table=table_extended):
    seq_opt_withN = []
    residual = len(seq)%3
    for i in range(3):
        seq_opt_withN.append('N'*i + seq + 'N'*((3-i-residual)%3))

    seq_aa = {}
    for seq in seq_opt_withN:
        opt_list = []
        for i in range(0, len(seq), 3):
            aa_opts = used_table[seq[i:i + 3]]
            opt_list.append(aa_opts)
        possible_seqs = [i for i in itertools_prod_for_strs(opt_list) if i in cds_aa]

        if len(possible_seqs)>0:
            seq_aa[seq] = itertools_prod_for_strs(opt_list)
    return seq_aa



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
    print(seq_opts_list)

    for seq_opt in seq_opts_list:
        ntaa_dict =translate_withN(seq_opt, plasmid_aa_cds, used_table=table_extended)
        if len (ntaa_dict) == 0 :
            continue
        org_NT_options=org_NT_options+list(ntaa_dict.keys())
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




#
# def codon_options(codon, table=table):
#     codon_options = []
#     aa_list = []
#     for i in codon:
#         if i in ['A', 'G', 'C', 'T']:
#             codon_options.append([i])
#         else:
#             codon_options.append(['A', 'G', 'C', 'T'])
#     for i1 in codon_options[0]:
#         for i2 in codon_options[1]:
#             for i3 in codon_options[2]:
#                 aa_list.append(table[i1 + i2 + i3])
#     return unique(aa_list)
#
#
# all_codons = ['AAA', 'AAG', 'AAC', 'AAT', 'AAN', 'AGA', 'AGG', 'AGC', 'AGT', 'AGN', 'ACA', 'ACG', 'ACC', 'ACT', 'ACN', 'ATA', 'ATG', 'ATC', 'ATT', 'ATN', 'ANA', 'ANG', 'ANC', 'ANT', 'ANN', 'GAA', 'GAG', 'GAC', 'GAT', 'GAN', 'GGA', 'GGG', 'GGC', 'GGT', 'GGN', 'GCA', 'GCG', 'GCC', 'GCT', 'GCN', 'GTA', 'GTG', 'GTC', 'GTT', 'GTN', 'GNA', 'GNG', 'GNC', 'GNT', 'GNN', 'CAA', 'CAG', 'CAC', 'CAT', 'CAN', 'CGA', 'CGG', 'CGC', 'CGT', 'CGN', 'CCA', 'CCG', 'CCC', 'CCT', 'CCN', 'CTA', 'CTG', 'CTC', 'CTT', 'CTN', 'CNA', 'CNG', 'CNC', 'CNT', 'CNN', 'TAA', 'TAG', 'TAC', 'TAT', 'TAN', 'TGA', 'TGG', 'TGC', 'TGT', 'TGN', 'TCA', 'TCG', 'TCC', 'TCT', 'TCN', 'TTA', 'TTG', 'TTC', 'TTT', 'TTN', 'TNA', 'TNG', 'TNC', 'TNT', 'TNN', 'NAA', 'NAG', 'NAC', 'NAT', 'NAN', 'NGA', 'NGG', 'NGC', 'NGT', 'NGN', 'NCA', 'NCG', 'NCC', 'NCT', 'NCN', 'NTA', 'NTG', 'NTC', 'NTT', 'NTN', 'NNA', 'NNG', 'NNC', 'NNT', 'NNN']
# all_codon_dict = {}
# for codon in all_codons:
#     all_codon_dict[codon] = codon_options(codon)
#
# print(all_codon_dict)


