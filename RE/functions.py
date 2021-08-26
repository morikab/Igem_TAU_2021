import re
from collections import defaultdict
from dnachisel import DnaOptimizationProblem, AvoidPattern, EnforceTranslation
import time
import itertools
import string




##open REbase from website : http://rebase.neb.com/rebase/rebase.f12.html
## REbase varible contain every line in the REbase.txt
# """            Every enzyme data called Block

# every block of enzyme is :
# 1. Enzyme Name
# 2. Recognition Site
# 3. Microrganism                       """
f = open(r'REbase_16.8.txt')
content = f.readlines()
REbase = [x.strip() for x in content]


# used variables
ambiguous_code = {'R': ['G', 'A'], 'Y': ['C', 'T'], 'M': ['A', 'C'], 'K': ['G', 'T'],
      'S': ['G', 'C'], 'W': ['A', 'T'], 'B': ['G', 'C', 'T'],
      'D': ['A', 'G', 'T'], 'H': ['A', 'C', 'T'], 'V': ['A', 'G', 'C'],
      'N': ['A', 'G', 'C', 'T'], 'A':['A'], 'G':['G'], 'C':['C'], 'T':['T']
                  }

nt_to_aa = {
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
# --------------------------------------------------------






def unique(list1):
    return sorted(set(list1))


def flatten(t):
    return [item for sublist in t for item in sublist]

def translate(seq, table=nt_to_aa):  # DONE
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
            string += list_of_lists[sublist_idx][str_idex]
        final_str_list.append(string)
    return final_str_list

def write_fasta(fid, list_seq, list_name):
    ofile = open(fid + '.fasta', "w+")
    for i in range(len(list_seq)):
        ofile.write(">" + list_name[i] + "\n" + list_seq[i] + "\n")
    ofile.close()
# --------------------------------------------------------



# creation of an extended_dict- translate ambiguous genetic code
extended_dict = {}
for i0, i0_val in ambiguous_code.items():
    for i1, i1_val in ambiguous_code.items():
        for i2, i2_val in ambiguous_code.items():
            codon_options = itertools_prod_for_strs([i0_val, i1_val, i2_val])
            amino_acid_opts = unique([nt_to_aa[aa] for aa in codon_options])
            extended_dict[i0+i1+i2] = amino_acid_opts
# --------------------------------------------------------




# extracting sequences from
def relevant_seq(re_raw_seq):
    re_raw_seq = re_raw_seq.replace('^', '')
    re_raw_seq = re.sub(r'\([^)]*\)', '', re_raw_seq)
    return re_raw_seq


def unambiguous_seqs(seq, amb_to_nt=ambiguous_code):
    'translates ambigious code into all seq options- str to list of strs, based on the amb_to_nt dictionary'
    list_of_nt_lists = [amb_to_nt[i] for i in seq]
    return itertools_prod_for_strs(list_of_nt_lists)




def translate_ambiguous_3RFs(nt_seq, cds_aa):
    seq_opt_withN = []
    residual = len(nt_seq) % 3
    for i in range(3):
        seq_opt_withN.append('N' * i + nt_seq + 'N' * ((3 - i - residual) % 3))

    seq_to_aa = {}
    for ambiguous_seq in seq_opt_withN:
        for nt_seq in  unambiguous_seqs(ambiguous_seq):
            aa_seq = translate(nt_seq)
            if aa_seq in cds_aa:
                seq_to_aa[nt_seq] = aa_seq
    return seq_to_aa



def REbase_org(org, cds_aa):  # DONE
    RE_org_ind = [i for i, k in enumerate(REbase)
                  if org in k]  ## index of the 3rd section of the "block"
    enzyme_dict = {}  ## REbase dict.
    for k in RE_org_ind:
        site = REbase[k + 2][3:]
        if '?' in site or ',' in site:
            continue
        enzyme_name = REbase[k - 2][3:]
        ambiguous_site = relevant_seq(site)
        if ambiguous_site in [i['ambiguous_site'] for i in enzyme_dict.values()]:
            continue
        nt_to_aa = translate_ambiguous_3RFs(ambiguous_site, cds_aa)
        if nt_to_aa:
            enzyme_dict[enzyme_name] = {
                'ambiguous_site':ambiguous_site,
                'nt_to_aa': nt_to_aa}
    return enzyme_dict




def insert_site_CDS(cds_nt, enzyme_dict):
    """this function insert req.site to cds according to ntaa_dict but it can insert more than site in same index
and insert over what have been done before"""
    cds_aa = translate(cds_nt)
    # start_end_idex_dict = {}
    for sub_dict in enzyme_dict.values():
        ntaa_dict = sub_dict['nt_to_aa']
        for aasite, ntsites in ntaa_dict.items():
            if aasite in cds_aa:
                stratindx = cds_aa.find(aasite)
                endindx = stratindx + len(aasite)
                cds_nt = list(cds_nt)
                cds_nt[3 * stratindx:3 * endindx] = list(ntsites[0])
                cds_nt = ''.join(cds_nt)
                # start_end_idex_dict[ntsites[0]] = {'start': stratindx,
                #                                    'end': stratindx + len(aasite),
                #                                    'aa': aasite}
    return cds_nt


def remove_site_from_plasmid(cds_nt, enzyme_dict):
    # plasmid must be in the correct reading frame
    cds_nt = str(cds_nt)
    tested_re_nt_list = []
    for sub_dict in enzyme_dict.values():
        tested_re_nt_list = tested_re_nt_list + list(sub_dict['nt_to_aa'].keys())
    site_constraints = [AvoidPattern(nt) for nt in tested_re_nt_list] + [EnforceTranslation()]
    problem = DnaOptimizationProblem(sequence=cds_nt, constraints=site_constraints)
    problem.resolve_constraints()
    return problem.sequence


def sites_in_cds(enzyme_dict, cds_nt):
    found_sites = {}
    for enzyme_name, enzyme_val in enzyme_dict.items():
        if [i for i in  list(enzyme_val['nt_to_aa'].keys()) if i in cds_nt]:
            found_sites[enzyme_name] = enzyme_val['ambiguous_site']
    return found_sites


def parse_inp1(usr_inp1):
    optimized_org_names = []
    deoptimized_org_names = []
    for org_name, org_dict in usr_inp1.items():
        if org_dict['optimized']:
            optimized_org_names.append(org_name)
        else:
            deoptimized_org_names.append(org_name)
    return optimized_org_names, deoptimized_org_names


