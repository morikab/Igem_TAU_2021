import re
from dnachisel import DnaOptimizationProblem, AvoidPattern, EnforceTranslation
import itertools
import os
from modules.shared_functions_and_vars import *


##open REbase from website : http://rebase.neb.com/rebase/rebase.f12.html
## REbase varible contain every line in the REbase.txt
# """            Every enzyme data called Block

# every block of enzyme is :
# 1. Enzyme Name
# 2. Recognition Site
# 3. Microrganism                       """

base_path = os.path.join(os.path.dirname(__file__))

f = open(os.path.join(base_path, 'REbase_16.8.txt'))
content = f.readlines()
REbase = [x.strip() for x in content]


# --------------------------------------------------------









def flatten(t):
    return [item for sublist in t for item in sublist]


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
    'removing irrelevant chars from RE sites from the REDB'
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
    '''
    parsing the REbase and creating a dictionary with the relevant sites from
    enzymes originated in org
    :param org: scientific name of investigated organism
    :param cds_aa: aa sequence of the cds
    :return: a dict in the following format:
                enzyme_dict[enzyme_name] = {
                'ambiguous_site':ambiguous_site, #contaiinng characters from the ambiguous_code dictionary
                'nt_to_aa': nt_to_aa} #non ambigiosus sites translated (in all 3 RFs)
    '''
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




def insert_site_CDS(RE_dict, NTseq):
    """inpout: -  RE_dict=dict. of all {organism:{siteName:{'ambiguous_site':NTsite,'nt_to_aa':{}}}}.
               -  NTseq= string of the sequence to insert in it. e.g. fullPlasmid,Promoter,CDS...
       output: """
    NTseq = translate(AAseq)
    busy=[]
    NTsite_org_inSeq={}
    overlappedSitesNT=[]
    for orgName,orginfo in RE_dict.items():
        for siteName,siteseq in orginfo.items():
            NT2AA=siteseq['nt_to_aa']
            for NTsite,AAsite in NT2AA.items():
                for sup in list(re.finditer(,NTseq)):
                    findedindx=[*range(sup.start(),sup.end())]
                    for indx in findedindex:
                        if indx not in busy:
                            busy.append(indx)
                        else:
                            if indx not in overlapped_indx:
                                overlapped_indx[indx]={}
                            overlapped_indx[indx][siteName]+=[AAsite]
                            overlappedSitesNT+=[NTsite]
                    NTsite_org_inSeq[NTsite][orgName]+=[NTsite_org_inSeq[NTsite][orgName],[sup.start(),sup.end()]]

        NTsite_org_inSeq[NTsite][orgName]=

    for NTsite,siteOrgInfo in NTsite_org_inSeq.items():
        if NTsite not in overlappedSitesNT:
            stratindx=siteOrgInfo
            cds_nt = list(cds_nt)
            cds_nt[3 * stratindx:3 * endindx] = list(ntsites[0])
            cds_nt = ''.join(cds_nt)
            
                            
                                


        site_org_inSeq[NTsite][
                    

    
    startEnd_sitesinSeq={} # {Organism[start end start end ....... ] 
    
    
    for sub_dict in enzyme_dict.values():
        ntaa_dict = sub_dict['nt_to_aa']
        for aasite, ntsites in ntaa_dict.items():
            if aasite in cds_aa:
                stratindx = cds_aa.find(aasite)
                endindx = stratindx + len(aasite)

                ##INSERT
##                cds_nt = list(cds_nt)
##                cds_nt[3 * stratindx:3 * endindx] = list(ntsites[0])
##                cds_nt = ''.join(cds_nt)
    return cds_nt


def remove_site_from_plasmid(cds_nt, enzyme_dict):
    'plasmid must be in the correct reading frame'
    cds_nt = str(cds_nt)
    tested_re_nt_list = []
    for sub_dict in enzyme_dict.values():
        tested_re_nt_list = tested_re_nt_list + list(sub_dict['nt_to_aa'].keys())
    site_constraints = [AvoidPattern(nt) for nt in tested_re_nt_list] + [EnforceTranslation()]
    problem = DnaOptimizationProblem(sequence=cds_nt, constraints=site_constraints)
    problem.resolve_constraints()
    return problem.sequence


def sites_in_cds(enzyme_dict, cds_nt):
    '''
    find out if and what sites are found in the cds_nt
    :param enzyme_dict: REbase_org output dict
    :param cds_nt: the str of the final sequence
    :return: {enzyme_name:ambiguous_site }
    '''
    found_sites_dict = {}
    for enzyme_name, enzyme_val in enzyme_dict.items():
        if [i for i in  list(enzyme_val['nt_to_aa'].keys()) if i in cds_nt]:
            found_sites_dict[enzyme_name] = enzyme_val['ambiguous_site']
    return found_sites_dict





