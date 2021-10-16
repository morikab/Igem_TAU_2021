import re
from dnachisel import DnaOptimizationProblem, AvoidPattern, EnforceTranslation
import itertools
import os
from modules.shared_functions_and_vars import *
import numpy as np

##open REbase from website : http://rebase.neb.com/rebase/rebase.f12.html
## REbase varible contain every line in the REbase.txt
# """            Every enzyme data called Block

# every block of enzyme is :
# 1. Enzyme Name
# 2. Recognition Site
# 3. Microrganism                       """

base_path = os.path.join(os.path.dirname(__file__))

f = open(os.path.join(base_path, 'REbase.txt'))
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

inputSeqNT='GGGCCCGGGCCCTATGGCTAA'
RE_dict={'Bacillus subtilis': {'Bsu1192I': {'ambiguous_site': 'CCGG', 'nt_to_aa': {'CCGGAA': 'PE', 'CCGGAG': 'PE', 'CCGGAC': 'PD', 'CCGGAT': 'PD', 'CCGGGA': 'PG', 'CCGGGG': 'PG', 'CCGGGC': 'PG', 'CCGGGT': 'PG', 'CCGGCA': 'PA', 'CCGGCG': 'PA', 'CCGGCC': 'PA', 'CCGGCT': 'PA', 'CCGGTA': 'PV', 'CCGGTG': 'PV', 'CCGGTC': 'PV', 'CCGGTT': 'PV', 'ACCGGA': 'TG', 'ACCGGG': 'TG', 'ACCGGC': 'TG', 'ACCGGT': 'TG', 'CCCGGA': 'PG', 'CCCGGG': 'PG', 'CCCGGC': 'PG', 'CCCGGT': 'PG', 'GGCCGG': 'GR', 'CTCCGG': 'LR'}}, 'Bsu1192II': {'ambiguous_site': 'CGCG', 'nt_to_aa': {'CGCGGA': 'RG', 'CGCGGG': 'RG', 'CGCGGC': 'RG', 'CGCGGT': 'RG', 'CGCGCA': 'RA', 'CGCGCG': 'RA', 'CGCGCC': 'RA', 'CGCGCT': 'RA', 'CGCGTA': 'RV', 'CGCGTG': 'RV', 'CGCGTC': 'RV', 'CGCGTT': 'RV', 'ACCGCG': 'TA', 'GACGCG': 'DA', 'GGCGCG': 'GA', 'CCCGCG': 'PA', 'TTCGCG': 'FA'}}, 'Bsu8565I': {'ambiguous_site': 'GGATCC', 'nt_to_aa': {'GGATCC': 'GS'}}, 'BsuBI': {'ambiguous_site': 'CTGCAG', 'nt_to_aa': {'CTGCAG': 'LQ', 'CCTGCAGAC': 'PAD', 'CCTGCAGAT': 'PAD'}}, 'BsuMI': {'ambiguous_site': 'CTCGAG', 'nt_to_aa': {'GACTCGAGC': 'DSS', 'GACTCGAGT': 'DSS', 'GCCTCGAGC': 'ASS', 'GCCTCGAGT': 'ASS'}},
'Bsu1854I': {'ambiguous_site': 'GRGCYC', 'nt_to_aa': {'GGGCCC': 'GP', 'GAGCTC': 'EL', 'GGGGCCCTA': 'GAL', 'GGGGCCCTG': 'GAL', 'GGGGCCCTC': 'GAL', 'GGGGCCCTT': 'GAL', 'GGGGCTCTA': 'GAL', 'GGGGCTCTG': 'GAL', 'GGGGCTCTC': 'GAL', 'GGGGCTCTT': 'GAL', 'GGAGCCCTA': 'GAL', 'GGAGCCCTG': 'GAL', 'GGAGCCCTC': 'GAL', 'GGAGCCCTT': 'GAL', 'GGAGCTCTA': 'GAL', 'GGAGCTCTG': 'GAL', 'GGAGCTCTC': 'GAL', 'GGAGCTCTT': 'GAL', 'GAGGGCTCA': 'EGS', 'GAGGGCTCG': 'EGS', 'GAGGGCTCC': 'EGS', 'GAGGGCTCT': 'EGS', 'GGGGGCCCA': 'GGP', 'GGGGGCCCG': 'GGP', 'GGGGGCCCC': 'GGP', 'GGGGGCCCT': 'GGP', 'GCGAGCTCA': 'ASS', 'GCGAGCTCG': 'ASS', 'GCGAGCTCC': 'ASS', 'GCGAGCTCT': 'ASS', 'CTGAGCCCA': 'LSP', 'CTGAGCCCG': 'LSP', 'CTGAGCCCC': 'LSP', 'CTGAGCCCT': 'LSP', 'TTGAGCCCA': 'LSP', 'TTGAGCCCG': 'LSP', 'TTGAGCCCC': 'LSP', 'TTGAGCCCT': 'LSP'}}, 'Bsu22I': {'ambiguous_site': 'TCCGGA', 'nt_to_aa': {'ATTCCGGAC': 'IPD', 'ATTCCGGAT': 'IPD', 'TATCCGGAA': 'YPE', 'TATCCGGAG': 'YPE', 'TTTCCGGAA': 'FPE', 'TTTCCGGAG': 'FPE'}}, 'Bsu36I': {'ambiguous_site': 'CCTNAGG', 'nt_to_aa': {'CCTGAGGAC': 'PED', 'CCTGAGGAT': 'PED', 'CCTGAGGGA': 'PEG', 'CCTGAGGGG': 'PEG', 'CCTGAGGGC': 'PEG', 'CCTGAGGGT': 'PEG'}},
'Bsu537I': {'ambiguous_site': 'GGTCTC', 'nt_to_aa': {'ATGGTCTCA': 'MVS', 'ATGGTCTCG': 'MVS', 'ATGGTCTCC': 'MVS', 'ATGGTCTCT': 'MVS'}}, 'Bsu54I': {'ambiguous_site': 'GGNCC', 'nt_to_aa': {'GGACCA': 'GP', 'GGACCG': 'GP', 'GGACCC': 'GP', 'GGACCT': 'GP', 'GGGCCA': 'GP', 'GGGCCG': 'GP', 'GGGCCC': 'GP', 'GGGCCT': 'GP', 'GGCCCA': 'GP', 'GGCCCG': 'GP', 'GGCCCC': 'GP', 'GGCCCT': 'GP', 'GGTCCA': 'GP', 'GGTCCG': 'GP', 'GGTCCC': 'GP', 'GGTCCT': 'GP', 'AGGGCC': 'RA', 'AGGCCC': 'RP', 'GGGACC': 'GT', 'GGGGCC': 'GA', 'GGGTCC': 'GS', 'CGGGCC': 'RA', 'CGGCCC': 'RP', 'AAGGTCCAC': 'KVH', 'AAGGTCCAT': 'KVH', 'GAGGGCCGA': 'EGR', 'GAGGGCCGG': 'EGR', 'GAGGGCCGC': 'EGR', 'GAGGGCCGT': 'EGR', 'GGGGGCCAC': 'GGH', 'GGGGGCCAT': 'GGH', 'GGGGGCCCA': 'GGP', 'GGGGGCCCG': 'GGP', 'GGGGGCCCC': 'GGP', 'GGGGGCCCT': 'GGP', 'GGGGCCCTA': 'GAL', 'GGGGCCCTG': 'GAL', 'GGGGCCCTC': 'GAL', 'GGGGCCCTT': 'GAL', 'CCGGTCCAA': 'PVQ', 'CCGGTCCAG': 'PVQ'}}, 'Bsu6I': {'ambiguous_site': 'CTCTTC', 'nt_to_aa': {'CCTCTTCCA': 'PLP', 'CCTCTTCCG': 'PLP', 'CCTCTTCCC': 'PLP', 'CCTCTTCCT': 'PLP', 'TCTCTTCAA': 'SLQ', 'TCTCTTCAG': 'SLQ', 'GACTCTTCA': 'DSS', 'GACTCTTCG': 'DSS', 'GACTCTTCC': 'DSS', 'GACTCTTCT': 'DSS', 'GCCTCTTCA': 'ASS', 'GCCTCTTCG': 'ASS', 'GCCTCTTCC': 'ASS', 'GCCTCTTCT': 'ASS'}}, 'Bsu7003I': {'ambiguous_site': 'GACGAG', 'nt_to_aa': {'GACGAG': 'DE'}}, 'Bsu1076I': {'ambiguous_site': 'GGCC', 'nt_to_aa': {'GGCCAC': 'GH', 'GGCCAT': 'GH', 'GGCCGA': 'GR', 'GGCCGG': 'GR', 'GGCCGC': 'GR', 'GGCCGT': 'GR', 'GGCCCA': 'GP', 'GGCCCG': 'GP', 'GGCCCC': 'GP', 'GGCCCT': 'GP', 'AGGCCA': 'RP', 'AGGCCG': 'RP', 'AGGCCC': 'RP', 'AGGCCT': 'RP', 'GGGCCA': 'GP', 'GGGCCG': 'GP', 'GGGCCC': 'GP', 'GGGCCT': 'GP', 'CGGCCA': 'RP', 'CGGCCG': 'RP', 'CGGCCC': 'RP', 'CGGCCT': 'RP', 'AAGGCC': 'KA', 'AGGGCC': 'RA', 'ACGGCC': 'TA', 'ATGGCC': 'MA', 'GAGGCC': 'EA', 'GGGGCC': 'GA', 'CGGGCC': 'RA', 'CCGGCC': 'PA'}}, 'M.Phi3TII': {'ambiguous_site': 'TCGA', 'nt_to_aa': {'TCGAAA': 'SK', 'TCGAAG': 'SK', 'TCGAGC': 'SS', 'TCGAGT': 'SS', 'TCGACA': 'ST', 'TCGACG': 'ST', 'TCGACC': 'ST', 'TCGACT': 'ST', 'ATCGAA': 'IE', 'ATCGAG': 'IE', 'GTCGAA': 'VE', 'GTCGAG': 'VE', 'CTCGAC': 'LD', 'CTCGAT': 'LD', 'TTCGAA': 'FE', 'TTCGAG': 'FE', 'GGTCGA': 'GR', 'CTTCGA': 'LR'}}, 'BisI': {'ambiguous_site': 'GCNGC', 'nt_to_aa': {'AGCAGC': 'SS', 'GGCAGC': 'GS', 'GGCGGC': 'GG', 'GGCCGC': 'GR', 'CGCGGC': 'RG', 'AAGCAGCGA': 'KQR', 'AAGCAGCGG': 'KQR', 'AAGCAGCGC': 'KQR', 'AAGCAGCGT': 'KQR', 'AAGCTGCGA': 'KLR', 'AAGCTGCGG': 'KLR', 'AAGCTGCGC': 'KLR', 'AAGCTGCGT': 'KLR', 'GGGCGGCAC': 'GRH', 'GGGCGGCAT': 'GRH', 'GGGCGGCCA': 'GRP', 'GGGCGGCCG': 'GRP', 'GGGCGGCCC': 'GRP', 'GGGCGGCCT': 'GRP', 'GGGCCGCTA': 'GPL', 'GGGCCGCTG': 'GPL', 'GGGCCGCTC': 'GPL', 'GGGCCGCTT': 'GPL', 'GTGCAGCTA': 'VQL', 'GTGCAGCTG': 'VQL', 'GTGCAGCTC': 'VQL', 'GTGCAGCTT': 'VQL', 'CAGCGGCTA': 'QRL', 'CAGCGGCTG': 'QRL', 'CAGCGGCTC': 'QRL', 'CAGCGGCTT': 'QRL', 'CAGCTGCCA': 'QLP', 'CAGCTGCCG': 'QLP', 'CAGCTGCCC': 'QLP', 'CAGCTGCCT': 'QLP', 'CCGCTGCCA': 'PLP', 'CCGCTGCCG': 'PLP', 'CCGCTGCCC': 'PLP', 'CCGCTGCCT': 'PLP', 'TCGCCGCAA': 'SPQ', 'TCGCCGCAG': 'SPQ', 'TCGCTGCAA': 'SLQ', 'TCGCTGCAG': 'SLQ'}}}, 'Sulfolobus acidocaldarius': {'SulI': {'ambiguous_site': 'GGCC', 'nt_to_aa': {'GGCCAC': 'GH', 'GGCCAT': 'GH', 'GGCCGA': 'GR', 'GGCCGG': 'GR', 'GGCCGC': 'GR', 'GGCCGT': 'GR', 'GGCCCA': 'GP', 'GGCCCG': 'GP', 'GGCCCC': 'GP', 'GGCCCT': 'GP', 'AGGCCA': 'RP', 'AGGCCG': 'RP', 'AGGCCC': 'RP', 'AGGCCT': 'RP', 'GGGCCA': 'GP', 'GGGCCG': 'GP', 'GGGCCC': 'GP', 'GGGCCT': 'GP', 'CGGCCA': 'RP', 'CGGCCG': 'RP', 'CGGCCC': 'RP', 'CGGCCT': 'RP', 'AAGGCC': 'KA', 'AGGGCC': 'RA', 'ACGGCC': 'TA', 'ATGGCC': 'MA', 'GAGGCC': 'EA', 'GGGGCC': 'GA', 'CGGGCC': 'RA', 'CCGGCC': 'PA'}}}}


def insert_site_CDS(RE_dict, inputSeqNT):
    """inpout: -  RE_dict=dict. of all {organism:{siteName:{'ambiguous_site':NTsite,'nt_to_aa':{}}}}.
               -  inputSeqNT= string of the sequence to insert in it. e.g. fullPlasmid,Promoter,CDS...
       output: """
    inputSeqAA = translate(inputSeqNT)
    busyNTseq=[]
    overlapped_NTindx={}
    NTsite_org_inSeq={}
    overlappedSitesNT=[]
    for orgName,orginfo in RE_dict.items():
        for siteName,siteinfo in orginfo.items():
            NT2AA=siteinfo['nt_to_aa']
            for NTsite,AAsite in NT2AA.items():
                for sup in list(re.finditer(AAsite,inputSeqAA)):
                    findedindx=[*range(sup.start()*3,sup.end()*3)]
                    for indx in findedindx:
                        if indx not in busyNTseq:
                            busyNTseq.append(indx)
                        else: # inxd in busyNTseq
                            if indx not in overlapped_NTindx:
                                overlapped_NTindx[indx]={}
                            if orgName not in overlapped_NTindx[indx]:
                                overlapped_NTindx[indx][orgName] = []
                            overlapped_NTindx[indx][orgName]+=[NTsite]
                            overlappedSitesNT+=[NTsite]
                    if NTsite not in NTsite_org_inSeq:
                        NTsite_org_inSeq[NTsite]={}
                    if orgName not in NTsite_org_inSeq[NTsite]:
                        NTsite_org_inSeq[NTsite][orgName]={}
                    if siteName not in NTsite_org_inSeq[NTsite][orgName]:
                        NTsite_org_inSeq[NTsite][orgName][siteName]=[]
                    NTsite_org_inSeq[NTsite][orgName][siteName]+=[(sup.start()*3,sup.end()*3)]

    orgInserted = {}
    outputInsert={}
    outputSeqNT=list(inputSeqNT)
    for NTsite,siteOrgInfo in NTsite_org_inSeq.items():
        if NTsite not in overlappedSitesNT:
            for orgName,siteInfo in siteOrgInfo.items():
                for siteName,positions in siteInfo.items():
                    for position in positions:
                        findedindx = [*range(position[0], position[1])]
                        outputSeqNT[position[0]:position[1]]=list(NTsite)
                        for finded in findedindx:
                            if finded in busyNTseq:
                                busyNTseq.remove(finded)
                    if siteName not in outputInsert:
                        outputInsert[siteName]=[]
                    outputInsert[siteName]+=[positions]
                    if orgName not in orgInserted:
                        orgInserted[orgName]=0
                    orgInserted[orgName]+=1

    for indx,sitesConflicts in overlapped_NTindx.items():
        if indx in busyNTseq:
            indxOrganisms=[*sitesConflicts]
            indxNorg=[sitesConflicts[orgName] for orgName in sitesConflicts]
            maxiOrg=[i for i,x in enumerate(indxNorg) if x==max(indxNorg)]
            maxOrg=indxOrganisms[maxiOrg[0]]
            firstmaxsite=list(sitesConflicts[maxOrg])[0]
            maxNamesite=firstmaxsite[0]
            maxNTsite=firstmaxsite
            valuess=NTsite_org_inSeq[maxNTsite][maxOrg].values()
            positions= [v for v in valuess ]
            for position in positions[0]:
                findedindx = [*range(position[0], position[1])]
                outputSeqNT[position[0]:position[1]] = list(maxNTsite)
                for finded in findedindx:
                    if finded in busyNTseq:
                        busyNTseq.remove(finded)
            if maxNamesite not in outputInsert:
                outputInsert[maxNamesite] = []
            outputInsert[maxNamesite] += [positions]
            if maxOrg not in orgInserted:
                orgInserted[maxOrg] = 0
            orgInserted[maxOrg] += 1
    outputSeqNT=''.join(outputSeqNT)
    return outputInsert,outputSeqNT


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