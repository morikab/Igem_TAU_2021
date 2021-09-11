import numpy as np
from scipy.stats import stats
import re
import pandas as pd
import os

#TODO: make multi organism version

def extract_pssm_from_xml(f_name):
    """
    Extracts pssms of motifs from STREME output
    @param f_name: output of STREME as an XML file
    @return: dictionary where keys are motif sequences and items are corresponding pssms
    """
    try:
        with open(f_name, 'r') as f:
            data = f.read()
        pattern = '<motifs>.*?</motifs>'
        motifs = re.search(pattern, data, re.DOTALL)
        f.close()
        pssms = dict()
        sep_pat = '(?<=<motif ).*?</motif>'
        name_pat = '(?<=id=\").*?(?=\" alt)'
        width_pat = '(?<=width=\").*?(?=\" initial)'
        pos_pat = '(?<=<pos ).*?(?=/>)'
        freq_pat = '(?<=[ACGT]=\")[0-9.e-]*?(?=\")'

        if motifs != None:
            motifs = motifs.group(0)
            sep_motifs = re.findall(sep_pat, motifs, re.DOTALL)
            for m in sep_motifs:
                full_name = re.search(name_pat, m, re.DOTALL).group(0)
                n = full_name.index('-')
                index = int(full_name[:n])
                id_num = full_name[n + 1:]
                width = int(re.search(width_pat, m, re.DOTALL).group(0))
                all_pos = re.findall(pos_pat, m, re.DOTALL)
                df = pd.DataFrame(index=['A', 'C', 'G', 'T'])
                for i, pos in enumerate(all_pos):
                    freqs = re.findall(freq_pat, pos)
                    df[i + 1] = np.array(freqs, dtype=float)
                pssms[id_num] = df

        return pssms

    except IOError as e:
        print(e)


def padding_opt(v1, v2):
    """
    inserts uniform distributions at the edges for and calculates correlation between the two flattened pssms
    :param v1: the larger flattened vector (as a list)
    :param v2: the shorter flattened vector (as a list)
    :return:the highest correlation and pval between the motifs
    """
    pos_len_dif = int((len(v1) - len(v2)) / 4)
    corr = -1
    pval = -1
    for i in range(pos_len_dif + 1):
        padded_v2 = i * 4 * [0.25] + v2 + (pos_len_dif - i) * 4 * [0.25]
        current_corr, current_pval = stats.spearmanr(v1, padded_v2)
        if current_corr > corr:
            corr = current_corr
            pval = current_pval
    return corr, pval


def compare_pssms(pssm1, pssm2):
    """
    calculates corelation between 2 pssms
    :param pssm1: pssm for first motif as df
    :param pssm2: pssm for second motif as df
    :return: corr and p-value using spearman correlation
    """
    pssm1_vec = list(pssm1.to_numpy().flatten())
    pssm2_vec = list(pssm2.to_numpy().flatten())
    if len(pssm2_vec) == len(pssm1_vec):
        corr, pval = stats.spearmanr(pssm1_vec, pssm2_vec)

    elif len(pssm2_vec) > len(pssm1_vec):
        corr, pval = padding_opt(pssm2_vec, pssm1_vec)
    else:
        corr, pval = padding_opt(pssm1_vec, pssm2_vec)
    return corr, pval


def compare_pssm_sets(pssm1_dict, pssm2_dict):
    """
    make a df of correlations between all motif pairs
    :param pssm1_dict: formatted with the motif name as the key anf pssm (as a df) as value
    :param pssm2_dict: same format as the first set of motifs
    :return: a df with pssm1 as rows and pssm2 as columns
    """
    corr_df = pd.DataFrame()
    pval_df = pd.DataFrame()
    for motif1, pssm1 in pssm1_dict.items():
        for motif2, pssm2 in pssm2_dict.items():
            corr, pval = compare_pssms(pssm1, pssm2)
            corr_df.loc[motif1, motif2] = corr
            pval_df.loc[motif1, motif2] = pval
    return corr_df, pval_df


def find_selective_and_intergenic(selective_dict, intergenic_dict, final_percent_of_motifs=50):
    """
    find selective motifs to use for ranking the promoter options
    :param selective_dict: formatted with the motif name as the key anf pssm (as a df) as value, from the mast xml file
    of 50% highly of the optimised against 50% highly of deoptimized
    :param intergenic_dict: same format, from mast xml of all promoters against intergenic sequences
    :param final_percent_of_motifs: the percent of motifs from the initial set of selective that will be returned (those with the highest correlation)
    :return: a dict of {selective motif name: max(corr with intergenic)} for only for the motifs that were selected according to
    their correlation value and final_percent_of_motifs
    """
    corr_df, pval_df = compare_pssm_sets(selective_dict, intergenic_dict)
    max_corr_dict = corr_df.max(axis=1).to_dict()

    final_percent_of_motifs = round(len(max_corr_dict)*final_percent_of_motifs/100)
    corr_vals = list(max_corr_dict.values())
    corr_vals.sort(reverse=True)
    th = corr_vals[final_percent_of_motifs-1]
    selected_motifs = {motif:corr for motif, corr in max_corr_dict.items() if corr>=th}
    return selected_motifs



