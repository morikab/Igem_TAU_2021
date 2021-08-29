from collections import Counter
from Bio.Seq import Seq
import numpy as np


from tqdm import tqdm

def geometric_mean(iterable):
    a = np.log(iterable)
    return np.exp(a.mean())

def antiCodon_L(tgcn_dict):
    """This function get tRNA list of microorganism and return countered tGCN (tRNA Gene Copy Number)"""

    antiCodon_list = []
    for row in tqdm(tgcn_dict.keys()):
        antiCodon = row[-7:-4]
        antiCodon_list += [antiCodon]
    tGCN = Counter(antiCodon_list)
    #print_counter(tGCN)
    return tGCN

def print_counter(c):
    for t in c:
        print(str(t) + '\t' + str(c[t]))
    return


def weight_cal(cds_dict, Sij, tGCN):
    """ get the sequence with Sij table and tGCN (from prev. func) and return list of codons' weight."""

    my_seq = Seq(''.join(list(cds_dict.values()))) # merging all coding sequences into one big sequence
                                                    # (because the rest of the code treats an input as one sequence)
    print(len(my_seq))
    W_dict = {}
    for i in tqdm((np.arange(0, len(my_seq)-1, 3))):  # where to start scanning len%3!=0 - waiting for liyam answer

        codon = my_seq[i:i+3]
        codon_pos3 = my_seq[i+2]
        W = 0
        for key in Sij.keys():

            if key[2] == codon_pos3:
                S = Sij[key]
                i_codon = codon.reverse_complement()
                anticodon_pos1 = key[0]
                anti_codon = anticodon_pos1 + i_codon[1:3]
                tGC_curr = tGCN.get(str(anti_codon))

                if not tGC_curr:
                    continue

                else:
                    W += (1-S)*tGC_curr

        W_dict[codon] = W
    factor = max(W_dict.values())
    for k in W_dict:
        W_dict[k] = W_dict[k]/factor

    return W_dict


def TAI_cal(W_dict):

    w_norm_np = np.array(list(W_dict.values()))
    #w_norm_np = w_np/np.max(w_np)

    W_positive = w_norm_np[w_norm_np > 0]

    TAI_positive = geometric_mean(W_positive)
    W_norm_switched = W_dict.copy()
    for k, key in enumerate(W_norm_switched.keys()):
        if w_norm_np[k] == 0:
            W_norm_switched[key] = TAI_positive

    return W_norm_switched

# ---------------------------------------------------

class TAI(object):

    def __init__(self, cds_dict, tgcn):
        """
        :param cds_path: cds file (fasta)
        """
        self.Sij = Sij

        W_dict = weight_cal(cds_dict, self.Sij, tgcn)
        TAI_weight = TAI_cal(W_dict)
        self.index = TAI_weight


# ---------------------------------------------------
seq_file = r'..\..\data\tRNA\K12.txt'

# antiCodon:Codon
# antiCodon:Codon
# Sij = {'A:T': 0,
#        'G:C': 0,
#        'T:A': 0,
#        'C:G': 0,
#        'G:T': 1,
#        'A:C': 0.254414541,
#        'A:A': 0.810328488,
#        'T:G': 1}

# Sij = {'A:T': 0,
#        'G:C': 0,
#        'T:A': 0,
#        'C:G': 0.41,
#        'G:T': 0.63,
#        'A:C': 0.9749,
#        'A:A': 0.68,
#        'T:G': 0.95}

# todo: is that from c serevisiae?
Sij = {'A:T': 0,
       'G:C': 0,
       'T:A': 0,
       'C:G': 0.41,
       'G:T': 0.28,
       'A:C': 0.9999,
       'A:A': 0.68,
       'T:G': 0.89}

# Sij = {'A:T': 0,
#        'G:C': 0,
#        'T:A': 0,
#        'C:G': 0.561,
#        'G:T': 0.28,
#        'A:C': 0.9999,
#        'A:A': 0.68,
#        'T:G': 0.89}

#Good S's!
"""
Sij = {'A:T': 0,
       'G:C': 0,
       'T:A': 0,
       'C:G': 0,
       'G:T': 0.41,
       'A:C': 0.63,
       'A:A': 0.9749,
       'T:G': 0.68,
       'C:A': 0.95}
"""
"""
# Check if replacing is legal! U T / I A
tRNA_fa = r'..\..\data\tRNA\eschColi_K_12_MG1655-tRNAs.fa'
# SeC (TCA) -> in the codon table STOP
#  ACT not exceeded
# ACA

TAI_W = get_TAI(seq_file, Sij, tRNA_fa)
print(TAI_W)
"""

####################################################
# Deprecated functions:


def get_TAI(cds, Sij, tGCN):
    """

    :param cds: coding sequence
    :param Sij:
    :param tGCN:
    :return:
    """
    W_dict = weight_cal(cds, Sij, tGCN)
    return W_dict
