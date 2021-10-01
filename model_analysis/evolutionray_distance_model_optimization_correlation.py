from running_modules_functions import *
import json
import random
from Bio.SeqIO import read, parse
from Bio import pairwise2
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
import time
from modules.shared_functions_and_vars import *


# ribosomal_msa_dict = {}
# for seq in parse('16s_sequences_msa_clustalo.fasta', 'fasta'):
#     ribosomal_msa_dict[seq.id.replace('_', ' ')] = str(seq.seq)
# with open('data_for_analysis/org_name_to_16s_msa.json', 'w') as fp:
#     json.dump(ribosomal_msa_dict, fp)

def diff_letters(a,b):
    return sum ( a[i] != b[i] for i in range(len(a)) )


def intersect(lst1, lst2):
    return set(lst1).intersection(lst2)

with open('data_for_analysis/org_name_to_dict.json', 'r') as fp:
    org_dict = json.load(fp)



with open('data_for_analysis/org_name_to_16s_msa.json', 'r') as fp:
    ribosomal_msa_dict = json.load(fp)

with open('data_for_analysis/org_name_to_16s.json', 'r') as fp:
    ribosomal_dict = json.load(fp)

# gene = read('../example_data/mcherry_original.fasta', 'fasta')
gene  = read('gcd_phosphurous_uptake.fasta', 'fasta')
gene = read('luciferase.fasta', 'fasta')
cds = str(gene.seq)
print(len(cds)/3)
aln_scores = []
opt_scores = []
msa_scores= []


tic = time.time()
for org1 in list(org_dict.keys()):
    for org2 in list(org_dict.keys()):
        # if org1.split(' ')[0] == org2.split(' ')[0]:
        #     continue
        if org1 == org2:
            continue
        optimized_dict = {org1: org_dict[org1]}
        deoptimized_dict = {org2: org_dict[org2]}
        optimization_index = final_run(cds, optimized_dict, deoptimized_dict)
        alignment_score = pairwise2.align.globalxx(
            ribosomal_dict[org1], ribosomal_dict[org2], score_only=True)
        aln_scores.append(alignment_score)
        msa_score = diff_letters(ribosomal_msa_dict[org1], ribosomal_msa_dict[org2])
        msa_scores.append(msa_score)
        opt_scores.append(optimization_index)
        # print(org1, org2, msa_scores, optimization_index)
toc = time.time()
print(toc-tic)

colors = np.where(pd.DataFrame(opt_scores)<0, 'C0', 'C1')
print(spearmanr(aln_scores, opt_scores, nan_policy='omit'))
plt.scatter(aln_scores, opt_scores, s=0.1, c=colors.ravel())
plt.title('comparison using pairwise alignment')
plt.xlabel('alignment score')
plt.ylabel('optimization index')
plt.ylim(-20, 20)
plt.show()

colors = np.where(pd.DataFrame(opt_scores)<0, 'C0', 'C1')

print(spearmanr(msa_scores, opt_scores, nan_policy='omit'))
plt.scatter(msa_scores, opt_scores, s=0.1, c=colors.ravel())
# plt.title('comparison using MSA')

positive_msa = [msa_scores[i] for i in range(len(opt_scores)) if opt_scores[i]>0]
negative_msa = [msa_scores[i] for i in range(len(opt_scores)) if opt_scores[i]<0]

positive_opt = [i for i in opt_scores if i>0]
negative_opt = [i for i in opt_scores if i<0]

print(len(positive_msa))
print(len(negative_opt))
optimized = plt.scatter(positive_msa, positive_opt, edgecolors='green' ,s=0.7)
nonoptimized = plt.scatter(negative_msa, negative_opt, edgecolors='red' ,s=0.7)
plt.xlabel('#Different chars for couple after multiple sequence alignment')
plt.ylabel('Optimization score')
plt.legend((optimized, nonoptimized), ('Optimized pairs', 'Non-optimized pairs'))
# plt.ylim(-20, 20)
plt.show()

print(spearmanr(positive_opt+[0]*len(negative_opt), positive_msa+negative_msa))
