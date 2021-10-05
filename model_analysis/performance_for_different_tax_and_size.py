from tax_splitting_func import choose_best_split
import json
from Bio.SeqIO import read
from running_modules_functions import final_run
import matplotlib.pyplot as  plt
import time
import numpy as np

# open organism dict
with open('data_for_analysis/org_name_to_dict.json') as f:
    org_dict = json.load(f)

gene = read('luciferase.fasta', 'fasta')
initial_seq = str(gene.seq)
runtimes_for_tax = 100
tax_to_color = {'Phylum':'Red', 'Family':'blue', 'Genus':'green'}
for tax_con, color in tax_to_color.items():
    print(tax_con, color)
    np_runs_for_tax = np.empty(shape=(len(org_dict)-2, runtimes_for_tax))
    for i in range(runtimes_for_tax):
        tic = time.time()
        scores_for_x = []
        for x in range(2,len(org_dict)):
            optimized_dict, deoptimized_dict = choose_best_split(org_dict, tax_con, x, run_times=50)
            opt_score = final_run(initial_seq, optimized_dict, deoptimized_dict)
            scores_for_x.append(opt_score)
        np_runs_for_tax[:,i]=scores_for_x
    opt_scores = np.nanmean(np_runs_for_tax, axis=1)
    sub_microbiome_size = np.array(list(range(2,len(org_dict))))
    plt.plot(sub_microbiome_size, opt_scores, color=color)

plt.title('Scale-up analysis for single codon optimization of translation')
plt.xlabel('Size of sub-microbiome')
plt.ylabel('Optimization score')
plt.legend(["Phylum", "Family", "Genus"], loc ="upper right") #todo: fix label problem
plt.show()