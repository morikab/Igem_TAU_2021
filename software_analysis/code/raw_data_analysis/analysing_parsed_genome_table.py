import pandas as pd
import matplotlib.pyplot as plt
import math
from collections import Counter
import operator
import json
from modules.shared_functions_and_vars import write_fasta
import random

'''
general plots and analysis to understand the output of the "parse_genome_file_job.py"
trying to understand the uniqueness of 16s sequences, and the effect of the redundancy in the CUB scores. 
'''


df = pd.read_csv('../../data/processed_genomes/cai_and_16s_for_genomes.csv')

rrna_seq_list = df['16s'].to_list()
only_str_16s = []
for idx, i in enumerate(rrna_seq_list):
    try:
        len(i)
        only_str_16s.append(i)
    except:
        print(i)
        print(idx)



# plt.hist([math.log10(len(i)) for i in only_str_16s], bins= 1000)
# plt.ylim([0,10])
# plt.xlim([3,3.3])
# plt.xlim([0,3])
# plt.show()


print(len(only_str_16s), len(sorted(set(only_str_16s))))
rrna_counts_dict = dict(Counter(only_str_16s))
rrna_counts = list(rrna_counts_dict.values())
most_repeated_seq = max(rrna_counts_dict.items(), key=operator.itemgetter(1))[0]
print('number of unique sequences', rrna_counts.count(1))
print('most repetitive sequence:', len(most_repeated_seq), ' number of repeats: ', max(rrna_counts))
# plt.hist(rrna_counts, bins=1000)
# plt.xlim([1,50])
# plt.show()
#
# plt.scatter(rrna_counts, [len(i) for i in seq_freq_dict.keys()])
# plt.ylabel('length')
# plt.xlabel('repeats')
# plt.show()


def entries_with_specified_16s(seq, cai_dict):
    repeated_seq_dict = {}
    for org, org_dict in cai_dict.items():
        try:
            if org_dict["16s"] == seq:
                repeated_seq_dict[org] = org_dict
        except:
            continue
    return repeated_seq_dict


def calc_std_for_seq_dict(repeated_seq_dict):
    most_repeated_seq_analysis = pd.DataFrame(repeated_seq_dict)
    most_repeated_seq_analysis.drop(['16s', '23s', '5s'], inplace=True)
    most_repeated_seq_analysis = most_repeated_seq_analysis.transpose()
    std_for_repeating_seq = [most_repeated_seq_analysis[col].std() for col in most_repeated_seq_analysis.columns]
    avg_std = sum(std_for_repeating_seq) / len(std_for_repeating_seq)

    return avg_std

def calculate_control(filtered_cai_dict, ctrl_size):
    org_list = list(filtered_cai_dict.keys())
    selected_org_list = [random.choice(list(org_list)) for i in range(ctrl_size)]
    selected_cai_dict = {i:filtered_cai_dict[i] for i in selected_org_list}
    avg_std = calc_std_for_seq_dict(selected_cai_dict)
    return avg_std



with open('../../data/processed_genomes/cai_and_16s_for_genomes.json', 'r') as f:
    cai_dict = json.load(f)

with open('../../data/processed_genomes/filtered/cai_and_16s_for_genomes_filtered.json', 'r') as f:
    filtered_cai_dict = json.load(f)


dict_for_plots = {}
dict_for_plots['repeats'] = []
dict_for_plots['avg_std_for_repeats'] = []
dict_for_plots['control_repeats'] = []
dict_for_plots['avg_std_for_control'] = []


n_repeats = 10

counter_idx = 0
for seq, count in rrna_counts_dict.items():
    if count <2:
        continue
    repeated_seq_dict = entries_with_specified_16s(seq, cai_dict)
    try:
        print(counter_idx)
        counter_idx +=1
        avg_std = calc_std_for_seq_dict(repeated_seq_dict)
        dict_for_plots['avg_std_for_repeats'].append(avg_std)
        dict_for_plots['repeats'].append(count)
        for i in range(n_repeats):
            dict_for_plots['control_repeats'].append(count)
            avg_std = calculate_control(filtered_cai_dict, ctrl_size = count)
            dict_for_plots['avg_std_for_control'].append(avg_std)

    except:
        print('skipped')
        continue
### plot of the avg std between
repeating_seqs = plt.scatter(dict_for_plots['repeats'], dict_for_plots['avg_std_for_repeats'])
ctrl_repeats = plt.scatter(dict_for_plots['control_repeats'], dict_for_plots['avg_std_for_control'])

plt.legend((repeating_seqs, ctrl_repeats),
           ('repeating_seqs', 'ctrl_repeats'),
           scatterpoints=1,
           loc='upper right',
           ncol=3,
           fontsize=8)

plt.xlabel('repeats')
plt.ylabel('average std')
plt.show()


#### writing fasta of seqs
# rrna_seq_for_fasta = {}
# for org, org_dict in cai_dict.items():
#     try:
#         rrna_seq_for_fasta[org] = org_dict['16s']
#     except:
#         continue
# write_fasta('../../data/processed_genomes/cai_and_16s_for_genomes',
#             list(rrna_seq_for_fasta.values()) , list(rrna_seq_for_fasta.keys()))


