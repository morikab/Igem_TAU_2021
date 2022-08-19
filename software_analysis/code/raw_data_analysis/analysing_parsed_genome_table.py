import pandas as pd
import matplotlib.pyplot as plt
import math
from collections import Counter
import operator
import json
from modules.shared_functions_and_vars import write_fasta

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


with open('../../data/processed_genomes/cai_and_16s_for_genomes.json', 'r') as f:
    cai_dict = json.load(f)

rrna_seq_for_fasta = {}
for org, org_dict in cai_dict.items():
    try:
        rrna_seq_for_fasta[org] = org_dict['16s']
    except:
        continue
write_fasta('../../data/processed_genomes/cai_and_16s_for_genomes',
            list(rrna_seq_for_fasta.values()) , list(rrna_seq_for_fasta.keys()))

counts_for_plt = []
std_for_plt = []
for seq, count in rrna_counts_dict.items():
    repeated_seq_dict = {}
    if count <2:
        continue
    for org, org_dict in cai_dict.items():
        try:
            if org_dict["16s"] == seq:
                repeated_seq_dict[org] = org_dict
        except:
            continue
    try:
        most_repeated_seq_analysis = pd.DataFrame(repeated_seq_dict)
        most_repeated_seq_analysis.drop(['16s', '23s', '5s'], inplace= True)
        most_repeated_seq_analysis = most_repeated_seq_analysis.transpose()
        std_for_repeating_seq = [most_repeated_seq_analysis[col].std() for col in most_repeated_seq_analysis.columns]
        avg_std = sum(std_for_repeating_seq)/len(std_for_repeating_seq)
        std_for_plt.append(avg_std)
        counts_for_plt.append(count)
    except:
        continue
#### plot of the avg std between
# plt.scatter(counts_for_plt, std_for_plt)
# plt.xlabel('repeats')
# plt.ylabel('average std')
# plt.show()




