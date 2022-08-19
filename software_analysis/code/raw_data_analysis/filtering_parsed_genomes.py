import pandas as pd
import matplotlib.pyplot as plt
import math
from collections import Counter
import operator
import json
from modules.shared_functions_and_vars import write_fasta

'''
filtering the data coming from the "parse_genome_files_fob.py" according to length (filtering out sqeuncing shorter than 100bp) 
and using a randomly selected genome from sequences that have more than one 16S 
'''

with open('../../data/processed_genomes/cai_and_16s_for_genomes.json', 'r') as f:
    cai_dict = json.load(f)


only_str_16s = []
for org_dict in cai_dict.values():
    try:
        if len(org_dict['16s']) >1000:
            print(len(org_dict['16s']))
            only_str_16s.append(org_dict['16s'])
    except:
        continue

rrna_counts_dict = dict(Counter(only_str_16s))
rrna_counts = list(rrna_counts_dict.values())

filtered_org_dict = {}
counts_for_plt = []
std_for_plt = []
for seq, count in rrna_counts_dict.items():
    for org, org_dict in cai_dict.items():
        try:
            if org_dict["16s"] == seq:
                filtered_org_dict[org] = org_dict
                break
        except:
            continue

print(len(filtered_org_dict))

out_dir = '../../data/processed_genomes/'
with open(out_dir + "cai_and_16s_for_genomes_filtered.json", 'w') as handle:
    json.dump(filtered_org_dict, handle)

csv_data = pd.DataFrame(filtered_org_dict).transpose()
csv_data.to_csv(out_dir + "cai_and_16s_for_genomes_filtered.csv")

fasta_dict = {key: value['16s'] for key, value in filtered_org_dict.items()}
write_fasta(out_dir + "cai_and_16s_for_genomes_filtered.fasta",
            list(fasta_dict.values()),
            list(fasta_dict.keys()))