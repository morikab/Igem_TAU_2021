import json
import pandas as pd
from Bio import  SeqIO
import matplotlib.pyplot as plt
'''
analysing the csv outputs produced by the cross_tls_with_genome_blast_job.py
currently work on commands that create a csv with only the best hit. 
'''


genomes_df = pd.read_csv('../../data/tls_genome_match/tls_genome_matches.csv')

with open('../../data/tls_genome_match/tls_genome_matches.json', 'r') as f:
    tls_dict = json.loads(f.read())

used_genomes = []
for tls, data in tls_dict.items():
    used_genomes += list(data['cai'].keys())

used_genomes = [int(i) for i in (sorted(set(used_genomes)))]
used_genomes.sort()
print('sequences found in blast', used_genomes) #right now, all microbiomes use a collective set of 35

# ### non aligned
rrna_example_seqs = SeqIO.parse('../../data/processed_genomes/16s_fasta_example.fasta', 'fasta')
rrna_example_seqs = SeqIO.to_dict(rrna_example_seqs)
rrna_example_seqs = {int(idx):str(obj.seq) for idx, obj in rrna_example_seqs.items()}

only_seqs = list(rrna_example_seqs.values())

repeats_for_used = []
repeats_for_unused = []
for i in range(len(only_seqs)):
    if i not in used_genomes:
        repeats_for_unused.append(only_seqs.count(only_seqs[i]))
    else:
        repeats_for_used.append(only_seqs.count(only_seqs[i]))

print('repeats for unused: ', sum(repeats_for_unused)/len(repeats_for_unused))
print('repeats for used: ', sum(repeats_for_used)/len(repeats_for_used))

used_rrna = {idx:seq for idx, seq in rrna_example_seqs.items() if idx in used_genomes}
print('avg length of found 16s: ', sum([len(i) for i in used_rrna.values() ])/len(used_genomes))
unused_rrna = {idx:seq for idx, seq in rrna_example_seqs.items() if idx not in used_genomes}
print('avg length of not found 16s: ', sum([len(i) for i in unused_rrna.values() ])/len(unused_rrna))

plt.scatter([len(i) for i in unused_rrna.values() ], repeats_for_unused, s=100, alpha=0.5, linewidths=0)
plt.scatter([len(i) for i in used_rrna.values()], repeats_for_used, s=100, alpha=0.5, linewidths=0)
plt.xlabel('sequence length')
plt.ylabel('number of repeats')
plt.show()

# plt.hist([len(i) for i in used_rrna.values()], alpha= 0.5, bins=100, density=True, label='found seqs')
# plt.hist([len(i) for i in unused_rrna.values() ], alpha= 0.5, bins=100, density=True, label='non found seqs')
# plt.xlabel('16s length')
# plt.ylabel('normalized distribution')
# plt.legend()
# plt.show()




# ###aligned
# rrna_example_seqs_aligned = SeqIO.parse('../../data/processed_genomes/16s_fasta_example_aligned.fasta', 'fasta')
# rrna_example_seqs_aligned = SeqIO.to_dict(rrna_example_seqs_aligned)
# rrna_example_seqs_aligned = {int(idx):str(obj.seq).strip('-') for idx, obj in rrna_example_seqs_aligned.items()}
#
#
# used_rrna = {idx:seq for idx, seq in rrna_example_seqs_aligned.items() if idx in used_genomes}
# print('avg #gaps of found 16s: ', sum([i.count('-') for i in used_rrna.values() ])/len(used_genomes))
# unused_rrna = {idx:seq for idx, seq in rrna_example_seqs_aligned.items() if idx not in used_genomes}
# print('avg #gaps of not found 16s: ', sum([i.count('-') for i in unused_rrna.values() ])/len(unused_rrna))
#
#
# plt.hist([i.count('-') for i in used_rrna.values()], alpha= 0.5, bins=100, density=True, label='found seqs')
# plt.hist([i.count('-') for i in unused_rrna.values() ], alpha= 0.5, bins=100, density=True, label='non found seqs')
# plt.xlabel('#gaps in msa')
# plt.ylabel('normalized distribution')
# plt.legend()
# plt.show()
#
