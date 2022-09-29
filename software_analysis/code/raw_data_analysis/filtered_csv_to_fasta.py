from modules.shared_functions_and_vars import write_fasta
import pandas as pd


org_dict = pd.read_csv('../../data/processed_genomes/filtered/cai_and_16s_for_genomes_filtered.csv', index_col=0).to_dict(orient='index')


names = []
seqs = []
for org, val in org_dict.items():
    seqid = org.split('.')[0]
    names.append(seqid)
    seq = val['16s']
    seqs.append(seq)

out_dir = '../../data/processed_genomes/filtered/'
write_fasta(out_dir + "cai_and_16s_for_genomes_filtered",
            list(seqs),
            list(names))

