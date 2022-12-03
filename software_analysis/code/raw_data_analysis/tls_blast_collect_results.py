import pandas as pd

from os import listdir
from os.path import isfile, join

genomes_path = '../../data/processed_genomes/filtered/cai_and_16s_for_genomes_filtered.csv'
tls_blast_path = '../../data/genbank_tls/'
genomes_df = pd.read_csv(genomes_path)

tls_csv_file_names = [f for f in listdir(tls_blast_path) if isfile(join(tls_blast_path, f))]
print(tls_csv_file_names)


