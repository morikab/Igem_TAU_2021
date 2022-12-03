### moves all tls projects that are not 16s projects into a designated file
from os import listdir
from os.path import isfile, join

tls_blast_path = '../../data/genbank_tls/'

tls_files = [f for f in listdir(tls_blast_path) if isfile(join(tls_blast_path, f))]
tls_csv_file_names = [f[:8] for f in tls_files if '.csv' if f]
non_needed_files = [f for f in tls_files if f[:8] not in tls_files]

print('mkdir non_16s_tls_data')
for f in non_needed_files:
    print(f'mv {f} non_16s_tls_data')


