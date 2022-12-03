### moves all tls projects that are not 16s projects into a designated file
from os import listdir
from os.path import isfile, join
import subprocess

tls_blast_path = '../../data/genbank_tls/'

def run_cmd(cmd, verbose=False, *args, **kwargs):
    process = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        shell=True
    )
    std_out, std_err = process.communicate()
    if verbose:
        print(std_out.strip(), std_err)

tls_files = [f for f in listdir(tls_blast_path) if isfile(join(tls_blast_path, f))]
tls_csv_file_names = [f[:8] for f in tls_files if '.csv' if f]
non_needed_files = [f for f in tls_files if f[:8] not in tls_files]

run_cmd('mkdir non_16s_tls_data')
for f in non_needed_files:
    run_cmd('mv '+ f +' non_16s_tls_data')


