import os
from pathlib import Path
import shutil


global dna, start, opt_path, deopt_path, end, organism_dict
dna = "ACGT"
start = 'promoters_not_for_user'
o_path = 'opt_files'
d_path = 'deopt_files'
opt_path = os.path.join(start, o_path)
deopt_path = os.path.join(start, d_path)
end = '.fasta'
organism_dict = {'opt': [], 'deopt': []}


"""
Creates a new folder named dname.
If this folder already exists - overwrites it.

@param dname: name of a folder to create
"""
def create_folder(dname):
    if os.path.exists(dname):
        shutil.rmtree(dname)
    Path(dname).mkdir(parents=True)


