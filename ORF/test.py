import os
from Bio import SeqIO
from pathlib import Path

from Igem_TAU_2021.user_input import parse_input
from Igem_TAU_2021.ORF.orf_main import orf_main

# input:
base_path = os.path.join(Path(__file__).parents[1], 'example_data')
user_inp1_raw = {
    'opt1': {'genome_path': os.path.join(base_path, 'Escherichia coli.gb'),
             'optimized': True},
    'deopt1': {'genome_path': os.path.join(base_path, 'Bacillus subtilis.gb'),
               'optimized': False}
}

user_inp1 = parse_input(user_inp1_raw)
user_inp2 = SeqIO.read(os.path.join(base_path, 'mCherry_original.fasta'), "fasta")


optimized_sequence = orf_main(user_inp2, user_inp1)



print(optimized_sequence)
