import os
from Bio import SeqIO
from pathlib import Path

from Igem_TAU_2021.user_input import parse_input
from Igem_TAU_2021.ORF.orf_main import orf_main

# input:
base_path = os.path.join(Path(__file__).parents[1], 'example_data')
user_inp1_raw = {
    'sequence':SeqIO.read(os.path.join(base_path, 'mCherry_original.fasta'), "fasta"),
    'selected_promoters': os.path.join(base_path, 'ORF_optimized_sequences.fasta'), #or None, # or a fasta file of promoter name and promoter
    'opt1': {'genome_path': os.path.join(base_path, 'Escherichia coli.gb'),
             'expression_csv': None, #todo: add treatment for expression data
             'optimized': True},
    'deopt1': {'genome_path': os.path.join(base_path, 'Bacillus subtilis.gb'),
               'expression_csv': None,
               'optimized': False}
}
user_inp = parse_input(user_inp1_raw)


optimized_sequence = orf_main(user_inp)



print(optimized_sequence)
