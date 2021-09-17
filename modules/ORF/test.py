import os
from Bio import SeqIO
from pathlib import Path

from Igem_TAU_2021.user_IO.input_main import parse_input
from Igem_TAU_2021.ORF.orf_main import orf_main

# input:
base_path = os.path.join(Path(__file__).parents[1], 'example_data')
user_inp_raw = {
    'sequence': os.path.join(base_path, 'mCherry_original.fasta'),
    'selected_promoters': None,
    'organisms': {
        'opt1': {'genome_path': os.path.join(base_path, 'Escherichia coli.gb'),
                           'optimized': True,
                           'expression_csv': None},
                 'deopt1': {'genome_path': os.path.join(base_path, 'Bacillus subtilis.gb'),
                            'optimized': False,
                            'expression_csv': None},
                 'deopt2': {'genome_path': os.path.join(base_path, 'Sulfolobus acidocaldarius.gb'),
                          'optimized': False,
                          'expression_csv': None},
                 'opt2': {'genome_path': os.path.join(base_path, 'Mycobacterium tuberculosis.gb'),
                             'optimized': True,
                             'expression_csv': None},

                 'opt3':{'genome_path': os.path.join(base_path, 'Pantoea ananatis.gb'),
                             'optimized': True,
                             'expression_csv': None},
                 'opt4': {'genome_path': os.path.join(base_path, 'Azospirillum brasilense.gb'),
                         'optimized': True,
                         'expression_csv': None}
                }
    }



user_inp = parse_input(user_inp_raw)

optimized_sequence = orf_main(user_inp)



print(optimized_sequence)
