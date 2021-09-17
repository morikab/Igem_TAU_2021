import os
import glob
import shutil
import sys
sys.path.append('..')
import xml.etree.ElementTree as et
from Bio import SeqIO
import re
from shared_functions_and_vars import write_fasta
import user_IO
from promoters.intersect_motifs_2_org_final import *
from promoters.globals_and_shared_methods import *
import pandas as pd
import numpy as np



#####todo:
#1 create pipeline from a to z                                                                               1-------- finished
#2 return best sequence - not write into file                                                                2-------- finished
#3 if directory already exists - overwrite it, shouldn't exit with an error                                  3-------- finished               
#4 copy mast.html into IO                                                                                    4-------- finished
#5 make sure all files are created in the right folder!!! promoters // promoters_not_for_user                5-------- finished
#6 mast should return file name for subsequent calls                                                         6-------- finished
#7 separate into several files: cmd // xml1 before mast // xml2 after mast. global vars in a separate        7-------- finished
#8 change extract_pssm_from_xml to work with element tree                                                    8-------- finished
#9 make sure folder hierarchy works!!!!!!!                                                                   9-------- finished

# anything else written in main

# things to check: is inter/100_200 neccesary for deopt organisms? make sure: RC not used altogether in MAST
# (currently there's a flag to not allow RC) -- in this case can delete check
# modify only best promoter????



"""
organism_dict = {
    'opt': ['Escherichia_coli', 'Mycobacterium_tuberculosis', 'Pantoea_ananatis', 'Azospirillum_brasilense'],
    'deopt': ['Bacillus_subtilis', 'Sulfolobus_acidocaldarius']
    }
"""









"""

base_path = 'C:\\wsl_files\\iGEM\\final\\example_data'
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

                    'opt3': {'genome_path': os.path.join(base_path, 'Pantoea ananatis.gb'),
                             'optimized': True,
                             'expression_csv': None},

                    'opt4': {'genome_path': os.path.join(base_path, 'Azospirillum brasilense.gb'),
                             'optimized': True,
                             'expression_csv': None}
    }
}

if '__name__' == '__main__':
    input_dict = user_IO.UserInputModule.run_module(user_inp_raw) #keys: sequence, selected_prom, organisms
    promoter_file_name = create_files_for_meme(input_dict)
    run_streme()
    motif_file_name = create_final_motif_xml(0.2, 0.2)
    mast_output_directory = run_mast(motif_file_name, promoter_file_name)
    return modify_promoter(promoter_file_name, mast_output_directory)

"""
