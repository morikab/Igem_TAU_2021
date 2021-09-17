import pandas as pd
from Igem_TAU_2021.user_IO.input_functions import find_tgcn

E_coli_web_tgcn_path = '/Users/Ilbres/studies/IGEM 2021/Model/Igem_TAU_2021/data/tRNA/ecoli_tGCN.txt'
Bac_web_tgcn_path = '/Users/Ilbres/studies/IGEM 2021/Model/Igem_TAU_2021/data/tRNA/bacillus_tGCN.txt'

# compute tgcn from gb files:
ecoli_gb_path = '/Users/Ilbres/studies/IGEM 2021/Model/Igem_TAU_2021_new/Igem_TAU_2021/example_data/Escherichia coli.gb'
bac_gb_path = '/Users/Ilbres/studies/IGEM 2021/Model/Igem_TAU_2021_new/Igem_TAU_2021/example_data/Bacillus subtilis.gb'

Ecoli_gb_dict = find_tgcn(ecoli_gb_path)
Bac_gb_dict = find_tgcn(bac_gb_path)

# turn everything into a dataframe
ecoli_web_df = pd.read_csv(E_coli_web_tgcn_path, delimiter='\t', index_col='AntiCodonsList')
bac_web_df = pd.read_csv(Bac_web_tgcn_path, delimiter='\t', index_col='AntiCodonsList')

ecoli_gb_df = pd.DataFrame.from_dict(Ecoli_gb_dict, orient='index')
bac_gb_df = pd.DataFrame.from_dict(Bac_gb_dict, orient='index')

ecoli_tgcn_df = pd.merge(ecoli_gb_df, ecoli_web_df)
bac_tgcn_df = pd.merge(bac_gb_df, bac_web_df)

ecoli_tgcn_df.to_excel('/Users/Ilbres/studies/IGEM 2021/Model/Igem_TAU_2021_new/Igem_TAU_2021/example_data/ecoli_tgcn_compared')
bac_tgcn_df.to_excel('/Users/Ilbres/studies/IGEM 2021/Model/Igem_TAU_2021_new/Igem_TAU_2021/example_data/bacillus_tgcn_compared')

