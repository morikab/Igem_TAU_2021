### fixing the warning error described in parse genome files header
### returns the fixed df, and the filtered version + unique values
### filteres ot rows with 16s shorter than 1000bp and less than 450 proteins
import pandas as pd
from modules.shared_functions_and_vars import write_fasta


n_protein_th = 450
len_16s_th = 1000

df = pd.read_csv('../../data/processed_genomes/cai_and_16s_for_genomes_with_errors.csv', index_col= 0)
print(df)


for index, row in df.iterrows():
    rrna_lst = [row['16s'], row['23s'], row['5s']]
    try:
        rrna_lst.sort(key=len)
        df.at[index, '5s'] = rrna_lst[0]
        df.at[index, '16s'] = rrna_lst[1]
        n_prot = df.at[index, 'n_proteins']


        if len(rrna_lst[1])<len_16s_th:
            df.drop(index, axis=0, inplace=True)
            print(f'{index} dropped because len 16s gene = {len(rrna_lst[1])}, shorter than 1000bp')
        elif n_prot< n_protein_th:
            df.drop(index, axis=0, inplace=True)
            print(f'{index} dropped because n_protein = {n_prot}, smaller than 450')
        else:
            df.at[index, '23s'] = rrna_lst[2]

    except:
        print(f'{index} dropped because value type error')
        df.drop(index, axis=0, inplace=True)


df.to_csv('../../data/processed_genomes/cai_and_16s_for_genomes.csv', index_label= ' ')
print(df)

df_filtered = df.drop_duplicates(subset=['16s'], keep='first')
df_filtered.to_csv('../../data/processed_genomes/filtered/cai_and_16s_for_genomes_filtered.csv', index_label= ' ')
print(df_filtered)

write_fasta('../../data/processed_genomes/filtered/cai_and_16s_for_genomes_filtered',
            list(df_filtered['16s']),
            [i.split('.')[0] for i in list(df_filtered.index)])

# my findings are that flattening the dict mixes up the different RRNA types,
# I'm trying to fix this and reassign the correct values to each in this code (s5 is smallest, 16s middle and 23s largest
# this will also be used for filtering the 16s RRNA differences
# after this, I will need to run the metagenome matching thing again (the one that uses blast) after fixing the alignment and blastdb
