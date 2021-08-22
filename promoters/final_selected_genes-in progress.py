import os
import pandas as pd


# תעיפי את הפרומוטורים שלא ב 30% העליונים מבחינת ביטוי
# תעיפי את כל הפרומוטורים  שהם לא ב 15% העליונים לפי   TH=0.5 ולפי 400.
#
# ממה שנשאר, תדרגי לפי TH=0.5 לפי 200. תקחי את ה 4 העליונים .. (תשלחי לי את רשימת הפונקציות שלהם .. נראה שזה הגיוני).
#


def percent_highest_from_df(df, col_name, percent, ascending=True):
    final_excepted_row = round(df.shape[0]*percent/100)
    df.sort_values(col_name, inplace = True, ascending=ascending)
    final_df = df.iloc[:final_excepted_row, :]
    return final_df


def th_from_path(fid):
    return fid.split('=')[-1].split(' ')[0]


def extract_all_mast_table_fids(prom_len, org, base_fid ):
    dir = os.path.join(base_fid, str(prom_len)+'bp promoters')
    list_of_paths = [os.path.join(dir, i) for i in os.listdir(dir) if 'table' in i and org in i]#+ os.listdir(os.path.join(base_fid, '400bp promoters')))
    return list_of_paths

def open_sort_and_rank(fid):
    mast_df = pd.read_csv(fid, delimiter='\s+')
    mast_df = mast_df.drop(['DESCRIPTION', 'E-VALUE', 'LENGTH'], axis=1)
    mast_df.columns = ['gene_name', 'evalue']
    mast_df = mast_df[1:-1]
    mast_df.head()
    mast_df['RANK'] = mast_df.index
    mast_df.sort_values('gene_name', inplace=True)
    mast_df.reset_index(drop=True, inplace=True)
    return mast_df

# org = 'bacillus'
org = 'ecoli'
exp_lvl_df =pd.read_csv('outputs\\'+ org+'\\cds_data_'+org+'.csv')
base_fid = 'outputs\\streme_motifs\\50%\\motifs by th\\only_sense_strand'
prom200 = extract_all_mast_table_fids(200, org, base_fid)
prom400 = extract_all_mast_table_fids(400, org, base_fid)

exp_lvl_df.drop(['cds_start', 'cds_stop', 'cds_strand', 'cds_seq'], axis = 1, inplace= True)
exp_lvl_df['exp rank'] = exp_lvl_df.index
exp_lvl_df.sort_values('gene_name', inplace=True)
exp_lvl_df.reset_index(drop=True, inplace=True)

res_df = exp_lvl_df

th_options = []
for fid in prom200:
    th = th_from_path(fid)
    th_options.append(th)
    col_name = 'th=' + th + ', len=200'
    mast_df = open_sort_and_rank(fid)
    res_df[col_name + ', evalue'] = mast_df.evalue
    res_df[col_name + ', rank'] = mast_df.RANK

for fid in prom400:
    col_name = 'th=' + th_from_path(fid) + ', len=400'
    mast_df = open_sort_and_rank(fid)
    res_df[col_name + ', evalue'] = mast_df.evalue
    res_df[col_name + ', rank'] = mast_df.RANK



# print(res_df.columns)
filtered_by_exp_level = percent_highest_from_df(res_df, col_name='exp rank', percent=30)
filtered_by_prom400 = percent_highest_from_df(res_df, col_name='th='+str(max(th_options))+', len=400, rank', percent=15)
prom400_and_exp_level = pd.merge(filtered_by_prom400, filtered_by_exp_level)


final_filtered = prom400_and_exp_level.sort_values('th=' + str(max(th_options)) + ', len=200, rank')
final_filtered.to_csv(os.path.join(base_fid, org + ' filtered_options.csv'), index=False)
