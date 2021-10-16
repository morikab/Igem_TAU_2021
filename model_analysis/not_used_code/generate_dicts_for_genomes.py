from  modules.user_IO.input_functions import *
from os import listdir
from os.path import isfile, join
import time
from not_used_code.running_modules_functions import *
from Bio.SeqIO import read


def sample_run(gb_files):  # do not use!!
    optimized_org_dict = {}
    for gb_path in gb_files[0:2]:
        org_name, org_dict = _parse_single_input(gb_path)
        optimized_org_dict[org_name] = org_dict

    deoptimized_org_dict = {}
    for gb_path in gb_files[3:5]:
        org_name, org_dict = _parse_single_input(gb_path)
        deoptimized_org_dict[org_name] = org_dict
    print('\n', final_run(cds, optimized_org_dict, deoptimized_org_dict))


# for all files in genomes directory- apply function and save the output

def _parse_single_input(gb_path):
    exp_csv_fid = None
    gb_file = SeqIO.read(gb_path, format='gb')
    org_name = find_org_name(gb_file)
    print(f'\n{org_name}')
    prom200_dict, cds_dict, intergenic_dict, estimated_expression = extract_gene_data(gb_path, exp_csv_fid)
    cai_weights = calculate_cai_weights_for_input(cds_dict, estimated_expression, exp_csv_fid)
    cai_scores = general_geomean(sequence_lst=cds_dict.values(), weights=cai_weights)
    miu = np.mean(np.array(cai_scores))
    sigma = np.std(np.array(cai_scores))

    org_dict = {
        'cai_profile': cai_weights,  # {dna_codon:cai_score}
        'avg': miu,  # float, avg cai score
        'std' : sigma
        }
    return org_name, org_dict


def input_json_for_all_org(gb_files, json_file_name):
    tic = time.time()
    total_dict  = {}
    for gb_path in gb_files:
        head, tail = os.path.split(gb_path)
        inner_tic = time.time()
        org_name, org_dict = _parse_single_input(gb_path)
        total_dict[tail[:-3]] = org_dict

        #sanity check- make sure the genome name matches the gb org name
        if org_name != tail[:-3]:
            print(tail[:-3])
        #############

        print(f'takes {time.time()-inner_tic}')
    toc = time.time()
    print(f'create dict {toc-tic}')
    # with open(json_file_name, 'w') as fp:
    #     json.dump(total_dict, fp)
    tic = time.time()
    print(f'write json {tic-toc}')



gene = read('../example_data/mcherry_original.fasta', 'fasta')
cds = str(gene.seq)
base_path = join(os.path.dirname(__file__), 'genomes')
gb_files = [join('genomes', f) for f in listdir(base_path) if isfile(join(base_path, f))]
input_json_for_all_org(gb_files, json_file_name ='../data_for_analysis/org_name_to_dict.json')




#### creating 16S dict
# df_16s = pd.read_excel('Supplementary_Data_1.xlsx', sheet_name='org_name_to_16s')
# org_name_to_16s = dict(zip(df_16s['org_name'].to_list(), df_16s['16s'].to_list()))
#
# with open('org_name_to_16s.json', 'w') as fp:
#     json.dump(org_name_to_16s, fp)
