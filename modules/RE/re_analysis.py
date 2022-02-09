from RE.multi_org_functions import *
from RE.re_functions import insert_site_CDS, translate
import json
from Bio.SeqIO import read
from modules.testing_for_modules import *
from shared_functions_and_vars import *
import pandas as pd
base_path = os.path.join(os.path.dirname(__file__))

f = open(os.path.join(base_path, 'REbase.txt'))
content = f.readlines()
REbase = [x.strip() for x in content]


def REbase_org(org):  # DONE
    '''
    parsing the REbase and creating a dictionary with the relevant sites from
    enzymes originated in org
    :param org: scientific name of investigated organism
    :param cds_aa: aa sequence of the cds
    :return: a dict in the following format:
                enzyme_dict[enzyme_name] = {
                'ambiguous_site':ambiguous_site, #contaiinng characters from the ambiguous_code dictionary
                'nt_to_aa': nt_to_aa} #non ambigiosus sites translated (in all 3 RFs)
    '''
    org = org.split(' ')
    org = org[0]# todo: this uses only the fam name and not specie
    print(org)
    RE_org_ind = [i for i, k in enumerate(REbase)
                  if org in k]  ## index of the 3rd section of the "block"
    print(RE_org_ind)
    enzyme_dict = {}  ## REbase dict.
    for k in RE_org_ind:
        org_list = []
        site = REbase[k + 2][3:]
        if '?' in site or ',' in site:
            continue
        enzyme_name = REbase[k - 2][3:]
        org_list.append(enzyme_name)

    return enzyme_dict


def single_run(org_list, n_org, percent_opt, cds_nt):
    opt_org_names = org_list

    opt_genomes = random.sample(opt_org_names,
                                round(n_org*percent_opt))

    deopt_genomes = random.sample([i for i in opt_org_names if i not in opt_genomes],
                                   round(n_org*(1-percent_opt)))

    cds_aa = translate(cds_nt)
    optimized_RE_dict = multi_organism_RE_dict(opt_genomes, cds_aa)
    deoptimized_RE_dict = multi_organism_RE_dict(deopt_genomes, cds_aa)
    re_positions,  add_cds_nt= insert_site_CDS(deoptimized_RE_dict, cds_nt)
    final_cds_nt = multi_org_remove_site(optimized_RE_dict, add_cds_nt)


    opt_sites = multi_org_final_found_sites(optimized_RE_dict, cds_nt)
    deopt_sites = multi_org_final_found_sites(deoptimized_RE_dict, final_cds_nt)

    return sum(list(opt_sites.values())), sum(list(deopt_sites.values()))






with open('C:\\Users\\labuser\\Documents\\Liyam\\tamir projects\\igem article\\Igem_TAU_2021\\model_analysis\\data_for_analysis\\org_name_to_dict.json') as f:
    org_dict = json.load(f)
gene = read('C:\\Users\\labuser\\Documents\\Liyam\\tamir projects\\igem article\\Igem_TAU_2021\\model_analysis\\zorA anti-phage defense.fasta', 'fasta')
cds_nt = str(gene.seq)

n_unknown = 0
RE_org_ind = [k[3:] for i, k in enumerate(REbase) if '<3>' in k]
for index, org  in enumerate(RE_org_ind):
    if 'Unidentified bacterium' in org:
        n_unknown +=1
    else:
        org = ' '.join(str.split(org, ' ')[:2])

RE_org_list = unique(RE_org_ind)
print('total number of organisms: ', len(RE_org_list))
print('total number of organisms unknown organisms: ', len([i for i in RE_org_list if 'Unidentified bacterium' in i]))


avg_runs = {}
n_runs = 10
for n_org in range(10,1000, 10):
    opt_list  = []
    deopt_list = []
    for i in range(n_runs):
        n_sites_from_opt, n_sites_from_deopt = single_run(org_list= RE_org_list, n_org=n_org, percent_opt=0.5, cds_nt=cds_nt)
        opt_list.append(n_sites_from_opt)
        deopt_list.append(n_sites_from_deopt)

    avg_runs[n_org] = [sum(opt_list)/n_runs, sum(deopt_list)/n_runs]
    print(n_org, [sum(opt_list)/n_runs, sum(deopt_list)/n_runs])

runs_df = pd.DataFrame.from_dict(avg_runs, orient='index', columns=['sites from wanted', 'sites from unwanted'])
print(runs_df)
runs_df.to_csv('RE_analysis_results.csv', index_label='')






#
# ##### tests from the arabidopsis microbiome ####################################33
# #### initial test
# inp_dict = generate_testing_data(n_organisms=32, percent_optimized = 0.5, genome_path = 'C:\\Users\\labuser\\Documents\\Liyam\\tamir projects\\igem article\\Igem_TAU_2021\\model_analysis\\data_for_analysis\\genomes')
# all_org_names = list(org_dict.keys())
# for org in all_org_names:
#     REbase_org(org)
#
# ## relevant for the ZorA ORF
# optimized_RE_dict = multi_organism_RE_dict(all_org_names, cds_aa)
# print('total number of relevant sites: ', [len(i) for i in list(optimized_RE_dict.values()) if i])
# print('number of organisms with sites: ', sum([1 for i in list(optimized_RE_dict.values()) if i]))
# ####todo:  result; from the zorA gene, only one site is present and is from the Rhodanobacter denitrificans species
# #
# # print(optimized_RE_dict)
#
# # #### second test
# # for n_org in range(5, 30,1):
# #     print(n_org)
# #     for i in range(10):
# #         inp_dict = generate_testing_data(n_organisms=n_org, percent_optimized = 0.5, genome_path = 'C:\\Users\\labuser\\Documents\\Liyam\\tamir projects\\igem article\\Igem_TAU_2021\\model_analysis\\data_for_analysis\\genomes')
# #         org_dict = inp_dict['organisms']
# #         cds_nt = str(gene.seq)
# #         print(org_dict)
# #         optimized_org_names = list(org_dict.keys())[:5]
# #
# #         deoptimized_org_names = list(org_dict.keys())[6:13]
# #
# #         cds_aa = translate(cds_nt)
# #
# #         optimized_RE_dict = multi_organism_RE_dict(optimized_org_names, cds_aa)
# #         deoptimized_RE_dict = multi_organism_RE_dict(deoptimized_org_names, cds_aa)
# #         print(optimized_RE_dict)
# #         re_positions,  add_cds_nt= insert_site_CDS(deoptimized_RE_dict, cds_nt)
# #         final_cds_nt = multi_org_remove_site(optimized_RE_dict, add_cds_nt)
# #
# #
# #         opt_sites = multi_org_final_found_sites(optimized_RE_dict, cds_nt)
# #         deopt_sites = multi_org_final_found_sites(deoptimized_RE_dict, final_cds_nt)
# #         print('opt, deopt sites: ', sum(list(opt_sites.values())), sum(list(deopt_sites.values())))
# #
