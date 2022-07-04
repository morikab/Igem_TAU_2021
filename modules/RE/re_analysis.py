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


    n_opt_sites = sum(list(opt_sites.values()))
    n_opt_org = sum([1 for i in opt_sites.values() if i>0 ])
    n_deopt_sites = sum(list(deopt_sites.values()))
    n_deopt_org = sum([1  for i in deopt_sites.values() if i>0 ])
    return n_opt_sites, n_deopt_sites, n_opt_org, n_deopt_org






with open(r'C:\Users\97252\Documents\work\Igem_TAU_2021\model_analysis\data_for_analysis\org_name_to_dict.json') as f:
    org_dict = json.load(f)

with open(r'C:\Users\97252\Documents\work\Igem_TAU_2021\model_analysis\zorA anti-phage defense.fasta') as fasta:
    gene = read(fasta, 'fasta')
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

def different_sizes ():
    avg_runs = {}
    n_runs = 10
    for n_org in range(55, 100, 5):
        opt_list  = []
        deopt_list = []
        opt_org = []
        deopt_org = []
        for i in range(n_runs):
            n_opt_sites, n_deopt_sites, n_opt_org, n_deopt_org = single_run(org_list= RE_org_list, n_org=n_org, percent_opt=0.5, cds_nt=cds_nt)
            opt_list.append(n_opt_sites)
            deopt_list.append(n_deopt_sites)
            opt_org.append(n_opt_org)
            deopt_org.append(n_deopt_org)

        avg_runs[n_org] = [sum(i)/n_runs for i in [opt_list, deopt_list, opt_org, deopt_org]]
        print(n_org, avg_runs[n_org])

    runs_df = pd.DataFrame.from_dict(avg_runs, orient='index', columns=['sites from wanted', 'sites from unwanted', 'org from wanted', 'org from unwanted'])
    print(runs_df)
    runs_df.to_csv('different sizes.csv', index_label='')
    return None


def different_ratios ():
    avg_runs = {}
    n_runs = 10
    for percent_opt in [5, 90, 95]:
        percent_opt = percent_opt/100
        opt_list  = []
        deopt_list = []
        opt_org = []
        deopt_org = []
        for i in range(n_runs):
            n_opt_sites, n_deopt_sites, n_opt_org, n_deopt_org = single_run(org_list= RE_org_list, n_org=30, percent_opt=percent_opt, cds_nt=cds_nt)
            opt_list.append(n_opt_sites)
            deopt_list.append(n_deopt_sites)
            opt_org.append(n_opt_org)
            deopt_org.append(n_deopt_org)

        avg_runs[percent_opt] = [sum(i)/n_runs for i in [opt_list, deopt_list, opt_org, deopt_org]]
        print(percent_opt, avg_runs[percent_opt])

    runs_df = pd.DataFrame.from_dict(avg_runs, orient='index', columns=['sites from wanted', 'sites from unwanted', 'org from wanted', 'org from unwanted'])
    print(runs_df)
    runs_df.to_csv('analysis_results/different_ratios.csv', index_label='')

different_ratios ()

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
