from modules.ORF.orf_main import *
from modules.Zscore_calculation import ZscoreModule
from modules.ORF.calculating_cai import *
from statistics import mean

import pandas as pd
import numpy as np



def run_orf_module(optimized_org_dict, deoptimized_org_dict, target_gene):
    '''
    target_gene: DNA for gene to optimize
    optimized_org_dict: {org_name:org_dict}, containing only optimized organisms
    optimized_org_dict: {org_name:org_dict}, containing only optimized organisms
    ** org dict format: {
        'cai_profile': cai_weights,  # {dna_codon:cai_score}
        'cai_scores': cai_scores,  # {'gene_name': score}
        }
    returns: new ORF
    '''
    high_expression_organisms = [
        Organism(name=org_name, tai_weights={}, cai_weights=org_dict['cai_profile'], feature_to_generate='cai',
                 cai_std=org_dict['std'], tai_std=None )
        for org_name, org_dict in optimized_org_dict.items()]

    low_expression_organisms = [
        Organism(name=org_name, tai_weights={}, cai_weights=org_dict['cai_profile'], feature_to_generate='cai',
                 cai_std=org_dict['std'], tai_std=None )
        for org_name, org_dict in deoptimized_org_dict.items()]

    optimized_sequence = optimize_sequence(target_gene=target_gene,
                                           high_expression_organisms=high_expression_organisms,
                                           low_expression_organisms=low_expression_organisms,
                                           tuning_param=0.5,
                                           local_maximum=False)
    return optimized_sequence

#
def run_zscore_module(final_seq, initial_seq, optimized_org_dict, deoptimized_org_dict, optimization_type = 'cai'):

    weights = optimization_type + '_profile'
    opt_index_org = []
    deopt_index_org = []
    for org_name, org_dict in {**optimized_org_dict, **deoptimized_org_dict}.items():
        sigma = org_dict['std']
        index = general_geomean([initial_seq, final_seq], weights=org_dict[weights])
        initial_score = index[0]
        final_score = index[1]
        index_org = (final_score - initial_score) / sigma
        if org_name in optimized_org_dict.keys():
            opt_index_org.append(index_org)
        else:
            deopt_index_org.append(index_org)

    mean_opt_index = mean(opt_index_org)
    mean_deopt_index = mean(deopt_index_org)
    # norm_factor = max(mean_opt_index, mean_deopt_index)
    alfa = 0.5
    optimization_index = (alfa * mean_opt_index - (1 - alfa) * (mean_deopt_index))  # / norm_factor
    return mean_opt_index, mean_deopt_index, optimization_index



# def run_zscore_module(final_seq, initial_seq, optimized_org_dict, deoptimized_org_dict, optimization_type = 'cai'):
#
#     scores = optimization_type + '_scores'
#     weights = optimization_type + '_profile'
#
#     opt_organisms = []
#     deopt_organisms = []
#     opt_Zscores_original = []
#     opt_Zscores_eng = []
#     deopt_Zscores_original = []
#     deopt_Zscores_eng = []
#     opt_index_org = []
#     deopt_index_org = []
#     for org_name, org_dict in {**optimized_org_dict, **deoptimized_org_dict}.items():
#         miu = org_dict['avg']
#         sigma = org_dict['std']
#         index= general_geomean([initial_seq, final_seq], weights=org_dict[weights])
#         Zscores_original = (index[0] - miu) / sigma
#         Zscores_eng = (index[1] - miu) / sigma
#         index_org = (index[1] - index[0]) / sigma
#         if org_name in optimized_org_dict.keys():
#             opt_organisms.append(org_name)
#             opt_Zscores_original.append(Zscores_original)
#             opt_Zscores_eng.append(Zscores_eng)
#             opt_index_org.append(index_org)
#         else:
#             deopt_organisms.append(org_name)
#             deopt_Zscores_original.append(1 / Zscores_original)
#             deopt_Zscores_eng.append(1 / Zscores_eng)
#             deopt_index_org.append(index_org)
#
#     opt_Zscores_original = np.array([opt_Zscores_original])
#     opt_Zscores_eng = np.array([opt_Zscores_eng])
#     deopt_Zscores_original = np.array([deopt_Zscores_original]).T
#     deopt_Zscores_eng = np.array([deopt_Zscores_eng]).T
#     ratio_original = opt_Zscores_original * deopt_Zscores_original
#     ratio_eng = opt_Zscores_eng * deopt_Zscores_eng
#     Zscore_ratio = ratio_eng / ratio_original
#     mean_Zscore = np.mean(Zscore_ratio)
#     all_Zscores = pd.DataFrame(Zscore_ratio, index=deopt_organisms, columns=opt_organisms)
#     mean_opt_index = np.mean(np.array(opt_index_org))
#     mean_deopt_index = np.mean(np.array(deopt_index_org))
#     norm_factor = max(mean_opt_index, mean_deopt_index)
#     alfa = 0.5
#     optimization_index = alfa * (mean_opt_index / norm_factor) - (1 - alfa) * (mean_deopt_index / norm_factor)
#     return mean_Zscore, all_Zscores, mean_opt_index, mean_deopt_index, optimization_index

def final_run(initial_seq, optimized_org_dict, deoptimized_org_dict):
    optimized_sequence = run_orf_module(optimized_org_dict, deoptimized_org_dict, initial_seq)
    mean_opt_index, mean_deopt_index, optimization_index =  \
        run_zscore_module(optimized_sequence, initial_seq, optimized_org_dict, deoptimized_org_dict)
    return optimization_index