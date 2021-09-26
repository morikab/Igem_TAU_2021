from modules.ORF.calculating_cai import general_geomean
import numpy as np
import pandas as pd


class ZscoreModule(object):

    @staticmethod
    def run_module(final_seq, inp_dict, optimization_type = 'cai'):
        '''
            This function is aimed to calculate an optimization index
            :param final_seq - the optimized sequence of the gene according to ORF and RE model (str)
            :param optimization_type - 'cai' or 'tai', default='cai'
            :param inp_dict: in the following format:
                {
            'sequence': the original sequence of the gene (str)
            'prom_list': None,
            'organisms': {scientific organism name1 : { 'tgcn': #tgcn dict {codon:number of occurences}
                                        '200bp_promoters': # prom_dict {gene name and function: prom}, promoter model
                                        'third_most_HE': # '400bp_promoters': prom400_dict,  # prom_dict {gene name and function: prom}, promoter model
                                        'gene_cds': # cds dict {gene name and function : cds}, for ORF model
                                        'intergenic': # intergenic dict {position along the genome: intergenic sequence}, promoter model
                                        'expression_estimation_af_all_genes': # when the expression csv is not given- the CAI is used as expression levels
                                        'CAI_score_of_all_genes': # {'gene_name': expression} ORF and promoter
                                        'cai_profile': # {'codon_name': cai_score} ORF model
                                        'optimized': # is the sequence in the optimized or deoptimized group- bool
                                       }
                         scientific organism name2 : ....

                }
            :return: all_Zscores - a data frame (for log file), rows: names of deoptimized organisms ,
                                                                columns: names of optimized organisms,
                                                                data: Z score ratio between each pair (optimized-deoptimized)
                     mean_Zscore - the mean Zscore ratio (float)
            '''

        scores = optimization_type + '_scores'
        weights = optimization_type + '_profile'

        opt_organisms = []
        deopt_organisms = []
        opt_Zscores_original = []
        opt_Zscores_eng = []
        deopt_Zscores_original = []
        deopt_Zscores_eng = []
        opt_index_org = []
        deopt_index_org = []
        for key, val in inp_dict['organisms'].items():
            miu = np.mean(np.array(list(val[scores].values())))
            sigma = np.std(np.array(list(val[scores].values())))
            index, _ = general_geomean([inp_dict['sequence'], final_seq], weights=val[weights])
            Zscores_original = (index[0] - miu) / sigma
            Zscores_eng = (index[1] - miu) / sigma
            index_org = (index[1] - index[0]) / sigma
            if val['optimized']:
                opt_organisms.append(key)
                opt_Zscores_original.append(Zscores_original)
                opt_Zscores_eng.append(Zscores_eng)
                opt_index_org.append(index_org)
            else:
                deopt_organisms.append(key)
                deopt_Zscores_original.append(1 / Zscores_original)
                deopt_Zscores_eng.append(1 / Zscores_eng)
                deopt_index_org.append(index_org)

        opt_Zscores_original = np.array([opt_Zscores_original])
        opt_Zscores_eng = np.array([opt_Zscores_eng])
        deopt_Zscores_original = np.array([deopt_Zscores_original]).T
        deopt_Zscores_eng = np.array([deopt_Zscores_eng]).T
        ratio_original = opt_Zscores_original * deopt_Zscores_original
        ratio_eng = opt_Zscores_eng * deopt_Zscores_eng
        Zscore_ratio = ratio_eng / ratio_original
        mean_Zscore = np.mean(Zscore_ratio)
        all_Zscores = pd.DataFrame(Zscore_ratio, index=deopt_organisms, columns=opt_organisms)
        mean_opt_index = np.mean(np.array(opt_index_org))
        mean_deopt_index = np.mean(np.array(deopt_index_org))
        norm_factor = max(mean_opt_index, mean_deopt_index)
        alfa = inp_dict['tuning_param']
        optimization_index = alfa * (mean_opt_index / norm_factor) - (1-alfa) * (mean_deopt_index / norm_factor)
        return mean_Zscore, all_Zscores, mean_opt_index, mean_deopt_index, optimization_index
