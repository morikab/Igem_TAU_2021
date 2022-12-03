### make a random input for the analysis of genomes, not microbiomes
import typing
import pandas as pd
import random
from modules.shared_functions_and_vars import DEFAULT_ORGANISM_PRIORITY, nt_to_aa
from modules.models import UserInput

class Organism(object):
    def __init__(self,
                 df_row:pd.DataFrame,
                 name:str,
                 is_optimized : bool,
                 optimization_priority: float):
        self.name = name
        self.cai_profile = {key:df_row[key] for key in nt_to_aa.keys()}
        self.cai_avg = df_row['avg']
        self.cai_std = df_row['std']
        self.rrna = df_row['16s']
        self.n_proteins = df_row['n_proteins']
        self.is_optimized = is_optimized
        self.optimization_priority = optimization_priority



def make_organism_list(df,
                       is_optimized :dict,
                       optimization_priority = DEFAULT_ORGANISM_PRIORITY): #todo: check priority according to moran

    org_list = []
    for index, row in df.iterrows():
        org = Organism(row,
                       name = index,
                       is_optimized = is_optimized[index],
                       optimization_priority = optimization_priority)
        org_list.append(org)
    return org_list



def assign_opt_deopt(name_list:list, percent_optimized:float):
    opt_genomes = random.sample(name_list,
                                round(len(name_list)*percent_optimized))
    opt_status_dict = {
        i: i in opt_genomes
        for i in name_list
    }
    return opt_status_dict


def generate_testing_data(orf_seq : str,
                          n_organisms=15,
                          percent_optimized=0.5,
                          clusters_count= 2,
                          tuning_param=0.5,
                          genomes_path='../../../data/processed_genomes/cai_and_16s_for_genomes.csv'
                          ): #todo: check priority according to moran

    df = pd.read_csv(genomes_path, index_col=0)
    df = df.sample(n=n_organisms)
    opt_dict = assign_opt_deopt(list(df.index), percent_optimized)
    organisms_list = make_organism_list(df, opt_dict, optimization_priority= DEFAULT_ORGANISM_PRIORITY)

    model_input = UserInput(organisms=organisms_list,
                     sequence=orf_seq,
                     tuning_parameter=tuning_param,
                     optimization_method=None, #todo: make sure this leads to the latest optimization
                     clusters_count=clusters_count,
                     zip_directory='../zip_dir')


    return model_input

genomes_path = '../../../data/processed_genomes/cai_and_16s_for_genomes.csv'
single_inp = generate_testing_data(genomes_path = genomes_path,
                                             n_organisms=15,
                                             percent_optimized=0.5,
                                             clusters_count= 2,
                                             tuning_param=0.5,
                                             orf_seq= 'GTGGTGAGTAGCGATACGGTACTGATC')