from main import run_modules
from shared_functions_and_vars import DEFAULT_ORGANISM_PRIORITY, write_fasta
from pathlib import Path
import os
import pandas as pd
from user_IO.user_input import extract_gene_data, calculate_cai_weights_for_input

#TODO: consider swapping cai calculation to yuri's algorithm

ecoli_genome_fid = '../example_data/data_for_experiment/Escherichia coli.gb'
ecoli_exp = '../example_data/data_for_experiment/ecoli_mrna_level.csv'
ecoli_cds_dict, estimated_expression = extract_gene_data(ecoli_genome_fid, ecoli_exp)
ecoli_cai_weights = calculate_cai_weights_for_input(ecoli_cds_dict, estimated_expression, None)

bacillus_genome_fid = '../example_data/data_for_experiment/Bacillus subtilis.gb'
bacillus_exp = '../example_data/data_for_experiment/ecoli_mrna_level.csv'
bacillus_cds_dict, estimated_expression = extract_gene_data(bacillus_genome_fid, bacillus_exp)
bacillus_cai_weights = calculate_cai_weights_for_input(bacillus_cds_dict, estimated_expression, None)



cai_weights = {}
cai_weights['ecoli'] = ecoli_cai_weights
cai_weights['bacillus'] = bacillus_cai_weights
v = pd.DataFrame(cai_weights)
print(v)
v.to_csv('../example_data/ecoli and bacillus cai.csv')

organisms_runs_to_test = {
    # 'B new CAI, wo exp': {
    #         'B. subtilis':{'genome_path': bacillus_genome_fid,
    #                 'optimized': True,
    #                 'expression_csv': None,
    #                 'optimization_priority': DEFAULT_ORGANISM_PRIORITY},
    #     'E. coli': {'genome_path': ecoli_genome_fid,
    #                     'optimized': False,
    #                     'expression_csv': None,
    #                     'optimization_priority': DEFAULT_ORGANISM_PRIORITY}
    # },
    #
    # 'E new CAI, wo exp': {
    #     'B. subtilis': {'genome_path': bacillus_genome_fid,
    #                 'optimized': False,
    #                 'expression_csv': None,
    #                 'optimization_priority': DEFAULT_ORGANISM_PRIORITY},
    #     'E. coli': {'genome_path': ecoli_genome_fid,
    #                     'optimized': True,
    #                     'expression_csv': None,
    #                     'optimization_priority': DEFAULT_ORGANISM_PRIORITY}
    # },

    'B new CAI, w exp': {
        'B. subtilis': {'genome_path': bacillus_genome_fid,
                    'optimized': True,
                    'expression_csv': bacillus_exp,
                    'optimization_priority': DEFAULT_ORGANISM_PRIORITY},
        'E. coli': {'genome_path': ecoli_genome_fid,
                        'optimized': False,
                        'expression_csv': ecoli_exp,
                        'optimization_priority': DEFAULT_ORGANISM_PRIORITY}
    },

    'E new CAI, w exp': {
        'B. subtilis': {'genome_path': bacillus_genome_fid,
                    'optimized': False,
                    'expression_csv': bacillus_exp,
                    'optimization_priority': DEFAULT_ORGANISM_PRIORITY},
        'E. coli': {'genome_path': ecoli_genome_fid,
                        'optimized': True,
                        'expression_csv': ecoli_exp,
                        'optimization_priority': DEFAULT_ORGANISM_PRIORITY}
    }

}



if __name__ == "__main__":

    pd.DataFrame()
    current_directory = Path(__file__).parent.resolve()
    base_path = os.path.join(Path(current_directory).parent.resolve(), "example_data")
    genome_path = os.path.join(base_path, 'arabidopsis_microbiome')

    final_output = {}
    tested_dict= {
        'sequence': os.path.join(base_path, 'mCherry_original.fasta'),
        'tuning_param': 0.5,
        'organisms': {},
        'clusters_count': 1,
    }

    for run_name, run_dict in organisms_runs_to_test.items():
        print(f'***********{run_name}***********')
        tested_dict['organisms'] = run_dict
        run_output = run_modules(user_input_dict = tested_dict)
        final_output[run_name] = run_output
    #print(final_output)

    csv_data = pd.DataFrame(final_output).transpose()
    fasta_dict = {key: value['final_sequence'] for key, value in final_output.items()}


    # csv_data.to_csv('../example_data/sequences_for_exp.csv')
    # write_fasta('../example_data/sequences_for_exp.fasta',
    #             list(fasta_dict.values()),
    #             list(fasta_dict.keys()))
    # print(fasta_dict)
