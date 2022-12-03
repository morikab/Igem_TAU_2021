import traceback
from modules import sequence_family
from modules.ORF import ORFModule
from make_genome_groups import generate_testing_data
from modules.stats.evaluation import ZscoreModule

def run_modules( genomes_path = '../../../data/processed_genomes/cai_and_16s_for_genomes.csv',
                 n_organisms=15,
                 percent_optimized=0.5,
                 clusters_count= 2,
                 tuning_param=0.5,
                 orf_seq= 'GTGGTGAGTAGCGATACGGTACTGATC'):

    try:
        single_inp = generate_testing_data(genomes_path=genomes_path,
                                           n_organisms=n_organisms,
                                           percent_optimized = percent_optimized,
                                           clusters_count= clusters_count,
                                           tuning_param= tuning_param,
                                           orf_seq= orf_seq)

        clustered_user_inputs = sequence_family.SequenceFamilyModule.run_module(single_inp)


        for input_cluster in clustered_user_inputs:
            # TODO - what do we want to display for each run? We should store the results differently

            cds_nt_final_cai = ORFModule.run_module(input_cluster, 'cai', input_cluster.optimization_method)

            cai_mean_opt_index, cai_mean_deopt_index, cai_optimization_index, cai_weakest_score = \
                ZscoreModule.run_module(cds_nt_final_cai, input_cluster, optimization_type='cai')
            print(cai_optimization_index)


    except Exception:
        exception_str = traceback.format_exc()

        print("Encountered unknown error when running modules. Error message: %s", exception_str)

    return


run_modules( genomes_path = '../../../data/processed_genomes/cai_and_16s_for_genomes.csv',
                 n_organisms=15,
                 percent_optimized=0.5,
                 clusters_count= 2,
                 tuning_param=0.5,
                 orf_seq= 'GTGGTGAGTAGCGATACGGTACTGATC')