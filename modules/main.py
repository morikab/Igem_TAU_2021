import os
import time
from pathlib import Path

from modules.logger_factory import LoggerFactory
from modules import Zscore_calculation, user_IO, RE, ORF

tic = time.time()
logger = LoggerFactory.create_logger("main")

base_path = os.path.join(os.path.dirname(__file__), '..\example_data')
user_inp_raw = {
    'sequence': os.path.join(base_path, 'mCherry_original.fasta'),
    'selected_promoters': None,
    'organisms': {
#                    'opt1': {'genome_path': os.path.join(base_path, 'Escherichia coli.gb'),
#                             'optimized': True,
#                             'expression_csv': os.path.join(base_path, 'ecoli_mrna_level.csv')},

                    # 'deopt1': {'genome_path': os.path.join(base_path, 'Bacillus subtilis.gb'),
                    #            'optimized': False,
                    #            'expression_csv': os.path.join(base_path, 'bacillus_mrna_level.csv')},
                    #
                    # 'deopt2': {'genome_path': os.path.join(base_path, 'Sulfolobus acidocaldarius.gb'),
                    #           'optimized': False,
                    #           'expression_csv': None},
                    #
                    # 'opt2': {'genome_path': os.path.join(base_path, 'Mycobacterium tuberculosis.gb'),
                    #          'optimized': True,
                    #          'expression_csv': None},

                    'opt3': {'genome_path': os.path.join(base_path, 'Pantoea ananatis.gb'),
                             'optimized': True,
                             'expression_csv': None},

                    'opt4': {'genome_path': os.path.join(base_path, 'Azospirillum brasilense.gb'),
                             'optimized': True,
                             'expression_csv': None}
    }
}


def run_modules(user_inp_raw):
    input_dict = user_IO.UserInputModule.run_module(user_inp_raw) #keys: sequence, selected_prom, organisms


    ### unit 1 ############################################
    orf_optimized_cds_nt_cai = ORF.ORFModule.run_module(input_dict, 'cai')
    # orf_optimized_cds_nt_tai = ORF.ORFModule.run_module(input_dict, 'tai')
    cds_nt_final_cai = RE.REModule.run_module(input_dict, orf_optimized_cds_nt_cai)  # todo: run both of them together to save time, or split creation of enzyme dict and the actual optimization (seems like a better solution)
    # cds_nt_final_tai = RE.REModule.run_module(input_dict, orf_optimized_cds_nt)
    mean_Zscore, all_Zscores = Zscore_calculation.ZscoreModule.run_module(cds_nt_final_cai, input_dict)#todo: update szcore to work on tai as well    print(mean_Zscore)
    #######################################################

    # ### unit 2 ############################################
    # best_seq, best_seq_name, evalue, synthetic_option = promoters.promoterModule.run_module(input_dict, D1 = 0.2, D2=0.2)
    # #######################################################

    # TODO - get zip_directory from the user
    zip_directory_path = os.path.join(str(Path(__file__).parent.resolve()), "artifacts")
    Path(zip_directory_path).mkdir(parents=True, exist_ok=True)
    final_output = user_IO.UserOutputModule.run_module(cds_sequence=cds_nt_final_cai,
                                                       zscore=mean_Zscore,
                                                       zip_directory=zip_directory_path)
    logger.info("Final output: %s", final_output)
    return final_output


if __name__ == "__main__":
    run_modules(user_inp_raw)
toc = time.time()

print('time: ', toc-tic)
