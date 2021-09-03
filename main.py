from logger_factory import LoggerFactory
import RE
import Zscore_calculation
import user_IO
import os

logger = LoggerFactory.create_logger("main")


base_path = os.path.join(os.path.dirname(__file__), 'example_data')
user_inp_raw = {
    'sequence': os.path.join(base_path, 'mCherry_original.fasta'),
    'selected_promoters': None,
    'organisms': {
#                    'opt1': {'genome_path': os.path.join(base_path, 'Escherichia coli.gb'),
#                             'optimized': True,
#                             'expression_csv': None},

                    'deopt1': {'genome_path': os.path.join(base_path, 'Bacillus subtilis.gb'),
                               'optimized': False,
                               'expression_csv': None},

                    # 'deopt2': {'genome_path': os.path.join(base_path, 'Sulfolobus acidocaldarius.gb'),
                    #           'optimized': False,
                    #           'expression_csv': None},
                    #
                    # 'opt2': {'genome_path': os.path.join(base_path, 'Mycobacterium tuberculosis.gb'),
                    #          'optimized': True,
                    #          'expression_csv': None},

                    # 'opt3': {'genome_path': os.path.join(base_path, 'Pantoea ananatis.gb'),
                    #          'optimized': True,
                    #          'expression_csv': None},

                    # 'opt4': {'genome_path': os.path.join(base_path, 'Azospirillum brasilense.gb'),
                    #          'optimized': True,
                    #          'expression_csv': None}
    }
}


def run_modules(user_inp_raw):
    input_dict = user_IO.UserInputModule.run_module(user_inp_raw)
#    cds_nt_final = RE.REModule.run_module(input_dict)
#    mean_Zscore, all_Zscores = Zscore_calculation.ZscoreModule.run_module(cds_nt_final, input_dict)
    mean_Zscore, all_Zscores = Zscore_calculation.ZscoreModule.run_module(input_dict['sequence'], input_dict)
    print(mean_Zscore)
    print(all_Zscores)
#    final_output = user_IO.UserOutputModule.run_module(cds_nt_final)
#    logger.info("Final output: %s", final_output)


if __name__ == "__main__":
    run_modules(user_inp_raw)
