from modules.logger_factory import LoggerFactory
from modules.promoters.before_meme_run import *
from modules.promoters.run_meme import *
from modules.promoters.motif_filtering_after_streme import *
from modules.promoters.promoter_choosing_after_mast import *

#todo:
# run_mast should return the f_name as an output
# as to all IO: write everything into a directory called "promoters_not_for_user", except the mast html file which
# should be copied into the "logs" directory
# divide code into different sub-modules
# make run_module return the best sequence and it's e value.
# re-rank promoters

logger = LoggerFactory.create_logger("promoters")


class promoterModule(object):

    @staticmethod
    def run_module(full_input_dict):
        logger.info('##########################')
        logger.info('# PROMOTERS INFORMATION #')
        logger.info('##########################')

        promoter_file_path = create_files_for_meme(full_input_dict)
        logger.info("promoters file path: %s", promoter_file_path)
        run_streme()
        tuning_param = full_input_dict['tuning_param']
        #D1 is optimization, D2 is selectivity
        # todo: make sure tamir agrees with this....
        D1 = 0.05 + tuning_param*0.3
        D2 = 0.05 + (1-tuning_param)*0.3
        motif_file_path = create_final_motif_xml(D1, D2)
        logger.info("motif file path: %s", motif_file_path)
        mast_output_folder = run_mast(motif_file_path, promoter_file_path)
        logger.info("mast output folder: %s", mast_output_folder)
        return modify_promoter(promoter_file_path, mast_output_folder)


