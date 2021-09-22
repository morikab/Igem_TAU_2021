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



class promoterModule(object):

    @staticmethod
    def run_module(full_input_dict):
        promoter_file_path = create_files_for_meme(full_input_dict)
        run_streme()
        tuning_param = full_input_dict['tuning_param']
        #D1 is optimization, D2 is selectivity
        D1 = tuning_param*0.5
        D2 = (1-tuning_param)*0.5
        motif_file_path = create_final_motif_xml(D1, D2)
        mast_output_folder = run_mast(motif_file_path, promoter_file_path)
        return modify_promoter(promoter_file_path, mast_output_folder)


