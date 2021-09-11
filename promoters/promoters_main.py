from promoters.cmd_commands import *

#todo:
# run_mast should return the f_name as an output
# as to all IO: write everything into a directory called "promoters_not_for_user", except the mast html file which
# should be copied into the "logs" directory
# divide code into different sub-modules
# make run_module return the best sequence and it's e value.
# re-rank promoters



class promoterModule(object):

    @staticmethod
    def run_module(full_input_dict, D1, D2 ):
        create_files_for_meme(full_input_dict)
        run_streme()
        create_final_motif_xml(D1, D2)
        run_mast('final_motif_set.xml', 'promoters_path.fasta')
        modify_promoter('promoters.fasta', '')

        # return best_seq, best_seq_name, evalue, synthetic_option

