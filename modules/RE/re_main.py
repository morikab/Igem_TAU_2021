from modules.shared_functions_and_vars import translate
from .multi_org_functions import *
from .re_functions import insert_site_CDS

# initialize the logger object
logger = LoggerFactory.create_logger("Coding Sequence")


class REModule(object):
    @staticmethod
    def get_name() -> str:
        return "Restriction Enzymes"

    @staticmethod
    def run_module(input_dict, cds_nt):
        logger.info('###############################')
        logger.info('# CODING SEQUENCE ENGINEERING #')
        logger.info('###############################')

        logger.info('In this model, the ORF of the gene is analysed and synthetic changes are introduced into it in '
                    'order to:'
                    '\n     1. optimize translation and ribosome density for all optimized organisms, '
                    '\n        while deoptimizing it for the other group simultaneously'
                    '\n     2. remove restriction sites recognized by restriction enzymes from the optimized group'
                    '\n     3. insert restriction sites of enzymes present in the deoptimized group')

        cds_aa = translate(cds_nt)

        optimized_org_names, deoptimized_org_names = parse_inp1(input_dict)
        logger.info('Optimized organisms:')
        optimized_RE_dict = multi_organism_RE_dict(optimized_org_names, cds_aa)
        logger.info('Deoptimized organisms:')
        deoptimized_RE_dict = multi_organism_RE_dict(deoptimized_org_names, cds_aa)

        re_positions,  add_cds_nt= insert_site_CDS(deoptimized_RE_dict, cds_nt)
        logger.info(f'Positions of inserted sites {re_positions}')

        final_cds_nt = multi_org_remove_site(optimized_RE_dict, add_cds_nt)

        logger.info(f'Initial sequence before translation and restriction enzyme optimization: '
                    f'{cds_nt}\n')
        total_sequence_analysis(optimized_RE_dict, deoptimized_RE_dict, cds_nt)

        logger.info(f'Final sequence after translation and restriction enzyme optimization: '
                    f'{final_cds_nt}\n')
        total_sequence_analysis(optimized_RE_dict, deoptimized_RE_dict, final_cds_nt)
        print(add_cds_nt)

        return final_cds_nt #final nucleotide sequence after ORF and RE optimization
