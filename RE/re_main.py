from RE.multi_org_functions import *

# initialize the logger object
logger = LoggerFactory.create_logger("RE")


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
        logger.info('Deptimized organisms:')
        deoptimized_RE_dict = multi_organism_RE_dict(deoptimized_org_names, cds_aa)

        add_cds_nt = multi_org_insert_site(deoptimized_RE_dict, cds_nt)

        final_cds_nt = multi_org_remove_site(optimized_RE_dict, add_cds_nt)

        logger.info(f'Initial sequence before translation and restriction enzyme optimization: '
                    f'{cds_nt}\n')
        total_sequence_analysis(optimized_RE_dict, deoptimized_RE_dict, cds_nt)

        logger.info(f'Final sequence after translation and restriction enzyme optimization: '
                    f'{final_cds_nt}\n')
        total_sequence_analysis(optimized_RE_dict, deoptimized_RE_dict, final_cds_nt)

        return final_cds_nt #final nucleotide sequence after ORF and RE optimization
