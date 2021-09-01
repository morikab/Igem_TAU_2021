from RE.multi_org_functions import *
import typing
from ORF.orf_main import orf_main

# initialize the logger object
logger = LoggerFactory.create_logger("RE")


class REModule(object):
    @staticmethod
    def get_name() -> str:
        return "Restriction Enzymes"

    @classmethod
    def run_module(cls, input: typing.Optional[typing.Dict] = None) -> typing.Dict:
        logger.info('###############################')
        logger.info('# CODING SEQUENCE ENGINEERING #')
        logger.info('###############################')

        logger.info('In this model, the ORF of the gene is analysed and synthetic changes are introduced into it in '
                    'order to:'
                    '     1. optimize translation and ribosome density for all optimized organisms, '
                    '        while deoptimizing it for the other group simultaneously'
                    '     2. remove restriction sites recognized by restriction enzymes from the optimized group'
                    '     3. insert restriction sites of enzymes present in the deoptimized group')
        # cds_nt = orf_main(user_inp)
        cds_nt = input['sequence']
        cds_aa = translate(cds_nt)

        # finding the optimixed and deoptimized RE dictionaries
        optimized_org_names, deoptimized_org_names = parse_inp1(input)
        logger.info('Optimized organisms:')
        optimized_RE_dict = multi_organism_RE_dict(optimized_org_names, cds_aa)
        logger.info('Deptimized organisms:')
        deoptimized_RE_dict = multi_organism_RE_dict(deoptimized_org_names, cds_aa)

        # inserting sites for all optimized organisms
        # todo: use multi_organism_RE_dict after completing and finish this- do not call it just cds_nt
        add_cds_nt = multi_org_insert_site(deoptimized_RE_dict, cds_nt)

        # removing sites for all deoptimized organism
        final_cds_nt = multi_org_remove_site(optimized_RE_dict, add_cds_nt)

        # checking result:

        logger.info(f'Initial sequence before translation and restriction enzyme optimization: '
                    f'{cds_nt}\n')
        total_sequence_analysis(optimized_RE_dict, deoptimized_RE_dict, cds_nt)

        logger.info(f'Final sequence after translation and restriction enzyme optimization: '
                    f'{final_cds_nt}\n')
        total_sequence_analysis(optimized_RE_dict, deoptimized_RE_dict, final_cds_nt)

        return {
            "cds": final_cds_nt,
            "cds_nt": cds_nt,
        }
