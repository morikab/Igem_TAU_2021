from logger_factory import LoggerFactory
import typing

# initialize the logger object
logger = LoggerFactory.create_logger("user_output")


class UserOutputModule(object):
    @staticmethod
    def get_name() -> str:
        return "User Output"

    @staticmethod
    def run_module(cds_sequence, zscore):
        logger.info('###########################')
        logger.info('# USER OUTPUT INFORMATION #')
        logger.info('###########################')

        # TODO - fix the dict according to spec + create zip out of the relevant logs + data
        return {
            'final sequence: ': cds_sequence,  # str
            'Zscore': zscore,  # int
            'final promoter': None,  # str
            'promoter score': None,  # int
            'promoter fasta': None,  # fasta file
        }
