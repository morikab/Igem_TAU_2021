from logger_factory import LoggerFactory
import typing

# initialize the logger object
logger = LoggerFactory.create_logger("user_output")


class UserOutputModule(object):
    @staticmethod
    def get_name() -> str:
        return "User Output"

    @classmethod
    def run_module(cls, module_input: typing.Optional[typing.Dict] = None) -> typing.Dict:
        logger.info('##########################')
        logger.info('# USER OUTPUT INFORMATION #')
        logger.info('##########################')

        # TODO - fix the dict according to spec + create zip out of the relevant logs + data
        return {
            'final sequence: ': module_input["cds"],  # str
            'Zscore': None,  # int
            'final promoter': None,  # str
            'promoter score': None,  # int
            'promoter fasta': None,  # fasta file
        }
