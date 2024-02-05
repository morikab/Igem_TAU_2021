from logger_factory.logger_factory import LoggerFactory
from modules import models
from modules.configuration import Configuration
from modules.run_summary import RunSummary

from .weak_folding_optimization import optimize_by_weak_folding

config = Configuration.get_config()
logger = LoggerFactory.get_logger()


class InitiationModule(object):
    @staticmethod
    def run_module(
            module_input: models.ModuleInput,
            run_summary: RunSummary,
    ) -> str:
        logger.info("##########################")
        logger.info("# initiation #")
        logger.info("##########################\n")

        return optimize_by_weak_folding(
            sequence=module_input.sequence,
            codons_num=config["INITIATION"]["NUMBER_OF_CODONS_TO_OPTIMIZE"],
            run_summary=run_summary,
        )
