import typing

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
    ) -> typing.Tuple[str, int]:
        logger.info("##########################")
        logger.info("# Initiation #")
        logger.info("##########################\n")

        codons_num = config["INITIATION"]["NUMBER_OF_CODONS_TO_OPTIMIZE"]
        logger.info(f"Running initiation optimization on {codons_num} codons at the start of the ORF")
        optimized_sequence = optimize_by_weak_folding(
            sequence=module_input.sequence,
            codons_num=codons_num,
            run_summary=run_summary,
        )
        return optimized_sequence, codons_num
