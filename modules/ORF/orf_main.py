from logger_factory.logger_factory import LoggerFactory
from modules import models

from .single_codon_optimization_method import optimize_sequence as optimize_sequence_by_single_codon
from .zscore_optimization_method import optimize_sequence_by_zscore_single_aa
from .zscore_optimization_method import optimize_sequence_by_zscore_bulk_aa

logger = LoggerFactory.get_logger()


class ORFModule(object):
    @staticmethod
    def run_module(
            user_input: models.UserInput,
            optimization_cub_score: models.OptimizationCubIndex,
            optimization_method: models.OptimizationMethod,
    ) -> str:
        logger.info("##########################")
        logger.info("# ORF #")
        logger.info("##########################\n")
        logger.info(F"Optimization score used is: {optimization_method}")

        target_gene = user_input.sequence

        if optimization_method.is_zscore_single_aa_optimization:
            return optimize_sequence_by_zscore_single_aa(
                sequence=target_gene,
                user_input=user_input,
                optimization_method=optimization_method,
                optimization_cub_index=optimization_cub_score,
            )

        if optimization_method.is_zscore_bulk_aa_optimization:
            return optimize_sequence_by_zscore_bulk_aa(
                sequence=target_gene,
                user_input=user_input,
                optimization_method=optimization_method,
                optimization_cub_index=optimization_cub_score,
            )

        if optimization_method.is_single_codon_optimization:
            return optimize_sequence_by_single_codon(target_gene=target_gene,
                                                     organisms=user_input.organisms,
                                                     optimization_cub_index=optimization_cub_score,
                                                     optimization_method=optimization_method,
                                                     tuning_param=user_input.tuning_parameter)

        raise ValueError(F"optimization method {optimization_method} is invalid")
