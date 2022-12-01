from functools import partial

from logger_factory.logger_factory import LoggerFactory
from modules import models

from .loss_function_optimization_method import optimize_sequence_by_loss_function
from .zscore_optimization_method import hill_climbing_optimize_by_zscore as hill_climbing_optimize_by_zscore
from .zscore_optimization_method import hill_climbing_optimize_aa_bulk_by_zscore as \
    hill_climbing_optimize_aa_bulk_by_zscore

logger = LoggerFactory.get_logger()

# todo: add a statistical analysis of how close the organisms are - like what is the best codon for eah AA
#  and are they close


class ORFModule(object):
    @staticmethod
    def run_module(
            user_input: models.UserInput,
            optimization_cub_score: models.OptimizationCubScore,
            optimization_method: models.OptimizationMethod,
    ) -> str:
        logger.info("##########################")
        logger.info("# ORF #")
        logger.info("##########################\n")
        logger.info(F"Optimization score used is: {optimization_method}")

        target_gene = user_input.sequence

        if optimization_method in (models.OptimizationMethod.hill_climbing_average,
                                   models.OptimizationMethod.hill_climbing_weakest_link):
            return hill_climbing_optimize_by_zscore(
                seq=target_gene,
                user_input=user_input,
                optimization_method=optimization_method,
                optimization_cub_score=optimization_cub_score,
            )

        if optimization_method == models.OptimizationMethod.hill_climbing_bulk_aa_average:
            return hill_climbing_optimize_aa_bulk_by_zscore(
                seq=target_gene,
                user_input=user_input,
                optimization_method=optimization_method,
                optimization_cub_score=optimization_cub_score,
            )

        if optimization_method.is_single_codon_optimization:
            return optimize_sequence_by_loss_function(target_gene=target_gene,
                                                      organisms=user_input.organisms,
                                                      optimization_cub_score=optimization_cub_score,
                                                      optimization_method=optimization_method,
                                                      tuning_param=user_input.tuning_parameter)

        raise ValueError('optimization method is invalid')
