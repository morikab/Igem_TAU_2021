from functools import partial

from logger_factory.logger_factory import LoggerFactory
from loss_function_optimization_method import optimize_sequence
from modules import models
from zscore_optimization_method import hill_climbing_optimize_by_zscore as hill_climbing_optimize_by_zscore
from zscore_optimization_method import hill_climbing_optimize_aa_bulk_by_zscore as \
    hill_climbing_optimize_aa_bulk_by_zscore

logger = LoggerFactory.get_logger()

# todo: add a statistical analysis of how close the organisms are - like what is the best codon for eah AA
# and are they close


class ORFModule(object):
    @staticmethod
    def run_module(
            user_input: models.UserInput,
            optimization_cub_score: models.OptimizationCubScore,
            optimization_method: models.OptimizationMethod,
            max_iter=50,
    ):
        logger.info('##########################')
        logger.info('# ORF #')
        logger.info('##########################')
        logger.info(F"Optimization score used is: {optimization_method}")

        target_gene = user_input.sequence

        if optimization_method in (models.OptimizationMethod.hill_climbing_average,
                                   models.OptimizationMethod.hill_climbing_weakest_link):
            return hill_climbing_optimize_by_zscore(
                seq=target_gene,
                user_input=user_input,
                max_iter=max_iter,
                optimization_method=optimization_method,
                optimization_cub_score=optimization_cub_score,
            )

        if optimization_method == models.OptimizationMethod.hill_climbing_bulk_aa_average:
            return hill_climbing_optimize_aa_bulk_by_zscore(
                seq=target_gene,
                user_input=user_input,
                max_iter=max_iter,
                optimization_method=optimization_method,
                optimization_cub_score=optimization_cub_score,
            )
        # TODO - &&&&&&&&&&&&&& continue from here # TODO - &&&&&&&&&&&&&& -
        #  1. Use one list of organisms instead multiple lists to make the code simpler
        #  2. Simplify the code to NOT work with multiple features
        #  3. Better create the split between local/global variations and the is_ratio parameter
        optimize_sequence_by_loss_function_method = partial(optimize_sequence,
                                                            target_gene=target_gene,
                                                            organisms=user_input.organisms,
                                                            tuning_param=user_input.tuning_parameter)

        if optimization_method == models.OptimizationMethod.single_codon_global_ratio:
            return optimize_sequence_by_loss_function_method(local_maximum=False, is_ratio=True)
        if optimization_method == models.OptimizationMethod.single_codon_local_ratio:
            return optimize_sequence_by_loss_function_method(local_maximum=True, is_ratio=True)
        if optimization_method == models.OptimizationMethod.single_codon_global_diff:
            return optimize_sequence_by_loss_function_method(local_maximum=False, is_ratio=False)
        if optimization_method == models.OptimizationMethod.single_codon_local_diff:
            return optimize_sequence_by_loss_function_method(local_maximum=True, is_ratio=False)

        raise ValueError('optimization method is invalid')
