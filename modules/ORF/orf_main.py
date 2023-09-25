import typing
from functools import partial

from logger_factory.logger_factory import LoggerFactory
from modules import models
from modules.configuration import Configuration
from modules.run_summary import RunSummary
from modules.shared_functions_and_vars import synonymous_codon_permutation

from .single_codon_optimization_method import optimize_sequence as optimize_sequence_by_single_codon
from .zscore_optimization_method import optimize_sequence_by_zscore_single_aa
from .zscore_optimization_method import optimize_sequence_by_zscore_bulk_aa

config = Configuration.get_config()
logger = LoggerFactory.get_logger()


class ORFModule(object):
    @staticmethod
    def run_module(
            user_input: models.UserInput,
            optimization_cub_index: models.OptimizationCubIndex,
            optimization_method: models.OptimizationMethod,
            run_summary: RunSummary,
    ) -> typing.Sequence[str]:
        logger.info("##########################")
        logger.info("# ORF #")
        logger.info("##########################\n")
        logger.info(F"Optimization method used is: {optimization_method}")

        if optimization_method.is_single_codon_optimization:
            return optimize_sequence_by_single_codon(target_gene=user_input.sequence,
                                                     organisms=user_input.organisms,
                                                     optimization_cub_index=optimization_cub_index,
                                                     optimization_method=optimization_method,
                                                     tuning_param=user_input.tuning_parameter,
                                                     run_summary=run_summary)

        if optimization_method.is_zscore_optimization:
            return ORFModule.optimize_sequence_by_zscore(
                user_input=user_input,
                optimization_cub_index=optimization_cub_index,
                optimization_method=optimization_method,
                run_summary=run_summary,
            )

        raise ValueError(F"optimization method {optimization_method} is invalid")

    @staticmethod
    def optimize_sequence_by_zscore(
            user_input: models.UserInput,
            optimization_cub_index: models.OptimizationCubIndex,
            optimization_method: models.OptimizationMethod,
            run_summary: RunSummary,
    ):
        original_sequence = user_input.sequence
        target_genes = [original_sequence] + [
            synonymous_codon_permutation(original_sequence) for _ in
            range(config["ORF"]["ZSCORE_INITIAL_PERMUTATIONS_NUM"])
        ]
        zscore_optimization = ORFModule.get_zscore_optimization_method(optimization_method)
        partial_zscore_optimization_method = partial(
            zscore_optimization,
            user_input=user_input,
            optimization_method=optimization_method,
            optimization_cub_index=optimization_cub_index,
            run_summary=run_summary,
        )

        results = []
        for target_gene in target_genes:
            results.append(partial_zscore_optimization_method(sequence=target_gene))

        return results      # TODO - continue from here!!!!! , modify evaluation scheme to take the best option s

    @staticmethod
    def get_zscore_optimization_method(optimization_method: models.OptimizationMethod) -> typing.Callable:
        if optimization_method.is_zscore_single_aa_optimization:
            return optimize_sequence_by_zscore_single_aa
        return optimize_sequence_by_zscore_bulk_aa
