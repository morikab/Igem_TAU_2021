import typing

from logger_factory.logger_factory import LoggerFactory
from modules import models
from modules.ORF.optimization import optimize_sequence
from modules.ORF.organism import Organism
from modules.ORF.Liyam_new_optimization_function import hill_climbing_optimize_by_zscore

logger = LoggerFactory.get_logger()

# todo: add a statistical analysis of how close the organisms are- like what is the best codon for eah AA
# and are they close


class ORFModule(object):
    @staticmethod
    def run_module(
            user_input: models.UserInput,
            cai_or_tai: str,
            optimization_method: typing.Optional[models.OptimizationMethod],
            max_iter=50,
    ):
        """
        :param user_input: input from GUI parser.
        :param cai_or_tai: string indicating whether to optimize by cai or tai.
        :param optimization_method: optimization method to use.
        :param max_iter: Number of iterations for hill climbing.
        :return: optimized sequence (Biopython Seq)
        """
        logger.info('##########################')
        logger.info('# ORF #')
        logger.info('##########################')

        optimization_method = optimization_method or models.OptimizationMethod.hill_climbing_average
        target_gene = user_input.sequence
        logger.info(optimization_method)

        if optimization_method in (models.OptimizationMethod.hill_climbing_average,
                                   models.OptimizationMethod.hill_climbing_weakest_link):
            return hill_climbing_optimize_by_zscore(seq=target_gene,
                                                    user_input=user_input,
                                                    cai_or_tai='cai',
                                                    max_iter=max_iter,
                                                    optimization_method=optimization_method)
        input_organisms = user_input.organisms
        # TODO - remove old organism object / remove the method entirely?
        high_expression_organisms = [
            Organism(name=organism.name, tai_weights=organism.tai_profile, cai_weights=organism.cai_profile,
                     feature_to_generate=cai_or_tai, cai_std=organism.cai_std, tai_std=organism.tai_std)
            for organism in input_organisms if organism.is_optimized
        ]

        low_expression_organisms = [
            Organism(name=organism.name, tai_weights=organism.tai_profile, cai_weights=organism.cai_profile,
                     feature_to_generate=cai_or_tai, cai_std=organism.cai_std, tai_std=organism.tai_std)
            for organism in input_organisms if not organism.is_optimized
        ]

        if optimization_method == models.OptimizationMethod.single_codon_global:
            return optimize_sequence(target_gene=target_gene,
                                     high_expression_organisms=high_expression_organisms,
                                     low_expression_organisms=low_expression_organisms,
                                     tuning_param=user_input.tuning_parameter,
                                     local_maximum=False)
        if optimization_method == models.OptimizationMethod.single_codon_local:
            return optimize_sequence(target_gene=target_gene,
                                     high_expression_organisms=high_expression_organisms,
                                     low_expression_organisms=low_expression_organisms,
                                     tuning_param=user_input.tuning_parameter,
                                     local_maximum=True)

        raise ValueError('optimization method is invalid')
