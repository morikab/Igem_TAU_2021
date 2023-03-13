import typing

from logger_factory.logger_factory import LoggerFactory
from numpy import average
from modules import models
from modules.ORF.calculating_cai import general_geomean


logger = LoggerFactory.get_logger()


class OptimizationModule(object):
    @classmethod
    def run_module(cls,
                   final_seq: str,
                   user_input: models.UserInput,
                   optimization_method: models.OptimizationMethod,
                   optimization_cub_score: models.OptimizationCubScore):
        optimization_cub_score_value = optimization_cub_score.value.tolower()

        std_key = F"{optimization_cub_score_value}_std"
        average_key = F"{optimization_cub_score_value}_avg"
        weights = F"{optimization_cub_score_value}_profile"

        optimized_organisms_scores = []
        optimized_organisms_weights = []
        deoptimized_organisms_scores = []
        deoptimized_organisms_weights = []
        # todo: add something related to the ratio between the two worst organisms
        for organism in user_input.organisms:
            sigma = getattr(organism, std_key)
            miu = getattr(organism, average_key)
            profile = getattr(organism, weights)
            index = general_geomean([user_input.sequence, final_seq], weights=profile)
            final_score = index[1]
            organism_score = (final_score - miu) / sigma
            logger.info(F"CUB score for organism {organism.name} is: {final_score}")
            if organism.is_optimized:
                optimized_organisms_scores.append(organism_score)
                optimized_organisms_weights.append(organism.optimization_priority)
            else:
                deoptimized_organisms_scores.append(organism_score)
                deoptimized_organisms_weights.append(organism.optimization_priority)

        alpha = user_input.tuning_parameter
        if optimization_method == models.OptimizationMethod.zscore_single_aa_average:
            return cls._calculate_average_score(optimized_organisms_scores=optimized_organisms_scores,
                                                deoptimized_organisms_scores=deoptimized_organisms_scores,
                                                optimized_organisms_weights=optimized_organisms_weights,
                                                deoptimized_organisms_weights=deoptimized_organisms_weights,
                                                alpha=alpha)

        if optimization_method == models.OptimizationMethod.zscore_single_aa_weakest_link:
            return cls._calculate_weakest_link_score(
                optimized_organisms_scores=optimized_organisms_scores,
                deoptimized_organisms_scores=deoptimized_organisms_scores,
                alpha=alpha,
            )

        raise NotImplementedError(F"Optimization method: {optimization_method}")

    @staticmethod
    def _calculate_average_score(
            optimized_organisms_scores: typing.List[float],
            deoptimized_organisms_scores: typing.List[float],
            optimized_organisms_weights: typing.List[float],
            deoptimized_organisms_weights: typing.List[float],
            alpha: float,
    ) -> float:
        mean_opt_index = average(optimized_organisms_scores, weights=optimized_organisms_weights)
        mean_deopt_index = average(deoptimized_organisms_scores, weights=deoptimized_organisms_weights)
        # norm_factor = max(mean_opt_index, - mean_deopt_index)
        return alpha * mean_opt_index - (1 - alpha) * mean_deopt_index  # /norm_factor

    @staticmethod
    def _calculate_weakest_link_score(
            optimized_organisms_scores: typing.List[float],
            deoptimized_organisms_scores: typing.List[float],
            alpha: float,
    ) -> float:
        return alpha * min(optimized_organisms_scores)-(1-alpha) * max(deoptimized_organisms_scores)
