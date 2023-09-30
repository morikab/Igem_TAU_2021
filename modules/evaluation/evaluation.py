import typing

from numpy import average
from modules import models as main_models
from modules.run_summary import RunSummary
from modules.ORF.calculating_cai import general_geomean

from . import models


class EvaluationModule(object):
    @staticmethod
    def run_module(final_sequence: str,
                   user_input: main_models.UserInput,
                   optimization_cub_index: main_models.OptimizationCubIndex,
                   run_summary: RunSummary) -> models.EvaluationModuleResult:
        optimization_cub_index_value = optimization_cub_index.value.lower()
        std = f"{optimization_cub_index_value}_std"
        weights = f"{optimization_cub_index_value}_profile"

        optimized_organisms_scores = []
        optimized_organisms_non_normalized_scores = []
        optimized_organisms_weights = []
        deoptimized_organisms_scores = []
        deoptimized_organisms_non_normalized_scores = []
        deoptimized_organisms_weights = []

        organisms_evaluation_summary = []
        for organism in user_input.organisms:
            sigma = getattr(organism, std)
            profile = getattr(organism, weights)
            index = general_geomean([user_input.sequence, final_sequence], weights=profile)
            initial_score = index[0]
            final_score = index[1]
            organism_score = (final_score - initial_score) / sigma
            if organism.is_optimized:
                optimized_organisms_scores.append(organism_score)
                optimized_organisms_non_normalized_scores.append(final_score - initial_score)
                optimized_organisms_weights.append(organism.optimization_priority)
            else:
                deoptimized_organisms_scores.append(organism_score)
                deoptimized_organisms_non_normalized_scores.append(final_score - initial_score)
                deoptimized_organisms_weights.append(organism.optimization_priority)

            organism_summary = {
                "name": organism.name,
                "is_wanted": organism.is_optimized,
                F"{optimization_cub_index_value}_initial_score": initial_score,
                F"{optimization_cub_index_value}_final_score": final_score,
                "dist_score": organism_score,
            }
            organisms_evaluation_summary.append(organism_summary)

        average_distance_score = EvaluationModule._calculate_average_distance_score(
            optimized_organisms_scores=optimized_organisms_scores,
            deoptimized_organisms_scores=deoptimized_organisms_scores,
            optimized_organisms_weights=optimized_organisms_weights,
            deoptimized_organisms_weights=deoptimized_organisms_weights,
            tuning_parameter=user_input.tuning_parameter,
        )

        average_distance_non_normalized_score = EvaluationModule._calculate_average_distance_score(
            optimized_organisms_scores=optimized_organisms_non_normalized_scores,
            deoptimized_organisms_scores=deoptimized_organisms_non_normalized_scores,
            optimized_organisms_weights=optimized_organisms_weights,
            deoptimized_organisms_weights=deoptimized_organisms_weights,
            tuning_parameter=user_input.tuning_parameter,
        )

        weakest_link_score = EvaluationModule._calculate_weakest_link_score(
            optimized_organisms_scores=optimized_organisms_scores,
            deoptimized_organisms_scores=deoptimized_organisms_scores,
            optimized_organisms_weights=optimized_organisms_weights,
            deoptimized_organisms_weights=deoptimized_organisms_weights,
            tuning_parameter=user_input.tuning_parameter,
        )

        evaluation_result = models.EvaluationModuleResult(
            sequence=final_sequence,
            average_distance_score=average_distance_score,
            average_distance_non_normalized_score=average_distance_non_normalized_score,
            weakest_link_score=weakest_link_score,
        )

        evaluation_summary = {
            "organisms": organisms_evaluation_summary,
            **evaluation_result.summary,
        }
        run_summary.add_to_run_summary("evaluation", evaluation_summary)

        return evaluation_result

    # --------------------------------------------------------------
    @staticmethod
    def _calculate_average_distance_score(optimized_organisms_scores: typing.Sequence[float],
                                          deoptimized_organisms_scores: typing.Sequence[float],
                                          optimized_organisms_weights: typing.Sequence[float],
                                          deoptimized_organisms_weights: typing.Sequence[float],
                                          tuning_parameter: float) -> float:
        mean_opt_index = average(optimized_organisms_scores, weights=optimized_organisms_weights)
        mean_deopt_index = average(deoptimized_organisms_scores, weights=deoptimized_organisms_weights)
        return tuning_parameter * mean_opt_index - (1-tuning_parameter) * mean_deopt_index

    # --------------------------------------------------------------
    @staticmethod
    def _calculate_weakest_link_score(optimized_organisms_scores: typing.Sequence[float],
                                      deoptimized_organisms_scores: typing.Sequence[float],
                                      optimized_organisms_weights: typing.Sequence[float],
                                      deoptimized_organisms_weights: typing.Sequence[float],
                                      tuning_parameter: float) -> float:
        weighted_optimized_organisms_scores = [
            optimized_organisms_scores[i] * optimized_organisms_weights[i] for i in
            range(len(optimized_organisms_scores))
        ]
        weighted_deoptimized_organisms_scores = [
            deoptimized_organisms_scores[i] * deoptimized_organisms_weights[i] for i in
            range(len(deoptimized_organisms_scores))
        ]
        return tuning_parameter * min(weighted_optimized_organisms_scores) - (1-tuning_parameter) * max(
            weighted_deoptimized_organisms_scores
        )
