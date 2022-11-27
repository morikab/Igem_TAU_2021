from numpy import average
from modules import models as main_models
from modules.ORF.calculating_cai import general_geomean

from . import models


# TODO - Add additional stats (% optimized, % deoptimized, etc.)
class EvaluationModule(object):
    @staticmethod
    def run_module(final_seq: str,
                   user_input: main_models.UserInput,
                   optimization_type: str = "cai") -> models.EvaluationModuleResult:
        std = optimization_type + '_std'
        weights = optimization_type + '_profile'

        optimized_organisms_scores = []
        optimized_organisms_weights = []
        deoptimized_organisms_scores = []
        deoptimized_organisms_weights = []
        # todo: add something related to the ratio between the two worst organisms
        for organism in user_input.organisms:
            sigma = getattr(organism, std)
            profile = getattr(organism, weights)
            index = general_geomean([user_input.sequence, final_seq], weights=profile)
            initial_score = index[0]
            final_score = index[1]
            organism_score = (final_score - initial_score) / sigma
            if organism.is_optimized:
                optimized_organisms_scores.append(organism_score)
                optimized_organisms_weights.append(organism.optimization_priority)
            else:
                deoptimized_organisms_scores.append(organism_score)
                deoptimized_organisms_weights.append(organism.optimization_priority)

        mean_opt_index = average(optimized_organisms_scores, weights=optimized_organisms_weights)
        mean_deopt_index = average(deoptimized_organisms_scores, weights=deoptimized_organisms_weights)
        alpha = user_input.tuning_parameter
        optimization_index = (alpha * mean_opt_index - (1-alpha) * mean_deopt_index)  #/norm_factor
        weakest_score = alpha*min(optimized_organisms_scores)-(1-alpha)*max(deoptimized_organisms_scores)

        return models.EvaluationModuleResult(
            sequence=final_seq,
            mean_opt_index=mean_opt_index,
            mean_deopt_index=mean_deopt_index,
            optimization_index=optimization_index,
            weakest_score=weakest_score,
        )
