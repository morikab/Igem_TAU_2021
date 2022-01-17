from modules import models
from modules.ORF.calculating_cai import general_geomean
from statistics import mean


# TODO - Add additional stats (% optimized, % deoptimized, etc.)
class ZscoreModule(object):
    @staticmethod
    def run_module(final_seq: str, user_input: models.UserInput, optimization_type: str = 'cai'):
        std = optimization_type + '_std'
        weights = optimization_type + '_profile'

        opt_index_org = []
        deopt_index_org = []
        # todo: add something related to the ratio between the two worst organisms
        for organism in user_input.organisms:
            sigma = getattr(organism, std)
            profile = getattr(organism, weights)
            index = general_geomean([user_input.sequence, final_seq], weights=profile)
            initial_score = index[0]
            final_score = index[1]
            index_org = (final_score - initial_score) / sigma
            if organism.is_optimized:
                opt_index_org.append(index_org)
            else:
                deopt_index_org.append(index_org)

        mean_opt_index = mean(opt_index_org)
        mean_deopt_index = mean(deopt_index_org)
        alpha = user_input.tuning_parameter
        optimization_index = (alpha * mean_opt_index - (1-alpha) * mean_deopt_index)  #/norm_factor
        weakest_score = alpha*min(opt_index_org)-(1-alpha)*max(deopt_index_org)

        return mean_opt_index, mean_deopt_index, optimization_index, weakest_score
