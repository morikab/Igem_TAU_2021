import typing
from collections import defaultdict
from numpy import average
from scipy.stats.mstats import gmean

from logger_factory.logger_factory import LoggerFactory
from modules import models
from modules import shared_functions_and_vars
from modules.configuration import Configuration
from modules.run_summary import RunSummary
from modules.timer import Timer

logger = LoggerFactory.get_logger()
config = Configuration.get_config()


# --------------------------------------------------------------
def optimize_sequence(target_gene: str,
                      organisms: typing.Sequence[models.Organism],
                      optimization_method: models.OptimizationMethod,
                      optimization_cub_index: models.OptimizationCubIndex,
                      tuning_param: float,
                      run_summary: RunSummary) -> str:
    with Timer() as timer:
        aa_to_optimal_codon_mapping = _find_optimal_codons(organisms=organisms,
                                                           tuning_param=tuning_param,
                                                           optimization_method=optimization_method,
                                                           optimization_cub_index=optimization_cub_index,
                                                           run_summary=run_summary)

        target_protein = shared_functions_and_vars.translate(target_gene)
        if target_protein.endswith("_") and optimization_cub_index.is_trna_adaptation_index:
            # There is no point in optimizing stop codon by tAI, so leaving the original codon
            aa_to_optimal_codon_mapping["_"] = target_gene[-3:]
        optimized_sequence = "".join([aa_to_optimal_codon_mapping[aa] for aa in target_protein])

    orf_summary = {
        "aa_to_optimal_codon": aa_to_optimal_codon_mapping,
        "run_time": timer.elapsed_time,
    }
    run_summary.add_to_run_summary("orf", orf_summary)

    return optimized_sequence


# --------------------------------------------------------------
def _get_organism_attribute_name_by_optimization_cub_index(optimization_cub_index: models.OptimizationCubIndex) -> str:
    if optimization_cub_index.is_codon_adaptation_index:
        return models.Organism.CAI_PROFILE_ATTRIBUTE_NAME
    elif optimization_cub_index.is_trna_adaptation_index:
        return models.Organism.TAI_PROFILE_ATTRIBUTE_NAME
    else:
        raise ValueError(F"Unknown optimization_cub_score {optimization_cub_index}")


# --------------------------------------------------------------
def _get_max_organism_attribute_value(
        organism: models.Organism,
        codons: typing.Sequence[str],
        organism_attribute_name: str,
) -> float:
    all_organism_attribute_values = getattr(organism, organism_attribute_name)
    # In case of missing entry per codon, assign 0 for the missing weight
    codon_attribute_values = [all_organism_attribute_values.get(codon, 0) for codon in codons]

    max_value = max(codon_attribute_values)
    if max_value == 0:
        max_value = average(list(all_organism_attribute_values.values()))

    return max_value


# --------------------------------------------------------------
def _calculate_organism_loss_per_codon(organism: models.Organism,
                                       codon: str,
                                       max_value: float,
                                       optimization_method: models.OptimizationMethod,
                                       organism_attribute_name: str) -> float:
    organism_codon_weight = getattr(organism, organism_attribute_name).get(codon, 0)

    def _optimized_organism_diff_based_loss_function() -> float:
        return (max_value - organism_codon_weight) ** 2

    def _deoptimized_organism_diff_based_loss_function() -> float:
        return organism_codon_weight ** 2

    def _optimized_organism_ratio_based_loss_function() -> float:
        return (max_value - organism_codon_weight + 1) ** 2

    def _deoptimized_organism_ratio_based_loss_function() -> float:
        return (max_value - organism_codon_weight + 1) ** 2

    def _optimized_organism_weakest_link_based_loss_function() -> float:
        return (max_value - organism_codon_weight) ** 2

    def _deoptimized_organism_weakest_link_based_loss_function() -> float:
        return organism_codon_weight ** 2

    optimization_method_to_loss_function_for_optimized_organisms = {
        models.OptimizationMethod.single_codon_diff: _optimized_organism_diff_based_loss_function,
        models.OptimizationMethod.single_codon_ratio: _optimized_organism_ratio_based_loss_function,
        models.OptimizationMethod.single_codon_weakest_link: _optimized_organism_weakest_link_based_loss_function,
    }
    optimization_method_to_loss_function_for_deoptimized_organisms = {
        models.OptimizationMethod.single_codon_diff: _deoptimized_organism_diff_based_loss_function,
        models.OptimizationMethod.single_codon_ratio: _deoptimized_organism_ratio_based_loss_function,
        models.OptimizationMethod.single_codon_weakest_link: _deoptimized_organism_weakest_link_based_loss_function,
    }

    loss_function_mapping = optimization_method_to_loss_function_for_optimized_organisms if organism.is_optimized else \
        optimization_method_to_loss_function_for_deoptimized_organisms

    if optimization_method not in loss_function_mapping:
        raise ValueError(F"Missing loss function mapping for optimization method: {optimization_method}")

    return loss_function_mapping[optimization_method]()


# --------------------------------------------------------------
def loss_function(organisms: typing.Sequence[models.Organism],
                  codons: typing.Sequence[str],
                  tuning_param: float,
                  optimization_method: models.OptimizationMethod,
                  optimization_cub_index: models.OptimizationCubIndex) -> typing.Dict[str, float]:
    """
    The function iterates through each organism and sums up loss for each codon.
    It returns a mapping from codon to its score.
    """
    loss_per_codon = defaultdict(float)
    organism_attribute_name = _get_organism_attribute_name_by_optimization_cub_index(optimization_cub_index)
    for codon in codons:
        optimized_organisms_loss = []
        optimized_organisms_weights = []
        deoptimized_organisms_loss = []
        deoptimized_organisms_weights = []
        logger.info(f"Loss summary for: {codon}")
        for organism in organisms:
            max_value = _get_max_organism_attribute_value(organism=organism,
                                                          codons=codons,
                                                          organism_attribute_name=organism_attribute_name)
            organism_loss = _calculate_organism_loss_per_codon(organism=organism,
                                                               codon=codon,
                                                               max_value=max_value,
                                                               optimization_method=optimization_method,
                                                               organism_attribute_name=organism_attribute_name)
            logger.info(f"loss for {organism.name}: {organism_loss}")
            # TODO - create dict for the loss per codon for each organism, then include the frequency of each AA in the sequence to show the cancelling effect that eventually leads to the negtive score (in single codon approach).
            if organism.is_optimized:
                optimized_organisms_loss.append(organism_loss)
                optimized_organisms_weights.append(organism.optimization_priority)
            else:
                deoptimized_organisms_loss.append(organism_loss)
                deoptimized_organisms_weights.append(organism.optimization_priority)

        loss_per_codon[codon] = _calculate_total_loss_per_codon(
            optimization_method=optimization_method,
            optimized_organisms_loss=optimized_organisms_loss,
            deoptimized_organisms_loss=deoptimized_organisms_loss,
            optimized_organisms_weights=optimized_organisms_weights,
            deoptimized_organisms_weights=deoptimized_organisms_weights,
            tuning_parameter=tuning_param,
        )

    return loss_per_codon


# --------------------------------------------------------------
def _calculate_total_loss_per_codon(optimization_method: models.OptimizationMethod,
                                    optimized_organisms_loss: typing.List[float],
                                    deoptimized_organisms_loss: typing.List[float],
                                    optimized_organisms_weights: typing.List[float],
                                    deoptimized_organisms_weights: typing.List[float],
                                    tuning_parameter: float) -> float:
    def _diff_total_loss() -> float:
        mean_opt_index = average(optimized_organisms_loss, weights=optimized_organisms_weights)
        mean_deopt_index = average(deoptimized_organisms_loss, weights=deoptimized_organisms_weights)
        return tuning_parameter * mean_opt_index + (1 - tuning_parameter) * mean_deopt_index

    def _ratio_total_loss() -> float:
        mean_opt_index = gmean(optimized_organisms_loss, weights=optimized_organisms_weights)
        mean_deopt_index = gmean(deoptimized_organisms_loss, weights=deoptimized_organisms_weights)

        return (mean_opt_index ** tuning_parameter) / (mean_deopt_index ** (1 - tuning_parameter))

    def _weakest_link_total_loss() -> float:
        weighted_optimized_organisms_scores = [optimized_organisms_loss[i] * optimized_organisms_weights[i] for i in
                                               range(len(optimized_organisms_loss))]
        weighted_deoptimized_organisms_scores = [
            deoptimized_organisms_loss[i] * deoptimized_organisms_weights[i] for
            i in range(len(deoptimized_organisms_loss))
        ]

        return (tuning_parameter * max(weighted_optimized_organisms_scores) +
                (1 - tuning_parameter) * max(weighted_deoptimized_organisms_scores))

    optimization_method_to_total_loss = {
        models.OptimizationMethod.single_codon_diff: _diff_total_loss,
        models.OptimizationMethod.single_codon_ratio: _ratio_total_loss,
        models.OptimizationMethod.single_codon_weakest_link: _weakest_link_total_loss,
    }

    if optimization_method not in optimization_method_to_total_loss:
        raise NotImplementedError(F"Optimization method: {optimization_method}")

    return optimization_method_to_total_loss[optimization_method]()


# --------------------------------------------------------------
def _find_optimal_codons(organisms: typing.Sequence[models.Organism],
                         tuning_param: float,
                         optimization_method: models.OptimizationMethod,
                         optimization_cub_index: models.OptimizationCubIndex,
                         run_summary: RunSummary) -> typing.Dict[str, str]:
    """
    :return: Dictionary in the format Amino Acid: Optimal codon.
    """
    optimal_codons = {}
    codons_loss_score = {}
    for aa, codons in shared_functions_and_vars.synonymous_codons.items():
        loss = loss_function(organisms=organisms,
                             codons=codons,
                             tuning_param=tuning_param,
                             optimization_method=optimization_method,
                             optimization_cub_index=optimization_cub_index)
        logger.info(F"Loss dict is: {loss}")
        optimal_codons[aa] = min(loss, key=loss.get)
        codons_loss_score[aa] = loss
        logger.info(F"Optimal codon for {aa} is: {optimal_codons[aa]}")

    run_summary.add_to_run_summary("orf_debug", codons_loss_score)
    return optimal_codons


# --------------------------------------------------------------
def _optimize_initiation(seq: str) -> str:
    """
    For now, the function returns the same sequence as in the input.
    """
    return seq
