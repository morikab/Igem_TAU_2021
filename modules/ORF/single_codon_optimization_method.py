import typing
from collections import defaultdict
from numpy import average

from logger_factory.logger_factory import LoggerFactory
from modules import models
from modules import shared_functions_and_vars
from modules.configuration import Configuration
from modules.run_summary import RunSummary

logger = LoggerFactory.get_logger()
config = Configuration.get_config()


# --------------------------------------------------------------
def optimize_sequence(target_gene: str,
                      organisms: typing.Sequence[models.Organism],
                      optimization_method: models.OptimizationMethod,
                      optimization_cub_index: models.OptimizationCubIndex,
                      tuning_param: float) -> str:
    aa_to_optimal_codon_mapping = _find_optimal_codons(organisms=organisms,
                                                       tuning_param=tuning_param,
                                                       optimization_method=optimization_method,
                                                       optimization_cub_index=optimization_cub_index)

    target_protein = shared_functions_and_vars.translate(target_gene)
    optimized_sequence = "".join([aa_to_optimal_codon_mapping[aa] for aa in target_protein])

    orf_summary = {
        "aa_to_optimal_codon": aa_to_optimal_codon_mapping,
    }
    RunSummary.add_to_run_summary("orf", orf_summary)

    return optimized_sequence


# --------------------------------------------------------------
def _get_organism_attribute_name_by_optimization_cub_index(optimization_cub_index: models.OptimizationCubIndex) -> str:
    if optimization_cub_index.is_codon_adaptation_score:
        return models.Organism.CAI_PROFILE_ATTRIBUTE_NAME
    elif optimization_cub_index.is_trna_adaptation_score:
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
    codon_attribute_values = [all_organism_attribute_values[codon] for codon in codons]

    max_value = max(codon_attribute_values)
    if max_value == 0:
        max_value = 0.000001

    return max_value


# --------------------------------------------------------------
def _calculate_organism_loss_per_codon(organism: models.Organism,
                                       codon: str,
                                       max_value: float,
                                       tuning_param: float,
                                       optimization_method: models.OptimizationMethod,
                                       organism_attribute_name: str) -> float:
    organism_codon_weight = getattr(organism, organism_attribute_name)[codon]

    def _optimized_organism_diff_based_loss_function() -> float:
        return (max_value - organism_codon_weight) ** 2

    def _deoptimized_organism_diff_based_loss_function() -> float:
        return (1 - max_value + organism_codon_weight) ** 2

    def _optimized_organism_ratio_based_loss_function() -> float:
        return (organism_codon_weight / max_value - 1) ** 2

    def _deoptimized_organism_ratio_based_loss_function() -> float:
        return (organism_codon_weight / max_value) ** 2

    optimization_method_to_loss_function_for_optimized_organisms = {
        models.OptimizationMethod.single_codon_diff: _optimized_organism_diff_based_loss_function,
        models.OptimizationMethod.single_codon_ratio: _optimized_organism_ratio_based_loss_function,
    }
    optimization_method_to_loss_function_for_deoptimized_organisms = {
        models.OptimizationMethod.single_codon_diff: _deoptimized_organism_diff_based_loss_function,
        models.OptimizationMethod.single_codon_ratio: _deoptimized_organism_ratio_based_loss_function,
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
        for organism in organisms:
            max_value = _get_max_organism_attribute_value(organism=organism,
                                                          codons=codons,
                                                          organism_attribute_name=organism_attribute_name)
            organism_loss = _calculate_organism_loss_per_codon(organism=organism,
                                                               codon=codon,
                                                               max_value=max_value,
                                                               tuning_param=tuning_param,
                                                               optimization_method=optimization_method,
                                                               organism_attribute_name=organism_attribute_name)
            if organism.is_optimized:
                optimized_organisms_loss.append(organism_loss)
                optimized_organisms_weights.append(organism.optimization_priority)
            else:
                deoptimized_organisms_loss.append(organism_loss)
                deoptimized_organisms_weights.append(organism.optimization_priority)

        loss_per_codon[codon] = _calculate_total_loss_per_codon(
            optimized_organisms_loss=optimized_organisms_loss,
            deoptimized_organisms_loss=deoptimized_organisms_loss,
            optimized_organisms_weights=optimized_organisms_weights,
            deoptimized_organisms_weights=deoptimized_organisms_weights,
            tuning_parameter=tuning_param,
        )

    return loss_per_codon


# --------------------------------------------------------------
def _calculate_total_loss_per_codon(optimized_organisms_loss: typing.List[float],
                                    deoptimized_organisms_loss: typing.List[float],
                                    optimized_organisms_weights: typing.List[float],
                                    deoptimized_organisms_weights: typing.List[float],
                                    tuning_parameter: float) -> float:
    mean_opt_index = average(optimized_organisms_loss, weights=optimized_organisms_weights)
    mean_deopt_index = average(deoptimized_organisms_loss, weights=deoptimized_organisms_weights)
    return tuning_parameter * mean_opt_index + (1 - tuning_parameter) * mean_deopt_index


# --------------------------------------------------------------
def _find_optimal_codons(organisms: typing.Sequence[models.Organism],
                         tuning_param: float,
                         optimization_method: models.OptimizationMethod,
                         optimization_cub_index: models.OptimizationCubIndex) -> typing.Dict[str, str]:
    """
    :return: Dictionary in the format Amino Acid: Optimal codon.
    """
    optimal_codons = {}

    for aa, codons in shared_functions_and_vars.synonymous_codons.items():
        loss = loss_function(organisms=organisms,
                             codons=codons,
                             tuning_param=tuning_param,
                             optimization_method=optimization_method,
                             optimization_cub_index=optimization_cub_index)
        logger.info(F"Loss dict is: {loss}")
        optimal_codons[aa] = min(loss, key=loss.get)
        logger.info(F"Optimal codon for {aa} is: {optimal_codons[aa]}")

    return optimal_codons


# --------------------------------------------------------------
def _optimize_initiation(seq: str) -> str:
    """
    For now, the function returns the same sequence as in the input.
    """
    return seq
