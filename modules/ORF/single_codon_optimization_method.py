import typing
from collections import defaultdict

from logger_factory.logger_factory import LoggerFactory
from modules.configuration import Configuration
from modules import models
from modules import shared_functions_and_vars

logger = LoggerFactory.get_logger()
config = Configuration.get_config()


# --------------------------------------------------------------
def optimize_sequence(target_gene: str,
                      organisms: typing.Sequence[models.Organism],
                      optimization_method: models.OptimizationMethod,
                      optimization_cub_index: models.OptimizationCubIndex,
                      tuning_param: float) -> str:
    """
    The function calculates the difference between the features of each codon.

    :return: Optimized gene sequence according to the organisms' features
    """
    aa_to_optimal_codon_mapping = find_optimal_codons(organisms=organisms,
                                                      tuning_param=tuning_param,
                                                      optimization_method=optimization_method,
                                                      optimization_cub_index=optimization_cub_index)
    target_protein = shared_functions_and_vars.translate(target_gene)
    optimized_sequence = "".join([aa_to_optimal_codon_mapping[aa] for aa in target_protein])

    return optimized_sequence


# --------------------------------------------------------------
# def loss_function_old(organisms: typing.Sequence[models.Organism],
#                       codons: typing.Sequence[str],
#                       tuning_param: float,
#                       local_maximum: bool,
#                       is_ratio: bool) -> typing.Dict[str, float]:
#     """
#     The function iterates through each feature in each organism, and sums up loss for each codon.
#     It returns a mapping from codon to its score.
#     """
#     loss = {}
#     iterate_through_feature_method = partial(iterate_through_feature,
#                                              codons=codons,
#                                              loss=loss,
#                                              tuning_param=tuning_param,
#                                              is_ratio=is_ratio)
#     if local_maximum:
#         for organism in organisms:
#             # TODO - validate if the new logic gets same results as the new implementation
#             # TODO - check if the loss is indeed updated with each iteration when using partial methods
#             loss = iterate_through_feature_method([organism], high_expression=organism.is_optimized)
#         # for high_expression_organism in high_expression_organisms:
#         #     loss = iterate_through_feature_method([high_expression_organism], high_expression=True)
#         #
#         # for low_expression_organism in low_expression_organisms:
#         #     loss = iterate_through_feature_method([low_expression_organism], high_expression=False)
#     else:
#         high_expression_organisms = [organism for organism in organisms if organism.is_optimized]
#         low_expression_organisms = [organism for organism in organisms if not organism.is_optimized]
#         loss = iterate_through_feature_method(high_expression_organisms, high_expression=True)
#         loss = iterate_through_feature_method(low_expression_organisms, high_expression=False)
#
#     return loss


# --------------------------------------------------------------
def get_organism_attribute_name_by_optimization_cub_index(optimization_cub_index: models.OptimizationCubIndex) -> str:
    if optimization_cub_index.is_codon_adaptation_score:
        return models.Organism.CAI_PROFILE_ATTRIBUTE_NAME
    elif optimization_cub_index.is_trna_adaptation_score:
        return models.Organism.TAI_PROFILE_ATTRIBUTE_NAME
    else:
        raise ValueError(F"Unknown optimization_cub_score {optimization_cub_index}")


# --------------------------------------------------------------
def get_max_organisms_attribute_value(
    organisms: typing.Sequence[models.Organism],
    codons: typing.Sequence[str],
    organism_attribute_name: str,
) -> float:
    all_codon_attribute_values = []

    # TODO - continue from here ...
    for organism in organisms:
        organism_attribute_values = getattr(organism, organism_attribute_name)
        all_codon_attribute_values.extend([organism_attribute_values[codon] for codon in codons])

    max_value = max(all_codon_attribute_values)

    if max_value == 0:
        max_value = 0.000001

    return max_value
# --------------------------------------------------------------

# TODO - delete this method and use a single get_max_organisms_attribute_value
def get_max_value_for_organism(organism: models.Organism,
                               organism_attribute_name: str,
                               codons: typing.Sequence[str],
                               optimization_method: models.OptimizationMethod,
                               max_value_for_wanted_hosts: float,
                               max_value_for_unwanted_hosts: float) -> float:
    if optimization_method in (models.OptimizationMethod.single_codon_diff,
                               models.OptimizationMethod.single_codon_ratio):
        return max_value_for_wanted_hosts if organism.is_optimized else max_value_for_unwanted_hosts

    return get_max_organisms_attribute_value(organisms=[organism],
                                             codons=codons,
                                             organism_attribute_name=organism_attribute_name)


# --------------------------------------------------------------
def calculate_organism_loss_per_codon(organism: models.Organism,
                                      codon: str,
                                      max_value: float,
                                      tuning_param: float,
                                      optimization_method: models.OptimizationMethod,
                                      organism_attribute_name: str) -> float:
    organism_codon_weight = getattr(organism, organism_attribute_name)[codon]

    def _optimized_organism_diff_based_loss_function() -> float:
        return tuning_param * ((max_value - organism_codon_weight) ** 2)

    def _deoptimized_organism_diff_based_loss_function() -> float:
        return (1 - tuning_param) * ((1 - max_value + organism_codon_weight) ** 2)

    def _optimized_organism_ratio_based_loss_function() -> float:
        return tuning_param * ((organism_codon_weight / max_value - 1) ** 2)

    def _deoptimized_organism_ratio_based_loss_function() -> float:
        return (1 - tuning_param) * ((organism_codon_weight / max_value) ** 2)

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
    loss = defaultdict(int)
    organism_attribute_name = get_organism_attribute_name_by_optimization_cub_index(optimization_cub_index)
    max_value_for_high_expression_organisms = get_max_organisms_attribute_value(
        organisms=[organism for organism in organisms if organism.is_optimized],
        codons=codons,
        organism_attribute_name=organism_attribute_name,
    )
    max_value_for_low_expression_organisms = get_max_organisms_attribute_value(
        organisms=[organism for organism in organisms if not organism.is_optimized],
        codons=codons,
        organism_attribute_name=organism_attribute_name,
    )

    for organism in organisms:
        max_value = get_max_value_for_organism(organism=organism,
                                               codons=codons,
                                               organism_attribute_name=organism_attribute_name,
                                               optimization_method=optimization_method,
                                               max_value_for_wanted_hosts=max_value_for_high_expression_organisms,
                                               max_value_for_unwanted_hosts=max_value_for_low_expression_organisms)
        for codon in codons:
            loss[codon] += calculate_organism_loss_per_codon(organism=organism,
                                                             codon=codon,
                                                             max_value=max_value,
                                                             tuning_param=tuning_param,
                                                             optimization_method=optimization_method,
                                                             organism_attribute_name=organism_attribute_name)

    return loss


# --------------------------------------------------------------
# def iterate_through_feature(organisms: typing.Sequence[models.Organism],
#                             codons: typing.Sequence[str],
#                             loss: typing.Dict[str, float],
#                             tuning_param: float,
#                             high_expression: bool,
#                             is_ratio: bool) -> typing.Dict[str, float]:
#     """
#     The function calculates loss for each codon for the organism and adds it to the loss
#     from the previously calculated organisms.
#     It returns an updated loss dictionary for each codon.
#     """
#     # 1. iterate all features
#     # 2. Find the max value for the feature
#     #   2.1 Local - organisms is just one organism at a time
#     #   2.2 Global - organisms is all the optimized/deoptimized organisms given for the run
#
#     for feature_name in [feature.index_name for feature in organisms[0].features]:
#         new_loss = defaultdict(int)
#         max_value = find_max_value_per_feature(organisms, feature_name, codons)
#         logger.info(F"Max feature value is: {max_value}")
#         for organism in organisms:
#             feature = [feature for feature in organism.features if feature.index_name == feature_name]
#             f = feature[0]
#             logger.info(F"Feature ratio is: {f.ratio}")
#
#             for codon in codons:
#                 # Changed in order not to override the codon loss score in each iteration with the value of the low expression codons
#                 # loss[codon] = 0
#                 try:    # todo: temporal change. When synonymous codons dict is done, erase 'try-except'
#                     # optimized organisms should have small loss
#                     if high_expression:
#                         # loss[codon] += (tuning_param * f.ratio * ((f.weights[codon] / max_value - 1) ** 2))
#                         if is_ratio:
#                             new_loss[codon] += tuning_param * f.ratio * ((f.weights[codon] / max_value - 1) ** 2)
#                         else:
#                             new_loss[codon] += tuning_param * f.ratio * (max_value - f.weights[codon])
#                     else:
#                         # loss[codon] += (1 - tuning_param) * f.ratio * ((f.weights[codon] / max_value) ** 2)
#                         if is_ratio:
#                             new_loss[codon] += (1 - tuning_param) * f.ratio * ((f.weights[codon] / max_value) ** 2)
#                         else:
#                             new_loss[codon] += (1 - tuning_param) * f.ratio * (1 - max_value + f.weights[codon])
#                 except:
#                     continue
#
#         for codon in new_loss.keys():
#             if codon in loss:
#                 loss[codon] += new_loss[codon]
#             else:
#                 loss[codon] = new_loss[codon]
#
#     return loss
# # --------------------------------------------------------------
#
#
# def find_max_value_per_feature(organisms, feature_name, codons):
#     values = []
#     for organism in organisms:
#         for feature in organism.features:
#             if feature.index_name == feature_name:
#                 try:    # todo: temporal change. When synonymous codons dict is done, erase 'try-except'
#                     values.extend([feature.weights[codon] for codon in codons])
#                 except:
#                     values.append(0)
#     max_value = max(values)
#
#     if max_value == 0:
#         max_value = 0.000001
#
#     return max_value
# --------------------------------------------------------------
def find_optimal_codons(organisms: typing.Sequence[models.Organism],
                        tuning_param: float,
                        optimization_method: models.OptimizationMethod,
                        optimization_cub_index: models.OptimizationCubIndex,
                        evaluation_function: typing.Callable = loss_function) -> typing.Dict[str, str]:
    """
    :return: Dictionary in the format Amino Acid: Optimal codon.
    """
    optimal_codons = {}

    for aa, codons in shared_functions_and_vars.synonymous_codons.items():
        # TODO - remove the proxy function pointer.
        loss = evaluation_function(organisms=organisms,
                                   codons=codons,
                                   tuning_param=tuning_param,
                                   optimization_method=optimization_method,
                                   optimization_cub_index=optimization_cub_index)
        logger.info(F"Loss dict is: {loss}")
        optimal_codons[aa] = min(loss, key=loss.get)
        logger.info(F"Optimal codon for {aa} is: {optimal_codons[aa]}")

    return optimal_codons


# --------------------------------------------------------------
def optimize_initiation(seq: str) -> str:
    """
    For now, the function returns the same sequence as in the input.
    """
    return seq
