from collections import defaultdict
from functools import partial

from logger_factory.logger_factory import LoggerFactory
from modules.shared_functions_and_vars import translate, synonymous_codons

logger = LoggerFactory.get_logger()


# --------------------------------------------------------------
def optimize_sequence(target_gene,
                      high_expression_organisms,
                      low_expression_organisms,
                      local_maximum,
                      tuning_param,
                      is_ratio,
                      n_initiation_codons=12):
    """

    :param target_gene: Gene object, which is to be optimized
    :param high_expression_organisms: list of Organism objects. The organisms where we want to express the target gene in
    :param low_expression_organisms: list of Organism objects. The organisms where we go not want the expression
    :param tuning_param: a number from 0 to 1, which describes the priority in optimization,
    0 - optimization only for no-expression organisms, 1 - only for expression organisms

    :param n_initiation_codons: number of codons of the sequence which need to be optimized due to the initiation rules
    :return: Optimised gene sequence according to the organisms' features

    The function calculates the difference between the features of each codon. Each feature has its own weight (ratio)
    """
    optimized_sequence = ''
    # optimal codons->dict(AA:codon)
    optimal_codons = find_optimal_codons(high_expression_organisms,
                                         low_expression_organisms,
                                         tuning_param,
                                         local_maximum,
                                         is_ratio=is_ratio)

    target_protein = translate(target_gene)

    # optimize the initiation

    optimized_sequence += optimize_initiation(target_gene[:n_initiation_codons*3])

    for aa in target_protein[n_initiation_codons:]:
        optimized_sequence += optimal_codons[aa]

    return optimized_sequence

# --------------------------------------------------------------
def calc_diff(expression_organisms, low_expression_organisms, codons):
    """

    :param expression_organisms: list of expression organisms
    :param low_expression_organisms: list of no_expression organisms
    :param codons: list of codons to calculate the difference for
    :return: dict (codon:score).

    The function finds the difference between the features for each codon in codons list.
    Each feature has ratio, which is actually a weight of the feature in the optimization.



    The function works only for one expression organism and one no-expression organism
    """

    expression_organisms_features = expression_organisms[0].features
    no_expression_organisms_features = low_expression_organisms[0].features

    diff = {}
    for codon in codons:
        diff[codon] = 0
        for i in range(len(expression_organisms_features)):
            if no_expression_organisms_features[i].weights[codon] == 0:
                diff[codon] += 1000
                continue

            diff[codon] += expression_organisms_features[i].ratio * \
                           (expression_organisms_features[i].weights[codon] / no_expression_organisms_features[i].weights[codon])

        # we need to turn the values upside down, because we are looking for minimal value in find_optimal_codons_function

        diff[codon] = 1/diff[codon]

    return diff


# --------------------------------------------------------------
def loss_function(high_expression_organisms, low_expression_organisms, codons, tuning_param, local_maximum, is_ratio):
    """
    :param high_expression_organisms: list of expression organisms
    :param low_expression_organisms: list of no_expression organisms
    :param codons: list of codons to calculate the loss for
    :param tuning_param:
    :return: loss - dict (codon:score)

    The function iterates through each feature in each organism, and sums up loss for each codon
    """
    loss = {}
    iterate_through_feature_method = partial(iterate_through_feature,
                                             codons=codons,
                                             loss=loss,
                                             tuning_param=tuning_param,
                                             is_ratio=is_ratio)
    # high_expression is optimized, low expression is deoptimized
    if local_maximum:
        for high_expression_organism in high_expression_organisms:
            loss = iterate_through_feature_method([high_expression_organism], high_expression=True)

        for low_expression_organism in low_expression_organisms:
            loss = iterate_through_feature_method([low_expression_organism], high_expression=False)
    else:
        loss = iterate_through_feature_method(high_expression_organisms, high_expression=True)
        loss = iterate_through_feature_method(low_expression_organisms, loss=loss, high_expression=False)

    return loss


# --------------------------------------------------------------
def iterate_through_feature(organisms, codons, loss, tuning_param, high_expression, is_ratio):
    """
    :param organisms: List of organism objects for which the sequence is optimized
    :param codons: list of codons to choose from
    :param loss: loss(dict) taken from a previous iteration
    :param high_expression: Whether current organism is being optimized for
    high expression or low expression
    :return: updated loss dictionary

    The funciton calculates loss for each codon for the organism,
    and adds it to the loss from the previously calculated organisms.
    """

    for feature_name in [feature.index_name for feature in organisms[0].features]:
        new_loss = defaultdict(int)
        max_value = find_max_value_per_feature(organisms, feature_name, codons)
        logger.info(F"Max feature value is: {max_value}")
        for organism in organisms:
            feature = [feature for feature in organism.features if feature.index_name == feature_name]
            f = feature[0]
            logger.info(F"Feature ratio is: {f.ratio}")

            for codon in codons:
                # Changed in order not to override the codon loss score in each iteration with the value of the low expression codons
                # loss[codon] = 0
                try:    # todo: temporal change. When synonymous codons dict is done, erase 'try-except'
                    # optimized organisms should have small loss
                    if high_expression:
                        # loss[codon] += (tuning_param * f.ratio * ((f.weights[codon] / max_value - 1) ** 2))
                        if is_ratio:
                            new_loss[codon] += tuning_param * f.ratio * ((f.weights[codon] / max_value - 1) ** 2)
                        else:
                            new_loss[codon] += tuning_param * f.ratio * (max_value - f.weights[codon])
                    else:
                        # loss[codon] += (1 - tuning_param) * f.ratio * ((f.weights[codon] / max_value) ** 2)
                        if is_ratio:
                            new_loss[codon] += (1 - tuning_param) * f.ratio * ((f.weights[codon] / max_value) ** 2)
                        else:
                            new_loss[codon] += (1 - tuning_param) * f.ratio * (1 - max_value + f.weights[codon])
                except:
                    continue

        for codon in new_loss.keys():
            if codon in loss:
                loss[codon] += new_loss[codon]
            else:
                loss[codon] = new_loss[codon]

    return loss
# --------------------------------------------------------------


def find_max_value_per_feature(organisms, feature_name, codons):
    values = []
    for organism in organisms:
        for feature in organism.features:
            if feature.index_name == feature_name:
                try:    # todo: temporal change. When synonymous codons dict is done, erase 'try-except'
                    values.extend([feature.weights[codon] for codon in codons])
                except:
                    values.append(0)
    max_value = max(values)

    if max_value == 0:
        max_value = 0.000001

    return max_value
# --------------------------------------------------------------


def find_optimal_codons(high_expression_organisms,
                        low_expression_organisms,
                        tuning_param,
                        local_maximum,
                        evaluation_function=loss_function,
                        is_ratio=True):
    """
    :param high_expression_organisms: list of Organism objects. The organisms where we want to express the target gene in
    :param low_expression_organisms: list of Organism objects. The organisms where we do not want the expression in
    :param tuning_param:
    :param evaluation_function: function which evaluates loss
    :return: Dictionary in the format Amino Acid: Optimal codon.
    """

    optimal_codons = {}

    for aa, codons in synonymous_codons.items():
        loss = evaluation_function(high_expression_organisms,
                                   low_expression_organisms,
                                   codons,
                                   tuning_param,
                                   local_maximum=local_maximum,
                                   is_ratio=is_ratio)
        logger.info(F"Loss dict is: {loss}")
        optimal_codons[aa] = min(loss, key=loss.get)
        logger.info(F"optimal codon for {aa} is: {optimal_codons[aa]}")

    return optimal_codons


# --------------------------------------------------------------
def optimize_initiation(seq):
    """

    :param seq: Seq object (only the initiation part)
    :return: optimized initiation sequence

    For now, the function returns the same sequence as in the input, tbd.
    """
    return seq
