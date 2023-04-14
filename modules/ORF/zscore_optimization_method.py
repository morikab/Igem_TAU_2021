import typing
from collections import defaultdict

from numpy import average

from logger_factory.logger_factory import LoggerFactory
from modules import models
from modules.configuration import Configuration
from modules.run_summary import RunSummary
from modules.shared_functions_and_vars import nt_to_aa
from modules.shared_functions_and_vars import synonymous_codons
from modules.ORF.calculating_cai import general_geomean


logger = LoggerFactory.get_logger()
config = Configuration.get_config()


# --------------------------------------------------------------
# In each round - check all single synonymous codon changes and calculate optimization score - take the best one
def optimize_sequence_by_zscore_single_aa(
        sequence: str,
        user_input: models.UserInput,
        optimization_cub_index: models.OptimizationCubIndex,
        optimization_method: models.OptimizationMethod,
        max_iterations: int = config["ORF"]["ZSCORE_MAX_ITERATIONS"],
):
    """
    Iterative codon optimization:
    In each iteration - for each codon, change all synonymous codons to a specific one and test the zscore of the new
    sequence after each iteration, select the sequence with the best zscore - if it was not changed since the last
    iteration, break. The maximum number of iterations allowed is "max_iter".
    """
    sequence_options = {}

    score = _calculate_zscore_for_sequence(
        sequence=sequence,
        user_input=user_input,
        optimization_cub_index=optimization_cub_index,
        optimization_method=optimization_method,
    )
    sequence_options[sequence] = score
    original_sequence = sequence

    aa_to_codon_mapping = defaultdict(str)
    iterations_count = 0
    # Single codon replacement
    for run in range(max_iterations):
        iterations_count = run + 1
        tested_sequence_to_codon = {}
        for codon in nt_to_aa.keys():
            tested_sequence, _ = _change_all_codons_of_aa(sequence, codon)
            tested_sequence_to_codon[tested_sequence] = codon

            score = _calculate_zscore_for_sequence(
                sequence=tested_sequence,
                user_input=user_input,
                optimization_cub_index=optimization_cub_index,
                optimization_method=optimization_method,
            )
            sequence_options[tested_sequence] = score

        new_sequence = max(sequence_options, key=sequence_options.get)
        selected_codon = tested_sequence_to_codon[new_sequence]
        logger.info(F"selected_codon in iteration {run} is: {selected_codon}")
        aa_to_codon_mapping[nt_to_aa[selected_codon]] = selected_codon

        if new_sequence == sequence:
            break
        else:
            sequence = new_sequence

    orf_summary = {

        "iterations_count": iterations_count,
        "aa_to_optimal_codon": aa_to_codon_mapping,
        "initial_sequence_optimization_score": sequence_options[original_sequence],
        "final_sequence_optimization_score": sequence_options[sequence],
    }

    RunSummary.add_to_run_summary("orf", orf_summary)

    return sequence


# --------------------------------------------------------------
def optimize_sequence_by_zscore_bulk_aa(sequence: str,
                                        user_input: models.UserInput,
                                        optimization_method: models.OptimizationMethod,
                                        optimization_cub_index: models.OptimizationCubIndex,
                                        max_iterations: int = config["ORF"]["ZSCORE_MAX_ITERATIONS"]):
    sequence_options = {}
    original_sequence = sequence

    # # Get codon split in the original sequence
    # codon_counts_in_seq = defaultdict(int)
    # split_seq = [original_sequence[i:i + 3].upper() for i in range(0, len(original_sequence), 3)]
    # for codon in split_seq:
    #     codon_counts_in_seq[codon] += 1
    #
    # logger.info(F"codon_counts_in_seq: {codon_counts_in_seq}")

    score = _calculate_zscore_for_sequence(
        sequence=sequence,
        user_input=user_input,
        optimization_cub_index=optimization_cub_index,
        optimization_method=optimization_method,
    )
    sequence_options[sequence] = score
    aa_to_selected_codon = {}
    iterations_count = 0

    for run in range(max_iterations):
        def _find_best_aa_synonymous_codon(codons_list, sequence_to_change: str) -> str:
            aa_codons_to_score = {}
            for aa_codon in codons_list:
                candidate_codon_sequence, candidate_codon_count = _change_all_codons_of_aa(sequence_to_change, aa_codon)
                # logger.info(F"Number of occurrences of codon {aa_codon} in sequence is {candidate_codon_count}")
                logger.info(F"Running for codon: {aa_codon}")
                aa_codons_to_score[aa_codon] = _calculate_zscore_for_sequence(
                    sequence=candidate_codon_sequence,
                    user_input=user_input,
                    optimization_cub_index=optimization_cub_index,
                    optimization_method=optimization_method,
                )
                # logger.info(F"z-score after changing codon {aa_codon} is: {aa_codons_to_score[aa_codon]}")
            selected_aa_codon = max(aa_codons_to_score, key=aa_codons_to_score.get)
            return selected_aa_codon
        logger.info(F"zscore of sequence in run {run} is: {sequence_options[sequence]}")
        aa_to_selected_codon = {}
        iterations_count = run + 1
        # Find the best synonymous_codon per aa
        for aa in synonymous_codons.keys():
            selected_codon = _find_best_aa_synonymous_codon(codons_list=synonymous_codons[aa],
                                                            sequence_to_change=sequence)
            aa_to_selected_codon[aa] = selected_codon

        logger.info(F"aa_to_selected_codon in iteration {run} is: {aa_to_selected_codon}")

        # create new sequence by replacing all synonymous codons
        new_sequence = sequence
        for aa in aa_to_selected_codon:
            new_sequence, _ = _change_all_codons_of_aa(new_sequence, aa_to_selected_codon[aa])

        # Calculate score after all replacements
        score = _calculate_zscore_for_sequence(
            sequence=new_sequence,
            user_input=user_input,
            optimization_cub_index=optimization_cub_index,
            optimization_method=optimization_method,
        )
        sequence_options[new_sequence] = score
        logger.info(F"New seq for iteration {run} with score of: {score}")

        if new_sequence == sequence:
            break
        else:
            sequence = new_sequence

    orf_summary = {
        "iterations_count": iterations_count,
        "aa_to_optimal_codon": aa_to_selected_codon,
        "initial_sequence_optimization_score": sequence_options[original_sequence],
        "final_sequence_optimization_score": sequence_options[sequence],
    }

    RunSummary.add_to_run_summary("orf", orf_summary)

    return sequence


# --------------------------------------------------------------
def _change_all_codons_of_aa(seq: str, selected_codon: str) -> typing.Tuple[str, int]:
    split_seq = [seq[i:i+3].upper() for i in range(0, len(seq), 3)]
    new_split_seq = []
    changed_codons_count = 0
    for codon in split_seq:
        if nt_to_aa[codon] == nt_to_aa[selected_codon]:
            new_split_seq.append(selected_codon)
            changed_codons_count += 1
        else:
            new_split_seq.append(codon)
    return ''.join(new_split_seq), changed_codons_count


# --------------------------------------------------------------
def _calculate_zscore_for_sequence(sequence: str,
                                   user_input: models.UserInput,
                                   optimization_method: models.OptimizationMethod,
                                   optimization_cub_index: models.OptimizationCubIndex):
    optimization_cub_index_value = optimization_cub_index.value.lower()

    std_key = F"{optimization_cub_index_value}_std"
    average_key = F"{optimization_cub_index_value}_avg"
    weights = F"{optimization_cub_index_value}_profile"

    optimized_organisms_scores = []
    optimized_organisms_weights = []
    deoptimized_organisms_scores = []
    deoptimized_organisms_weights = []

    for organism in user_input.organisms:
        sigma = getattr(organism, std_key)
        miu = getattr(organism, average_key)
        profile = getattr(organism, weights)
        index_score = general_geomean([sequence], weights=profile)[0]
        organism_score = (index_score - miu) / sigma
        logger.info(F"CUB score for organism {organism.name} is: {index_score}")
        if organism.is_optimized:
            optimized_organisms_scores.append(organism_score)
            optimized_organisms_weights.append(organism.optimization_priority)
        else:
            deoptimized_organisms_scores.append(organism_score)
            deoptimized_organisms_weights.append(organism.optimization_priority)

    alpha = user_input.tuning_parameter
    if optimization_method.is_zscore_average_score_optimization:
        return _calculate_average_score(optimized_organisms_scores=optimized_organisms_scores,
                                        deoptimized_organisms_scores=deoptimized_organisms_scores,
                                        optimized_organisms_weights=optimized_organisms_weights,
                                        deoptimized_organisms_weights=deoptimized_organisms_weights,
                                        tuning_parameter=alpha)

    if optimization_method.is_zscore_weakest_link_score_optimization:
        return _calculate_weakest_link_score(
            optimized_organisms_scores=optimized_organisms_scores,
            deoptimized_organisms_scores=deoptimized_organisms_scores,
            optimized_organisms_weights=optimized_organisms_weights,
            deoptimized_organisms_weights=deoptimized_organisms_weights,
            tuning_parameter=alpha,
        )

    raise NotImplementedError(F"Optimization method: {optimization_method}")


# --------------------------------------------------------------
def _calculate_average_score(
        optimized_organisms_scores: typing.List[float],
        deoptimized_organisms_scores: typing.List[float],
        optimized_organisms_weights: typing.List[float],
        deoptimized_organisms_weights: typing.List[float],
        tuning_parameter: float,
) -> float:
    mean_opt_index = average(optimized_organisms_scores, weights=optimized_organisms_weights)
    mean_deopt_index = average(deoptimized_organisms_scores, weights=deoptimized_organisms_weights)
    return tuning_parameter * mean_opt_index - (1 - tuning_parameter) * mean_deopt_index


# --------------------------------------------------------------
def _calculate_weakest_link_score(
        optimized_organisms_scores: typing.List[float],
        deoptimized_organisms_scores: typing.List[float],
        optimized_organisms_weights: typing.List[float],
        deoptimized_organisms_weights: typing.List[float],
        tuning_parameter: float,
) -> float:
    weighted_optimized_organisms_scores = [optimized_organisms_scores[i] * optimized_organisms_weights[i] for i in
                                           range(len(optimized_organisms_scores))]
    weighted_deoptimized_organisms_scores = [
        deoptimized_organisms_scores[i] * deoptimized_organisms_weights[i] for
        i in range(len(deoptimized_organisms_scores))
    ]

    return (tuning_parameter * min(weighted_optimized_organisms_scores) -
            (1 - tuning_parameter) * max(weighted_deoptimized_organisms_scores))
