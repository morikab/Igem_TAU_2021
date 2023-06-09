import typing
from collections import defaultdict

from numpy import average

from logger_factory.logger_factory import LoggerFactory
from modules import models
from modules.configuration import Configuration
from modules.run_summary import RunSummary
from modules.shared_functions_and_vars import nt_to_aa
from modules.shared_functions_and_vars import synonymous_codons
from modules.timer import Timer
from modules.ORF.calculating_cai import general_geomean


logger = LoggerFactory.get_logger()
config = Configuration.get_config()


# --------------------------------------------------------------
def optimize_sequence_by_zscore_bulk_and_single_aa(
        sequence: str,
        user_input: models.UserInput,
        optimization_cub_index: models.OptimizationCubIndex,
        optimization_method: models.OptimizationMethod,
        max_iterations: int = config["ORF"]["ZSCORE_MAX_ITERATIONS"],
):
    bulk_sequence = optimize_sequence_by_zscore_bulk_aa(
        sequence=sequence,
        user_input=user_input,
        optimization_method=optimization_method,
        optimization_cub_index=optimization_cub_index,
        max_iterations=max_iterations,
    )
    return optimize_sequence_by_zscore_single_aa(
        sequence=bulk_sequence,
        user_input=user_input,
        optimization_method=optimization_method,
        optimization_cub_index=optimization_cub_index,
        max_iterations=max_iterations,
    )


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

    with Timer() as timer:
        initial_sequence = sequence
        score = _calculate_zscore_for_sequence(
            sequence=sequence,
            user_input=user_input,
            optimization_cub_index=optimization_cub_index,
        )
        initial_sequence_score = None

        aa_to_codon_mapping = defaultdict(str)
        iterations_count = 0
        iterations_summary = []
        # Single codon replacement
        for run in range(max_iterations):
            iterations_count = run + 1
            # In the first iteration, include also the original sequence
            sequence_options = {sequence: score}
            tested_sequence_to_codon = defaultdict(list)
            for codon in nt_to_aa.keys():
                tested_sequence, _ = _change_all_codons_of_aa(sequence, codon)
                tested_sequence_to_codon[tested_sequence].append(codon)

                sequence_options[tested_sequence] = _calculate_zscore_for_sequence(
                    sequence=tested_sequence,
                    user_input=user_input,
                    optimization_cub_index=optimization_cub_index,
                )

            if optimization_method.is_zscore_ratio_score_optimization:
                min_zscore = min(score.min_zscore for score in sequence_options.values())
                max_zscore = max(score.max_zscore for score in sequence_options.values())

                for score in sequence_options.values():
                    score.normalize(min_zscore=min_zscore, max_zscore=max_zscore)

            sequence_to_total_score = {
                sequence_option: get_total_score(zscore=sequence_score,
                                                 optimization_method=optimization_method,
                                                 tuning_parameter=user_input.tuning_parameter) for
                sequence_option, sequence_score in sequence_options.items()
            }
            if initial_sequence_score is None:
                initial_sequence_score = sequence_to_total_score[initial_sequence]

            new_sequence = max(sequence_to_total_score, key=sequence_to_total_score.get)
            selected_codons = tested_sequence_to_codon.get(new_sequence)

            codon_to_score = {}
            for codon_sequence in tested_sequence_to_codon:
                for codon in tested_sequence_to_codon[codon_sequence]:
                    codon_to_score[codon] = sequence_to_total_score[codon_sequence]

            iteration_summary = {
                "selected_codons": [(selected_codon, nt_to_aa[selected_codon]) for selected_codon in selected_codons],
                "codon_to_score": codon_to_score,
                "sequence_score": sequence_to_total_score[new_sequence],
            }
            iterations_summary.append(iteration_summary)

            for selected_codon in selected_codons:
                # If the aa does not appear at all in the cds, this may give a faulty result (that should be consistent)
                aa_to_codon_mapping[nt_to_aa[selected_codon]] = selected_codon

            if new_sequence == sequence:
                break
            else:
                sequence = new_sequence

    orf_summary = {
        "iterations_count": iterations_count,
        "aa_to_optimal_codon": aa_to_codon_mapping,
        "initial_sequence_optimization_score": initial_sequence_score,
        "final_sequence_optimization_score": sequence_to_total_score[sequence],
        "run_time": timer.elapsed_time,
        "iterations_summary": iterations_summary,
    }
    # RunSummary.add_to_run_summary("orf", orf_summary)
    RunSummary.put_in_run_summary("orf", orf_summary)

    return sequence


# --------------------------------------------------------------
# FIXME - continue from here!!!! Fix the code for bulk aa variation as I did for the single aa variation
def optimize_sequence_by_zscore_bulk_aa(sequence: str,
                                        user_input: models.UserInput,
                                        optimization_method: models.OptimizationMethod,
                                        optimization_cub_index: models.OptimizationCubIndex,
                                        max_iterations: int = config["ORF"]["ZSCORE_MAX_ITERATIONS"]):

    with Timer() as timer:
        initial_sequence_score = _calculate_zscore_for_sequence(
            sequence=sequence,
            user_input=user_input,
            optimization_cub_index=optimization_cub_index,
        )
        aa_to_selected_codon = {}
        iterations_count = 0
        iterations_summary = []
        for run in range(max_iterations):
            def _find_best_aa_synonymous_codon(codons_list, sequence_to_change: str) -> str:
                aa_codons_to_score = {}
                for aa_codon in codons_list:
                    candidate_codon_sequence, candidate_codon_count = _change_all_codons_of_aa(sequence_to_change, aa_codon)
                    # logger.info(F"Number of occurrences of codon {aa_codon} in sequence is {candidate_codon_count}")
                    # logger.info(F"Running for codon: {aa_codon}")
                    aa_codons_to_score[aa_codon] = _calculate_zscore_for_sequence(
                        sequence=candidate_codon_sequence,
                        user_input=user_input,
                        optimization_cub_index=optimization_cub_index,
                    )
                    # logger.info(F"z-score after changing codon {aa_codon} is: {aa_codons_to_score[aa_codon]}")
                selected_aa_codon = max(aa_codons_to_score, key=aa_codons_to_score.get)
                return selected_aa_codon

            aa_to_selected_codon = {}
            iterations_count = run + 1
            # Find the best synonymous_codon per aa
            for aa in synonymous_codons.keys():
                selected_codon = _find_best_aa_synonymous_codon(codons_list=synonymous_codons[aa],
                                                                sequence_to_change=sequence)
                aa_to_selected_codon[aa] = selected_codon

            # create new sequence by replacing all synonymous codons
            new_sequence = sequence
            for aa in aa_to_selected_codon:
                new_sequence, _ = _change_all_codons_of_aa(new_sequence, aa_to_selected_codon[aa])

            # Calculate score after all replacements
            score = _calculate_zscore_for_sequence(
                sequence=new_sequence,
                user_input=user_input,
                optimization_cub_index=optimization_cub_index,
            )
            iteration_summary = {
                "aa_to_selected_codon": aa_to_selected_codon,
                "sequence_score": score,
            }
            iterations_summary.append(iteration_summary)

            if new_sequence == sequence:
                break
            else:
                sequence = new_sequence

    orf_summary = {
        "iterations_count": iterations_count,
        "aa_to_optimal_codon": aa_to_selected_codon,
        "initial_sequence_optimization_score": initial_sequence_score,
        "final_sequence_optimization_score": score,
        "run_time": timer.elapsed_time,
        "iterations_summary": iterations_summary,
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
                                   optimization_cub_index: models.OptimizationCubIndex):
    optimization_cub_index_value = optimization_cub_index.value.lower()

    std_key = F"{optimization_cub_index_value}_std"
    average_key = F"{optimization_cub_index_value}_avg"
    weights = F"{optimization_cub_index_value}_profile"

    wanted_hosts_scores = []
    wanted_hosts_weights = []
    unwanted_hosts_scores = []
    unwanted_hosts_weights = []

    for organism in user_input.organisms:
        sigma = getattr(organism, std_key)
        miu = getattr(organism, average_key)
        profile = getattr(organism, weights)
        index_score = general_geomean([sequence], weights=profile)[0]
        organism_score = (index_score - miu) / sigma
        if organism.is_optimized:
            wanted_hosts_scores.append(organism_score)
            wanted_hosts_weights.append(organism.optimization_priority)
        else:
            unwanted_hosts_scores.append(organism_score)
            unwanted_hosts_weights.append(organism.optimization_priority)

        return models.SequenceZscores(
            wanted_hosts_scores=wanted_hosts_scores,
            wanted_hosts_weights=wanted_hosts_weights,
            unwanted_hosts_scores=unwanted_hosts_scores,
            unwanted_hosts_weights=unwanted_hosts_weights,
        )


# --------------------------------------------------------------
def get_total_score(zscore: models.SequenceZscores,
                    optimization_method: models.OptimizationMethod,
                    tuning_parameter: float) -> float:
    if optimization_method.is_zscore_diff_score_optimization:
        return _calculate_zscore_diff_score(optimized_organisms_scores=zscore.wanted_hosts_scores,
                                            deoptimized_organisms_scores=zscore.unwanted_hosts_scores,
                                            optimized_organisms_weights=zscore.wanted_hosts_weights,
                                            deoptimized_organisms_weights=zscore.unwanted_hosts_weights,
                                            tuning_parameter=tuning_parameter)

    if optimization_method.is_zscore_ratio_score_optimization:
        return _calculate_zscore_ratio_score(optimized_organisms_scores=zscore.wanted_hosts_scores,
                                             deoptimized_organisms_scores=zscore.unwanted_hosts_scores,
                                             optimized_organisms_weights=zscore.wanted_hosts_weights,
                                             deoptimized_organisms_weights=zscore.unwanted_hosts_weights,
                                             tuning_parameter=tuning_parameter)

    if optimization_method.is_zscore_weakest_link_score_optimization:
        return _calculate_zscore_weakest_link_score(
            optimized_organisms_scores=zscore.wanted_hosts_scores,
            deoptimized_organisms_scores=zscore.unwanted_hosts_scores,
            optimized_organisms_weights=zscore.wanted_hosts_weights,
            deoptimized_organisms_weights=zscore.unwanted_hosts_weights,
            tuning_parameter=tuning_parameter,
        )


# --------------------------------------------------------------

def _calculate_zscore_for_sequence_old(sequence: str,
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
        # logger.info(F"CUB score for organism {organism.name} is: {index_score}")
        if organism.is_optimized:
            optimized_organisms_scores.append(organism_score)
            optimized_organisms_weights.append(organism.optimization_priority)
        else:
            deoptimized_organisms_scores.append(organism_score)
            deoptimized_organisms_weights.append(organism.optimization_priority)

    alpha = user_input.tuning_parameter
    if optimization_method.is_zscore_diff_score_optimization:
        return _calculate_zscore_diff_score(optimized_organisms_scores=optimized_organisms_scores,
                                            deoptimized_organisms_scores=deoptimized_organisms_scores,
                                            optimized_organisms_weights=optimized_organisms_weights,
                                            deoptimized_organisms_weights=deoptimized_organisms_weights,
                                            tuning_parameter=alpha)

    if optimization_method.is_zscore_ratio_score_optimization:
        return _calculate_zscore_ratio_score(optimized_organisms_scores=optimized_organisms_scores,
                                             deoptimized_organisms_scores=deoptimized_organisms_scores,
                                             optimized_organisms_weights=optimized_organisms_weights,
                                             deoptimized_organisms_weights=deoptimized_organisms_weights,
                                             tuning_parameter=alpha)

    if optimization_method.is_zscore_weakest_link_score_optimization:
        return _calculate_zscore_weakest_link_score(
            optimized_organisms_scores=optimized_organisms_scores,
            deoptimized_organisms_scores=deoptimized_organisms_scores,
            optimized_organisms_weights=optimized_organisms_weights,
            deoptimized_organisms_weights=deoptimized_organisms_weights,
            tuning_parameter=alpha,
        )


    raise NotImplementedError(F"Optimization method: {optimization_method}")


# --------------------------------------------------------------
def _calculate_zscore_diff_score(
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
def _calculate_zscore_ratio_score(
        optimized_organisms_scores: typing.List[float],
        deoptimized_organisms_scores: typing.List[float],
        optimized_organisms_weights: typing.List[float],
        deoptimized_organisms_weights: typing.List[float],
        tuning_parameter: float,
) -> float:
    mean_opt_index = average(optimized_organisms_scores, weights=optimized_organisms_weights)
    mean_deopt_index = average(deoptimized_organisms_scores, weights=deoptimized_organisms_weights)
    # FIXME - think how to handle: 1. mechane that is 0 2. normalize the result to treat the direction correctly
    return (mean_opt_index ** tuning_parameter) / (mean_deopt_index ** (1 - tuning_parameter))


# --------------------------------------------------------------
def _calculate_zscore_weakest_link_score(
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
