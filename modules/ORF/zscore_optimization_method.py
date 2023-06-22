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


# TODO - need to modify single aa variation with the normalizd ratio variation..
# --------------------------------------------------------------
# In each round - check all single synonymous codon changes and calculate optimization score - take the best one
def optimize_sequence_by_zscore_single_aa(
        sequence: str,
        user_input: models.UserInput,
        optimization_cub_index: models.OptimizationCubIndex,
        optimization_method: models.OptimizationMethod,
        run_summary: RunSummary,
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
        previous_sequence_score = _calculate_zscore_for_sequence(
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
            # Include also the sequence from the previous iteration
            sequence_to_zscore = {sequence: previous_sequence_score}
            tested_sequence_to_codon = defaultdict(list)
            for codon in nt_to_aa.keys():
                tested_sequence, _ = _change_all_codons_of_aa(sequence, codon)
                tested_sequence_to_codon[tested_sequence].append(codon)

                sequence_to_zscore[tested_sequence] = _calculate_zscore_for_sequence(
                    sequence=tested_sequence,
                    user_input=user_input,
                    optimization_cub_index=optimization_cub_index,
                )

            if optimization_method.is_zscore_ratio_score_optimization:
                min_zscore = min(score.min_zscore for score in sequence_to_zscore.values())
                max_zscore = max(score.max_zscore for score in sequence_to_zscore.values())

                for score in sequence_to_zscore.values():
                    score.normalize(min_zscore=min_zscore, max_zscore=max_zscore)

            sequence_to_total_score = {
                sequence_option: get_total_score(zscore=sequence_score,
                                                 optimization_method=optimization_method,
                                                 tuning_parameter=user_input.tuning_parameter) for
                sequence_option, sequence_score in sequence_to_zscore.items()
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
                previous_sequence_score = sequence_to_total_score[sequence]

    orf_summary = {
        "iterations_count": iterations_count,
        "aa_to_optimal_codon": aa_to_codon_mapping,
        "initial_sequence_optimization_score": initial_sequence_score,
        "final_sequence_optimization_score": sequence_to_total_score[sequence],
        "run_time": timer.elapsed_time,
        "iterations_summary": iterations_summary,
    }
    run_summary.add_to_run_summary("orf", orf_summary)

    return sequence


# --------------------------------------------------------------
def optimize_sequence_by_zscore_bulk_aa(sequence: str,
                                        user_input: models.UserInput,
                                        optimization_method: models.OptimizationMethod,
                                        optimization_cub_index: models.OptimizationCubIndex,
                                        run_summary: RunSummary,
                                        max_iterations: int = config["ORF"]["ZSCORE_MAX_ITERATIONS"]):

    with Timer() as timer:
        initial_sequence_zscore = _calculate_zscore_for_sequence(
            sequence=sequence,
            user_input=user_input,
            optimization_cub_index=optimization_cub_index,
        )
        initial_sequence_score = None
        aa_to_selected_codon = {}
        iterations_count = 0
        iterations_summary = []
        for run in range(max_iterations):
            aa_to_selected_codon = {}
            iterations_count = run + 1
            codons_to_zscore = {}
            for codon in nt_to_aa.keys():
                candidate_codon_sequence, candidate_codon_count = _change_all_codons_of_aa(sequence, codon)
                codons_to_zscore[codon] = _calculate_zscore_for_sequence(
                    sequence=candidate_codon_sequence,
                    user_input=user_input,
                    optimization_cub_index=optimization_cub_index,
                )
            if optimization_method.is_zscore_ratio_score_optimization:
                min_zscore = min(score.min_zscore for score in codons_to_zscore.values())
                max_zscore = max(score.max_zscore for score in codons_to_zscore.values())

                for zscore in codons_to_zscore.values():
                    zscore.normalize(min_zscore=min_zscore, max_zscore=max_zscore)
                if initial_sequence_score is None:
                    initial_min_zscore = min(min_zscore, initial_sequence_zscore.min_zscore)
                    initial_max_zscore = max(max_zscore, initial_sequence_zscore.max_zscore)
                    initial_sequence_zscore.normalize(min_zscore=initial_min_zscore, max_zscore=initial_max_zscore)

            codons_to_total_score = {
                codon: get_total_score(zscore=zscore,
                                       optimization_method=optimization_method,
                                       tuning_parameter=user_input.tuning_parameter) for
                codon, zscore in codons_to_zscore.items()
            }

            if initial_sequence_score is None:
                initial_sequence_score = get_total_score(zscore=initial_sequence_zscore,
                                                         optimization_method=optimization_method,
                                                         tuning_parameter=user_input.tuning_parameter)
                previous_sequence_score = initial_sequence_score

            for aa in synonymous_codons.keys():
                selected_codon = _find_best_synonymous_codon_for_aa(codons_list=synonymous_codons[aa],
                                                                    codons_to_score=codons_to_total_score)
                aa_to_selected_codon[aa] = selected_codon

            # create new sequence by replacing all synonymous codons
            new_sequence = sequence
            for aa in aa_to_selected_codon:
                new_sequence, _ = _change_all_codons_of_aa(new_sequence, aa_to_selected_codon[aa])

            # Calculate score after all replacements
            zscore = _calculate_zscore_for_sequence(
                sequence=new_sequence,
                user_input=user_input,
                optimization_cub_index=optimization_cub_index,
            )

            if optimization_method.is_zscore_ratio_score_optimization:
                updated_min_zscore = min(min_zscore, zscore.min_zscore)
                updated_max_zscore = max(max_zscore, zscore.max_zscore)
                zscore.normalize(min_zscore=updated_min_zscore, max_zscore=updated_max_zscore)

            score = get_total_score(zscore=zscore,
                                    optimization_method=optimization_method,
                                    tuning_parameter=user_input.tuning_parameter)

            iteration_summary = {
                "aa_to_selected_codon": aa_to_selected_codon,
                "sequence_score": score,
            }
            iterations_summary.append(iteration_summary)

            if new_sequence == sequence or previous_sequence_score > score:
                break
            else:
                sequence = new_sequence
                previous_sequence_score = score

    orf_summary = {
        "iterations_count": iterations_count,
        "aa_to_optimal_codon": aa_to_selected_codon,
        "initial_sequence_optimization_score": initial_sequence_score,
        "final_sequence_optimization_score": score,
        "run_time": timer.elapsed_time,
        "iterations_summary": iterations_summary,
    }

    run_summary.add_to_run_summary("orf", orf_summary)

    return sequence


# --------------------------------------------------------------
def _find_best_synonymous_codon_for_aa(codons_to_score: typing.Dict[str, float],
                                       codons_list: typing.Sequence[str]) -> str:
    aa_codons_to_score = {codon: score for codon, score in codons_to_score.items() if codon in codons_list}
    selected_aa_codon = max(aa_codons_to_score, key=aa_codons_to_score.get)
    return selected_aa_codon


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
        initial_wanted_hosts_scores=wanted_hosts_scores,
        wanted_hosts_weights=wanted_hosts_weights,
        initial_unwanted_hosts_scores=unwanted_hosts_scores,
        unwanted_hosts_weights=unwanted_hosts_weights,
    )


# --------------------------------------------------------------
def get_total_score(zscore: models.SequenceZscores,
                    optimization_method: models.OptimizationMethod,
                    tuning_parameter: float) -> float:
    if optimization_method.is_zscore_diff_score_optimization:
        return _calculate_zscore_diff_score(zscore=zscore,
                                            tuning_parameter=tuning_parameter)

    if optimization_method.is_zscore_ratio_score_optimization:
        return _calculate_zscore_ratio_score(zscore=zscore,
                                             tuning_parameter=tuning_parameter)

    if optimization_method.is_zscore_weakest_link_score_optimization:
        return _calculate_zscore_weakest_link_score(zscore=zscore,
                                                    tuning_parameter=tuning_parameter,
                                                    )


# --------------------------------------------------------------

# def _calculate_zscore_for_sequence_old(sequence: str,
#                                        user_input: models.UserInput,
#                                        optimization_method: models.OptimizationMethod,
#                                        optimization_cub_index: models.OptimizationCubIndex):
#     optimization_cub_index_value = optimization_cub_index.value.lower()
#
#     std_key = F"{optimization_cub_index_value}_std"
#     average_key = F"{optimization_cub_index_value}_avg"
#     weights = F"{optimization_cub_index_value}_profile"
#
#     optimized_organisms_scores = []
#     optimized_organisms_weights = []
#     deoptimized_organisms_scores = []
#     deoptimized_organisms_weights = []
#
#     for organism in user_input.organisms:
#         sigma = getattr(organism, std_key)
#         miu = getattr(organism, average_key)
#         profile = getattr(organism, weights)
#         index_score = general_geomean([sequence], weights=profile)[0]
#         organism_score = (index_score - miu) / sigma
#         # logger.info(F"CUB score for organism {organism.name} is: {index_score}")
#         if organism.is_optimized:
#             optimized_organisms_scores.append(organism_score)
#             optimized_organisms_weights.append(organism.optimization_priority)
#         else:
#             deoptimized_organisms_scores.append(organism_score)
#             deoptimized_organisms_weights.append(organism.optimization_priority)
#
#     alpha = user_input.tuning_parameter
#     if optimization_method.is_zscore_diff_score_optimization:
#         return _calculate_zscore_diff_score(optimized_organisms_scores=optimized_organisms_scores,
#                                             deoptimized_organisms_scores=deoptimized_organisms_scores,
#                                             optimized_organisms_weights=optimized_organisms_weights,
#                                             deoptimized_organisms_weights=deoptimized_organisms_weights,
#                                             tuning_parameter=alpha)
#
#     if optimization_method.is_zscore_ratio_score_optimization:
#         return _calculate_zscore_ratio_score(optimized_organisms_scores=optimized_organisms_scores,
#                                              deoptimized_organisms_scores=deoptimized_organisms_scores,
#                                              optimized_organisms_weights=optimized_organisms_weights,
#                                              deoptimized_organisms_weights=deoptimized_organisms_weights,
#                                              tuning_parameter=alpha)
#
#     if optimization_method.is_zscore_weakest_link_score_optimization:
#         return _calculate_zscore_weakest_link_score(
#             optimized_organisms_scores=optimized_organisms_scores,
#             deoptimized_organisms_scores=deoptimized_organisms_scores,
#             optimized_organisms_weights=optimized_organisms_weights,
#             deoptimized_organisms_weights=deoptimized_organisms_weights,
#             tuning_parameter=alpha,
#         )
#
#
#     raise NotImplementedError(F"Optimization method: {optimization_method}")
#

# --------------------------------------------------------------
def _calculate_zscore_diff_score(zscore: models.SequenceZscores,
                                 tuning_parameter: float) -> float:
    mean_opt_index = average(zscore.wanted_hosts_scores, weights=zscore.wanted_hosts_weights)
    mean_deopt_index = average(zscore.unwanted_hosts_scores, weights=zscore.unwanted_hosts_weights)
    return tuning_parameter * mean_opt_index - (1 - tuning_parameter) * mean_deopt_index


# --------------------------------------------------------------
def _calculate_zscore_ratio_score(zscore: models.SequenceZscores,
                                  tuning_parameter: float) -> float:
    mean_opt_index = average(zscore.wanted_hosts_scores, weights=zscore.wanted_hosts_weights)
    mean_deopt_index = average(zscore.unwanted_hosts_scores, weights=zscore.unwanted_hosts_weights)
    epsilon = 10 ** -9
    mean_deopt_index = mean_deopt_index if mean_deopt_index != 0 else epsilon

    return (mean_opt_index ** tuning_parameter) / (mean_deopt_index ** (1 - tuning_parameter))


# --------------------------------------------------------------
def _calculate_zscore_weakest_link_score(zscore: models.SequenceZscores,
                                         tuning_parameter: float) -> float:
    weighted_optimized_organisms_scores = [zscore.wanted_hosts_scores[i] * zscore.wanted_hosts_weights[i] for i in
                                           range(len(zscore.wanted_hosts_scores))]
    weighted_deoptimized_organisms_scores = [
        zscore.unwanted_hosts_scores[i] * zscore.unwanted_hosts_weights[i] for
        i in range(len( zscore.unwanted_hosts_scores))
    ]

    return (tuning_parameter * min(weighted_optimized_organisms_scores) -
            (1 - tuning_parameter) * max(weighted_deoptimized_organisms_scores))
