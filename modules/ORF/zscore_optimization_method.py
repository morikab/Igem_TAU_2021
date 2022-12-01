import typing

from collections import defaultdict
from logger_factory.logger_factory import LoggerFactory
from modules import models
from modules.configuration import Configuration
from modules.stats.optimization import OptimizationModule
from modules.shared_functions_and_vars import nt_to_aa
from modules.shared_functions_and_vars import synonymous_codons


logger = LoggerFactory.get_logger()
config = Configuration.get_config()


def change_all_codons_of_aa(seq: str, selected_codon: str) -> typing.Tuple[str, int]:
    '''
    change all synonymous codons in the str seq to the selected codon, return str as well
    '''
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


# In each round - check all single synonymous codon changes and calculate optimization score - take the best one
def hill_climbing_optimize_by_zscore(seq: str,
                                     user_input: models.UserInput,
                                     optimization_cub_score: models.OptimizationCubScore,
                                     optimization_method: models.OptimizationMethod,
                                     max_iter: int = config["ORF"]["HILL_CLIMBING_MAX_ITERATIONS"]):
    """
    hill climbing function for performing codon optimization
    in each iteration - for each codon, change all synonymous codons to a specific one and test the zscore of the new
    sequence after each iteration, select the sequence with the best zscore - if it was not changed since the last
    iteration, break. The maximum number of iterations allowed is "max_iter"
    return: optimized sequence
    """
    seq_options = {}
    score = OptimizationModule.run_module(
        final_seq=seq,
        user_input=user_input,
        optimization_cub_score=optimization_cub_score,
        optimization_method=optimization_method,
    )
    seq_options[seq] = score

    # Single codon replacement
    for run in range(max_iter):
        # TODO - we change in each iteration a single codon. We may consider changing at most X codons at a time to
        #  reduce risk of falling to a local maxima.
        tested_seq_to_codon = {}
        for codon in nt_to_aa.keys():
            tested_seq, _ = change_all_codons_of_aa(seq, codon)
            tested_seq_to_codon[tested_seq] = codon

            score = OptimizationModule.run_module(
                final_seq=tested_seq,
                user_input=user_input,
                optimization_cub_score=optimization_cub_score,
                optimization_method=optimization_method,
            )
            seq_options[tested_seq] = score

        new_seq = max(seq_options, key=seq_options.get)
        logger.info(F"selected_codon in iteration {run} is: {tested_seq_to_codon[new_seq]}")

        if new_seq == seq:
            break
        else:
            seq = new_seq
    return seq


#
def hill_climbing_optimize_aa_bulk_by_zscore(seq: str,
                                             user_input: models.UserInput,
                                             optimization_method: models.OptimizationMethod,
                                             optimization_cub_score: models.OptimizationCubScore,
                                             max_iter: int = config["ORF"]["HILL_CLIMBING_MAX_ITERATIONS"]):
    """
    Hill climbing function for performing codon optimization
    in each iteration - in each round, check all single synonymous codon changes, calculate optimization score and
    take the best one
    """
    optimization_method = models.OptimizationMethod.hill_climbing_average
    seq_options = {}
    original_seq = seq

    # Get codon split in the original sequence
    codon_counts_in_seq = defaultdict(int)
    split_seq = [original_seq[i:i + 3].upper() for i in range(0, len(original_seq), 3)]
    for codon in split_seq:
        codon_counts_in_seq[codon] += 1

    logger.info(F"codon_counts_in_seq: {codon_counts_in_seq}")

    score = OptimizationModule.run_module(
        final_seq=seq,
        user_input=user_input,
        optimization_cub_score=optimization_cub_score,
        optimization_method=optimization_method,
    )
    seq_options[seq] = score
    for run in range(max_iter):
        def find_best_aa_synonymous_codon(codons_list, seq_to_change: str) -> str:
            aa_seq_options = {}
            for aa_codon in codons_list:
                option_seq, option_codon_count = change_all_codons_of_aa(seq_to_change, aa_codon)
                # logger.info(F"Number of occurrences of codon {aa_codon} in sequence is {option_codon_count}")
                logger.info(F"Running for codon: {aa_codon}")
                aa_seq_options[aa_codon] = OptimizationModule.run_module(
                    final_seq=option_seq,
                    user_input=user_input,
                    optimization_cub_score=optimization_cub_score,
                    optimization_method=optimization_method,
                )
                # logger.info(F"z-score after changing codon {aa_codon} is: {aa_seq_options[aa_codon]}")
            aa_new_seq = max(aa_seq_options, key=aa_seq_options.get)
            return aa_new_seq
        logger.info(F"zscore of sequence in run {run} is: {seq_options[seq]}")
        aa_to_selected_codon = {}
        # Find the best synonymous_codon per aa
        for aa in synonymous_codons.keys():
            selected_aa_codon = find_best_aa_synonymous_codon(codons_list=synonymous_codons[aa], seq_to_change=seq)
            aa_to_selected_codon[aa] = selected_aa_codon

        logger.info(F"aa_to_selected_codon in iteration {run} is: {aa_to_selected_codon}")

        # create new seq by replacing all synonymous codons
        new_seq = seq
        for aa in aa_to_selected_codon:
            new_seq, _ = change_all_codons_of_aa(new_seq, aa_to_selected_codon[aa])

        # Calculate score after all replacements
        score = OptimizationModule.run_module(
            final_seq=new_seq,
            user_input=user_input,
            optimization_cub_score=optimization_cub_score,
            optimization_method=optimization_method,
        )
        seq_options[new_seq] = score
        logger.info(F"New seq for iteration {run} with score of: {score}")

        if new_seq == seq:
            break
        else:
            seq = new_seq

    logger.info(F"Original score is: {seq_options[original_seq]}")
    return seq
