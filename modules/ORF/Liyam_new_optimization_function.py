import typing
from logger_factory.logger_factory import LoggerFactory
from modules import models
from modules.stats.optimization import OptimizationModule
from modules.shared_functions_and_vars import nt_to_aa
from modules.shared_functions_and_vars import synonymous_codons


logger = LoggerFactory.get_logger()


def change_all_codons_of_aa(seq: str, selected_codon: str):
    '''
    change all synonymous codons in the str seq to the selected codon, return str as well
    '''
    split_seq = [seq[i:i+3].upper() for i in range(0, len(seq), 3)]
    new_split_seq = []
    for codon in split_seq:
        if nt_to_aa[codon] == nt_to_aa[selected_codon]:
            new_split_seq.append(selected_codon)
        else:
            new_split_seq.append(codon)
    return ''.join(new_split_seq)


# in each round - check all single synonymous codon changes and calculate optimization score - take the best one
def hill_climbing_optimize_by_zscore(seq: str,
                                     user_input: models.UserInput,
                                     cai_or_tai: str,
                                     max_iter: int,
                                     optimization_method: models.OptimizationMethod):
    """
    hill climbing function for performing codon optimization
    in each iteration - for each codon, change all synonymous codons to a specific one and test the zscore of the new
    sequence after each iteration, select the sequence with the best zscore - if it was not changed since the last
    iteration, break. The maximum number of iterations allowed is "max_iter"
    @seq: str, tested seq
    @inp_dict: input dict after usr_inp code
    @opt_type: 'cai' or 'tai'
    @max_iter: maximal number of iterations to perform
    return: seq, which is the optimized sequence
    """

    seq_options = {}
    score = OptimizationModule.run_module(
        final_seq=seq,
        user_input=user_input,
        cai_or_tai=cai_or_tai,
        optimization_method=optimization_method,
    )
    seq_options[seq] = score

    # Single codon replacement
    for run in range(max_iter):
        # TODO - we change in each iteration a single codon. We may consider changing at most X codons at a time to
        #  reduce risk of falling to a local maxima.
        tested_seq_to_codon = {}
        for codon in nt_to_aa.keys():
            tested_seq = change_all_codons_of_aa(seq, codon)
            tested_seq_to_codon[tested_seq] = codon

            score = OptimizationModule.run_module(
                final_seq=tested_seq,
                user_input=user_input,
                cai_or_tai=cai_or_tai,
                optimization_method=optimization_method,
            )
            seq_options[tested_seq] = score

        new_seq = max(seq_options, key=seq_options.get)
        if new_seq == seq:
            break
        else:
            seq = new_seq
    return seq


# in each round - check all single synonymous codon changes and calculate optimization score - take the best one
def hill_climbing_optimize_aa_bulk_by_zscore(seq: str,
                                             user_input: models.UserInput,
                                             cai_or_tai: str,
                                             max_iter: int,
                                             optimization_method: models.OptimizationMethod):
    """
    hill climbing function for performing codon optimization
    in each iteration - for each codon, change all synonymous codons to a specific one and test the zscore of the new
    sequence after each iteration, select the sequence with the best zscore - if it was not changed since the last
    iteration, break. The maximum number of iterations allowed is "max_iter"
    @seq: str, tested seq
    @inp_dict: input dict after usr_inp code
    @opt_type: 'cai' or 'tai'
    @max_iter: maximal number of iterations to perform
    return: seq, which is the optimized sequence
    """
    optimization_method = models.OptimizationMethod.hill_climbing_average
    seq_options = {}
    score = OptimizationModule.run_module(
        final_seq=seq,
        user_input=user_input,
        cai_or_tai=cai_or_tai,
        optimization_method=optimization_method,
    )
    seq_options[seq] = score

    for run in range(max_iter):
        def find_best_aa_synonymous_codon(codons_list, seq_to_change: str) -> str:
            aa_seq_options = {}
            for aa_codon in codons_list:
                option_seq = change_all_codons_of_aa(seq_to_change, aa_codon)
                aa_seq_options[aa_codon] = OptimizationModule.run_module(
                    final_seq=option_seq,
                    user_input=user_input,
                    cai_or_tai=cai_or_tai,
                    optimization_method=optimization_method,
                )
            aa_new_seq = max(aa_seq_options, key=aa_seq_options.get)
            return aa_new_seq

        aa_to_selected_codon = {}
        # Find best synonymous_codon per aa
        for aa in synonymous_codons.keys():
            selected_aa_codon = find_best_aa_synonymous_codon(codons_list=synonymous_codons[aa], seq_to_change=seq)
            aa_to_selected_codon[aa] = selected_aa_codon

        # create new seq by replacing all synonymous codons
        new_seq = seq
        for aa in aa_to_selected_codon:
            new_seq = change_all_codons_of_aa(new_seq, aa_to_selected_codon[aa])

        # Calculate score after all replacements
        score = OptimizationModule.run_module(
            final_seq=new_seq,
            user_input=user_input,
            cai_or_tai=cai_or_tai,
            optimization_method=optimization_method,
        )
        seq_options[new_seq] = score
        logger.info(F"New seq for iteration {run} with score of: {score}")

        if new_seq == seq:
            break
        else:
            seq = new_seq

    return seq
