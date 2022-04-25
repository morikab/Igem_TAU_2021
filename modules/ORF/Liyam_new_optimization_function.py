from modules import models
from modules.stats.optimization import OptimizationModule
from modules.shared_functions_and_vars import nt_to_aa


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
    for run in range(max_iter):
        # TODO - we change in each iteration a single codon. We may consider changing at most X codons at a time to
        #  reduce risk of falling to a local maxima.
        for codon in nt_to_aa.keys():
            tested_seq = change_all_codons_of_aa(seq, codon)

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
