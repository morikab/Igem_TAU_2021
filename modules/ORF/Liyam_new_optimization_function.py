from modules import Zscore_calculation
from modules.shared_functions_and_vars import nt_to_aa


def change_all_codons_of_aa(seq, selected_codon):
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


# in each round- check all single synonymous codon changes and calculate the Zscore - take the best one
def hill_climbing_optimize_by_zscore(seq, inp_dict, cai_or_tai, max_iter, optimization_type):
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
    mean_opt_index, mean_deopt_index, zscore, weakest_score = Zscore_calculation.ZscoreModule.run_module(seq,
                                                                                                         inp_dict,
                                                                                                         cai_or_tai)
    if 'average' in cai_or_tai:
        seq_options[seq] = zscore
    else:
        seq_options[seq] = weakest_score
    for run in range(max_iter):
        # TODO - we change in each iteration a single codon. We may consider changing at most X codoons at a time to
        #  reduce risk of falling to a local maxima.
        for codon in nt_to_aa.keys():
            tested_seq = change_all_codons_of_aa(seq, codon)

            mean_opt_index, mean_deopt_index, zscore, weakest_score = \
                Zscore_calculation.ZscoreModule.run_module(tested_seq, inp_dict, cai_or_tai)
            if 'average' in optimization_type:
                seq_options[tested_seq] = zscore
            else:
                seq_options[tested_seq] = weakest_score

        new_seq = max(seq_options, key=seq_options.get)
        if new_seq == seq:
            break
        else:
            seq = new_seq
    return seq
