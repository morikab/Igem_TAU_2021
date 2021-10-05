from modules import Zscore_calculation
from modules.shared_functions_and_vars import nt_to_aa




def change_all_codons_of_aa(seq, selected_codon, nt_to_aa = nt_to_aa):
    '''
    change all synonymous codons in the str seq to the selected codon, return str as well
    '''
    split_seq = [seq[i:i+3].upper() for i in range(0, len(seq), 3)]
    print(split_seq)
    new_split_seq = []
    for codon in split_seq:
        if nt_to_aa[codon]==nt_to_aa[selected_codon]:
            new_split_seq.append(selected_codon)
        else:
            new_split_seq.append(codon)
    return ''.join(new_split_seq)


# in each round- check all single synonymous codon changes and calculate the Zscore- take the best one
def hill_climbing_optimize_by_zscore(seq, inp_dict, opt_type, max_iter, nt_to_aa = nt_to_aa):
    '''
    hill climbing function for performing codon optimization
    in each iteration- for each codon, change all synonymous codons to a specific one and test the zscore of the new sequence
    after each iteration, select the sequence with the best zscore- if it was not changed since the last iteration, break
    the maximum number of iterations allowed is "max_iter"
    @seq: str, tested seq
    @inp_dict: input dict after usr_inp code
    @opt_type: 'cai' or 'tai'
    @max_iter: maximal number of iterations to perform
    return: seq, which is the optimized sequence
    '''
    seq_options = {}
    mean_opt_index, mean_deopt_index, zscore = Zscore_calculation.ZscoreModule.run_module(seq, inp_dict, opt_type)
    seq_options[seq] = zscore
    for run in range(max_iter):
        for codon in nt_to_aa.keys():
            tested_seq = change_all_codons_of_aa(seq, codon)

            mean_opt_index, mean_deopt_index, zscore = \
                Zscore_calculation.ZscoreModule.run_module(tested_seq, inp_dict, opt_type)

            seq_options[codon] = zscore

        changed_codon = max(seq_options, key=seq_options.get)
        print(changed_codon, nt_to_aa[changed_codon], seq_options[changed_codon])
        new_seq = change_all_codons_of_aa(seq, changed_codon)
        if new_seq == seq:
            break
        else:
            seq = new_seq
    return seq





