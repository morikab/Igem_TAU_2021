from Igem_TAU_2021.ORF.optimization import optimize_sequence
from Igem_TAU_2021.ORF.organism import Organism, Gene
<<<<<<< Updated upstream
=======

#todo: add a statistical analysis of how close the organisms are- like what is the best codon for eah AA
#and are they close
>>>>>>> Stashed changes

def orf_main(full_input_dict):
    """
    :param full_input_dict: input from GUI parser (dict). Format:
    full_inp_dict[org_name] = {
            'tgcn': tgcn_dict,  # tgcn dict {codon:number of occurences}
            '200bp_promoters': prom200_dict,  # prom_dict {gene name and function: prom}
            '400bp_promoters': prom400_dict,  # prom_dict {gene name and function: prom}
            'gene_cds': cds_dict,  # cds dict {gene name and function : cds}
            'intergenic': intergenic_dict,  # intergenic dict {position along the genome: intergenic sequence}
            'caiScore_dict': cai_dict,
            'cai_profile': cai_weights,
            'optimized': val['optimized']}  # is the sequence in the optimized or deoptimized group- bool

            @selected_prom : final used list of promoters for MAST
            @sequence : the ORF to optimize

    :return: optimized sequence (Biopython Seq)
    """

    target_gene = Gene(full_input_dict['sequence'])

    input_organisms = full_input_dict
    input_organisms.pop('sequence')
    input_organisms.pop('selected_prom')

    high_expression_organisms = [Organism(name=org_name, gene_cds=dict['gene_cds'], tgcn=dict['tgcn'], cai_weights=dict['cai_profile'])
                                 for org_name, dict in input_organisms.items() if dict['optimized']]

    low_expression_organisms = [Organism(name=org_name, gene_cds=dict['gene_cds'], tgcn=dict['tgcn'], cai_weights=dict['cai_profile'])
                                 for org_name, dict in input_organisms.items() if not dict['optimized']]

    optimized_sequence = optimize_sequence(target_gene=target_gene,
                      high_expression_organisms=high_expression_organisms, low_expression_organisms=low_expression_organisms)

    return optimized_sequence