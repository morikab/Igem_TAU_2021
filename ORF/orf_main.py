from Igem_TAU_2021.ORF.optimization import optimize_sequence
from Igem_TAU_2021.ORF.organism import Organism, Gene

def orf_main(target_gene_seq_record, full_input_dict):
    """
    :param target_gene_seq: target gene aa sequence (Seq record)
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

    :return: optimized sequence (Biopython Seq)
    """

    target_gene = Gene(target_gene_seq_record)

    high_expression_organisms = [Organism(name=org_name, gene_cds=dict['gene_cds'], tgcn=dict['tgcn'], cai_weights=dict['cai_profile'])
                                 for org_name, dict in full_input_dict.items() if dict['optimized']]

    low_expression_organisms = [Organism(name=org_name, gene_cds=dict['gene_cds'], tgcn=dict['tgcn'], cai_weights=dict['cai_profile'])
                                 for org_name, dict in full_input_dict.items() if not dict['optimized']]

    optimized_sequence = optimize_sequence(target_gene=target_gene,
                      high_expression_organisms=high_expression_organisms, low_expression_organisms=low_expression_organisms)

    return optimized_sequence