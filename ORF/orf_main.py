from src.optimization import optimize_sequence
from src.organism import Organism, Gene

def orf_main(target_gene_seq, input_dict):
    """
    :param target_gene_seq: target gene aa sequence (str)
    :param input_dict: input from GUI (dict). Format:
    the keys are the organism name and the values are a dict of all other params:
    {'ecoli': {genome_path: '.gb' - mandatory
                genes_HE: '.csv'
                features_to_generate: 'CAI'
                optimized: True - mandatory
                }
     'bacillus': {genome_path: '.gb'
                genes_HE: '.csv'
                features_to_generate: 'CAI'
                optimized: False
                }
     }
    :return: optimized sequence (Seq)
    """

    target_gene = Gene(target_gene_seq)

    high_expression_organisms = [Organism(name=key, genome_path=value['genome_path'], mrna_levels_path=value['genes_HE'])
                                 for key, value in input_dict.items() if value['optimized']]

    low_expression_organisms = [Organism(name=key, genome_path=value['genome_path'], mrna_levels_path=value['genes_HE'])
                                 for key, value in input_dict.items() if not value['optimized']]

    optimized_sequence = optimize_sequence(target_gene=target_gene,
                      high_expression_organisms=high_expression_organisms, low_expression_organisms=low_expression_organisms)

    return optimized_sequence