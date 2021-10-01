from modules.ORF.optimization import optimize_sequence
from modules.ORF.organism import Organism


#todo: add a statistical analysis of how close the organisms are- like what is the best codon for eah AA
#and are they close


class ORFModule(object):

    @staticmethod
    def run_module(full_input_dict, feature, local_maximum):
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

        target_gene = full_input_dict['sequence']
        input_organisms = full_input_dict['organisms']

        high_expression_organisms = [
            Organism(name=org_name, tai_weights=org_dict['tai_profile'], cai_weights=org_dict['cai_profile'],
                     feature_to_generate=feature, cai_std=org_dict['cai_std'],  tai_std=org_dict['tai_std'])
            for org_name, org_dict in input_organisms.items() if org_dict['optimized']]

        low_expression_organisms = [
            Organism(name=org_name, tai_weights=org_dict['tai_profile'], cai_weights=org_dict['cai_profile'],
                     feature_to_generate=feature, cai_std=org_dict['cai_std'],  tai_std=org_dict['tai_std'])
            for org_name, org_dict in input_organisms.items() if not org_dict['optimized']]

        optimized_sequence = optimize_sequence(target_gene=target_gene,
                                               high_expression_organisms=high_expression_organisms,
                                               low_expression_organisms=low_expression_organisms,
                                               tuning_param=full_input_dict['tuning_param'],
                                               local_maximum = local_maximum,
                                               )

        return optimized_sequence
