import os
import typing
from user_IO.input_functions import *

from logger_factory import LoggerFactory

# initialize the logger object
logger = LoggerFactory.create_logger("user_input")


class UserInputModule(object):
    @staticmethod
    def get_name() -> str:
        return "User Input"

    @classmethod
    def run_module(cls, user_inp_raw, input: typing.Optional[typing.Dict] = None) -> typing.Dict:
        logger.info('##########################')
        logger.info('# USER INPUT INFORMATION #')
        logger.info('##########################')
        #
        # # TODO - use input var from yarin
        return cls._parse_input(user_inp_raw)

    @classmethod
    def _parse_input(cls, usr_inp):
        '''
        :param usr_inp: in the following format
            {
        'ecoli': {'genome_path': '.gb' - mandatory
                    'genes_HE': '.csv'
                    'optimized': True - mandatory
                    }
         'bacillus': {'genome_path': '.gb'
                    'genes_HE': '.csv'
                    'optimized': False
                    }
            }
        :return: an extended dictionary with the following items:
        @selected_prom : final used list of promoters for MAST
        @sequence : the ORF to optimize
        and for each organism- the key is the scientific organism name, and the value is _parse_single_input for the organism's
        input
        '''
        full_inp_dict = {}
        for key, val in usr_inp['organisms'].items():
            if key in ['sequence', 'selected_promoters']:
                continue
            org_name, org_dict = cls._parse_single_input(val)
            full_inp_dict[
                org_name] = org_dict  # creating the sub dictionary for each organism- where the key is the scientific name and the value is the following dict:

        # add non org specific keys to dict
        orf_fasta_fid = usr_inp['sequence']
        orf_seq = str(SeqIO.read(orf_fasta_fid, 'fasta').seq)
        prom_fasta_fid = usr_inp['selected_promoters']
        selected_prom = {}
        logger.info(f'\n\nSequence to be optimized given in the following file {orf_fasta_fid}')
        logger.info(f'containing this sequence: {orf_seq}')
        if prom_fasta_fid is not None:  #
            selected_prom = fasta_to_dict(prom_fasta_fid)
            logger.info(f'Promoter options ranked are given in the following file {prom_fasta_fid}, '
                        f'which contains {len(selected_prom)} promoters')
        else:
            logger.info(
                f'External promoter options were not supplied. endogenous promoters will be used for optimization.'
                f'promoters from the 1/3 most highly expressed genes of all organisms are used- ')
            for org, org_dict in full_inp_dict.items():
                if org_dict['optimized']:
                    org_third_he_prom_dict = org_dict['third_most_HE']
                    for prom_name, prom_Seq in org_third_he_prom_dict.items():
                        selected_prom[prom_name + ' from organism: ' + org] = prom_Seq
                    logger.info(f'{len(org_third_he_prom_dict)} promoters are selected from {org}')
            logger.info(f'Resulting in a total of {len(selected_prom)} used for promoter selection and optimization')

        full_inp_dict['selected_prom'] = selected_prom
        full_inp_dict['sequence'] = orf_seq
        return full_inp_dict

    @staticmethod
    def _parse_single_input(val):
        '''
        create the relevant information for each organism
        :param val: the dictionary supplied for every organism:
        {'genome_path': '.gb'
                    'genes_HE': '.csv'
                    'optimized': False }
        :return:
        @org_name = scientific organism name
        @org_dict = {
            'tgcn': {anti codon:number of occurences}, for ORF model
            '200bp_promoters': {gene name and function: prom}, promoter model
            'third_most_HE': same as 200bp_promoters but only third most highly expressed promoters
            'gene_cds':  {gene name and function : cds}, for ORF model
            'intergenic':{position along the genome: intergenic sequence}, promoter model
            'expression_estimation_of_all_genes': {'gene name and function: estimated expression }when the expression csv is not given- the CAI is used as expression levels
            'CAI_score_of_all_genes': {'gene_name': cai} ORF and promoter
            'cai_profile': {codon:cai score},  # ORF model
            'optimized': bool- True if organism is optimized}
        '''
        gb_path = val['genome_path']
        gb_file = SeqIO.read(gb_path, format='gb')
        org_name = find_org_name(gb_file)
        logger.info(f'\nInformation about {org_name}:')
        if val['optimized']:
            logger.info('Organism is optimized')
        else:
            logger.info('Organism is deoptimized')
        tgcn_dict = find_tgcn(gb_path)
        logger.info(f'Number of tRNA genes found: {sum(tgcn_dict.values())}, for {len(tgcn_dict)} anticodons out of 61')
        prom200_dict, prom400_dict, cds_dict, intergenic_dict, cai_dict, cai_weights = extract_gene_data(gb_path)
        logger.info(f'Number of genes: {len(cds_dict)}, number of intregenic regions: {len(intergenic_dict)}')
        highly_expressed_gene_name_list = extract_highly_expressed_gene_names(cai_dict)

        org_dict = {
            'tgcn': tgcn_dict,  # tgcn dict {codon:number of occurences} for ORF model
            '200bp_promoters': prom200_dict,  # prom_dict {gene name and function: prom}, promoter model
            'third_most_HE': {name: seq for name, seq in prom200_dict.items()
                              if name in highly_expressed_gene_name_list},
            # '400bp_promoters': prom400_dict,  # prom_dict {gene name and function: prom}, promoter model
            'gene_cds': cds_dict,  # cds dict {gene name and function : cds}, for ORF model
            'intergenic': intergenic_dict,
            # intergenic dict {position along the genome: intergenic sequence}, promoter model
            'expression_estimation_of_all_genes': cai_dict,
            # when the expression csv is not given- the CAI is used as expression levels
            'CAI_score_of_all_genes': cai_dict,  # {'gene_name': expression} ORF and promoter
            'cai_profile': cai_weights,  # ORF model
            'optimized': val['optimized']}
        return org_name, org_dict
