import typing

from Bio import SeqIO

from modules import models
from modules.user_IO.input_functions import *
from modules.ORF.TAI import TAI
from modules.ORF.calculating_cai import general_geomean
from modules.logger_factory import LoggerFactory

# initialize the logger object
logger = LoggerFactory.create_logger("user_input")


class UserInputModule(object):
    @staticmethod
    def get_name() -> str:
        return "User Input"

    @classmethod
    def run_module(cls, user_inp_raw: typing.Dict) -> models.UserInput:
        logger.info('##########################')
        logger.info('# USER INPUT INFORMATION #')
        logger.info('##########################')
        return cls._parse_input(user_inp_raw)

    @classmethod
    def _parse_input(cls, user_input):
        """
        :param user_input: in the following format
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
            @organism: contains the org names as keys, and for each organism- the key is the scientific organism name,
            and the value is _parse_single_input for the organism's input
        """

        full_inp_dict = {}

        organisms_list = []
        organisms_names = set()

        for key, val in user_input['organisms'].items():
            try:
                organism = cls._parse_single_input(val)
            except:
                raise ValueError(f'Error in input: {key}, re-check your input')
            if organism.name in organisms_names:
                raise ValueError(f"Organism: {organism.name}'s genome is inserted twice, re-check your input.")
            organisms_list.append(organism)
            organisms_names.add(organism.name)

        # Read ORF sequence
        orf_fasta_fid = user_input['sequence']
        orf_seq = None
        if orf_fasta_fid is not None:
            try:
                orf_seq = str(SeqIO.read(orf_fasta_fid, 'fasta').seq)
            except:
                raise ValueError(
                    f'Error in protein .fasta file: {orf_fasta_fid}, make sure you inserted an undamaged .fasta file containing a single recored')
        logger.info(f'\n\nSequence to be optimized given in the following file {orf_fasta_fid}')
        logger.info(f'containing this sequence: {orf_seq}')

        tuning_parameter = user_input['tuning_param']

        return models.UserInput(organisms=organisms_list,
                                sequence=orf_seq,
                                tuning_parameter=tuning_parameter)

    @staticmethod
    def _parse_single_input(organism_input):
        """
        create the relevant information for each organism
        :param organism_input: the dictionary supplied for every organism:
            {'genome_path': '.gb', 'genes_HE': '.csv', 'optimized': False }
        :return:
        @organism_object: Parsed organism object
        """
        gb_path = organism_input['genome_path']
        exp_csv_fid = organism_input['expression_csv']
        try:
            gb_file = SeqIO.read(gb_path, format='gb')
        except:
            raise ValueError(
                f'Error in genome GenBank file: {gb_path}, make sure you inserted an undamaged .gb file containing '
                f'the full genome sequence and annotations'
            )

        organism_name = find_org_name(gb_file)
        logger.info(f'\nInformation about {organism_name}:')
        is_optimized = organism_input['optimized']
        if is_optimized:
            logger.info('Organism is optimized')
        else:
            logger.info('Organism is deoptimized')

        cds_dict, estimated_expression = extract_gene_data(gb_path, exp_csv_fid)
        logger.info(f'Number of genes: {len(cds_dict)}')
        gene_names = list(cds_dict.keys())
        cai_weights = calculate_cai_weights_for_input(cds_dict, estimated_expression, exp_csv_fid)
        cai_scores = general_geomean(sequence_lst=cds_dict.values(), weights=cai_weights)
        cai_scores_dict = {gene_names[i]: cai_scores[i] for i in range(len(gene_names))}

        tai_weights = None
        tai_scores_dict = {}
        try:
            tai_weights = TAI(tai_from_tgcnDB(organism_name)).index
            tai_scores = general_geomean(sequence_lst=cds_dict.values(), weights=tai_weights)
            tai_scores_dict = {gene_names[i]: tai_scores[i] for i in range(len(gene_names))}
        except:
            pass

        organism_object = models.Organism(name=organism_name,
                                          cai_profile=cai_weights,
                                          tai_profile=tai_weights,
                                          cai_scores=cai_scores_dict,
                                          tai_scores=tai_scores_dict,
                                          is_optimized=is_optimized)

        return organism_object
