import typing

from modules import models
from modules.run_summary import RunSummary
from modules.user_IO.input_functions import *
from modules.ORF.TAI import TAI
from modules.ORF.calculating_cai import general_geomean

logger = LoggerFactory.get_logger()


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
    def _parse_input(cls, user_input: typing.Dict[str, typing.Any]) -> models.UserInput:
        """
        :param user_input: in the following format
        {   'tuning_param': 0.5,
            'optimization_method': 'hill_climbing_average',
            'optimization_cub_score': 'CAI',
            'clustering_num': 3,
            'organisms: {
                'ecoli': {
                    'genome_path': '.gb',
                    'genes_HE': '.csv', ==> Optional
                    'optimized': True,
                },
                'bacillus': {
                    'genome_path': '.gb',
                    'genes_HE': '.csv', ==> Optional
                    'optimized': False,
                }, ...
            }
        }
        """
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

        # Normalize prioritization weights - from user's defined values (1-100) to normalized values in range (0, 1)
        total_optimized_weights = sum([organism.optimization_priority for organism in organisms_list if
                                       organism.is_optimized])
        total_deoptimized_weights = sum([organism.optimization_priority for organism in organisms_list if
                                         not organism.is_optimized])
        for organism in organisms_list:
            if organism.is_optimized:
                organism.optimization_priority /= total_optimized_weights
            else:
                organism.optimization_priority /= total_deoptimized_weights
            logger.info(f"{organism.name} has weight of {organism.optimization_priority}")

        # Read ORF sequence
        orf_fasta_fid = user_input['sequence']
        orf_seq = None
        if orf_fasta_fid is not None:
            try:
                orf_seq = str(SeqIO.read(orf_fasta_fid, 'fasta').seq)
            except:
                raise ValueError(
                    f'Error in protein .fasta file: {orf_fasta_fid}, make sure you inserted an undamaged .fasta file '
                    f'containing a single recored')
        logger.info(f'\n\nSequence to be optimized given in the following file {orf_fasta_fid}')
        logger.info(f'containing this sequence: {orf_seq}')

        tuning_parameter = user_input["tuning_param"]
        optimization_method = models.OptimizationMethod(user_input["optimization_method"]) if \
            user_input.get("optimization_method") else None
        # FIXME - change default to None and add the option to define cub score in ui
        optimization_cub_score = models.OptimizationCubScore(user_input["optimization_cub_score"]) if \
            user_input.get("optimization_cub_score") else models.OptimizationCubScore.codon_adaptation_index
        clusters_count = user_input["clusters_count"]
        zip_directory = user_input.get("output_path")

        user_input = models.UserInput(organisms=organisms_list,
                                      sequence=orf_seq,
                                      tuning_parameter=tuning_parameter,
                                      optimization_method=optimization_method,
                                      optimization_cub_score=optimization_cub_score,
                                      clusters_count=clusters_count,
                                      zip_directory=zip_directory)

        RunSummary.add_to_run_summary("user_input", user_input.summary)

        return user_input

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
        is_optimized = organism_input["optimized"]
        if is_optimized:
            logger.info("Organism is optimized")
        else:
            logger.info("Organism is deoptimized")
        optimization_priority = organism_input.get("optimization_priority") or DEFAULT_ORGANISM_PRIORITY

        cds_dict, estimated_expression = extract_gene_data(gb_path, exp_csv_fid)
        logger.info(f'Number of genes: {len(cds_dict)}')
        gene_names = list(cds_dict.keys())

        # TODO - extract weights by the value of optimization_cub_score to reduce running times of the code
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
                                          is_optimized=is_optimized,
                                          optimization_priority=optimization_priority)

        logger.info(F"name={organism_object.name}, cai_std={organism_object.cai_std}, cai_avg={organism_object.cai_avg}")
        return organism_object
