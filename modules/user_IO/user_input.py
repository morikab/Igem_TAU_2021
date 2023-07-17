import json

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
    def run_module(cls, user_inp_raw: typing.Dict, run_summary: RunSummary) -> models.UserInput:
        logger.info('##########################')
        logger.info('# USER INPUT INFORMATION #')
        logger.info('##########################')
        return cls._parse_input(user_input=user_inp_raw, run_summary=run_summary)

    @classmethod
    def _parse_input(cls, user_input: typing.Dict[str, typing.Any], run_summary: RunSummary) -> models.UserInput:
        """
        :param user_input: in the following format
        {   'tuning_param': 0.5,
            'optimization_method': 'hill_climbing_average',
            'optimization_cub_index': 'CAI',
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
        optimization_cub_index = models.OptimizationCubIndex(user_input["optimization_cub_index"]) if \
            user_input.get("optimization_cub_index") else None
        optimization_method = models.OptimizationMethod(user_input["optimization_method"]) if \
            user_input.get("optimization_method") else None
        tuning_parameter = user_input["tuning_param"]
        clusters_count = user_input["clusters_count"]
        output_path = user_input.get("output_path")

        organisms_list = cls._parse_organisms_list(organisms_input_list=user_input["organisms"],
                                                   optimization_cub_index=optimization_cub_index)

        orf_sequence = cls._parse_orf_sequence(user_input)
        logger.info(F"Open reading frame sequence for optimization is: {orf_sequence}")

        user_input = models.UserInput(organisms=organisms_list,
                                      sequence=orf_sequence,
                                      tuning_parameter=tuning_parameter,
                                      optimization_method=optimization_method,
                                      optimization_cub_index=optimization_cub_index,
                                      clusters_count=clusters_count,
                                      output_path=output_path)

        run_summary.add_to_run_summary("user_input", user_input.summary)

        return user_input

    @classmethod
    def _parse_orf_sequence(cls, user_input: typing.Dict[str, typing.Any]) -> str:
        sequence = user_input.get("sequence")
        if sequence is not None:
            return sequence

        orf_fasta_file_path = user_input["sequence_file_path"]
        logger.info(F"Sequence to be optimized given in the following file: {orf_fasta_file_path}")

        try:
            return str(SeqIO.read(orf_fasta_file_path, "fasta").seq)
        except:
            raise ValueError(
                F"Error in orf sequence .fasta file: {orf_fasta_file_path}. Make sure you inserted a valid .fasta file "
                F"containing a single record")

    @classmethod
    def _parse_organisms_list(cls,
                              organisms_input_list: typing.Dict[str, typing.Any],
                              optimization_cub_index: models.OptimizationCubIndex) -> typing.List[models.Organism]:
        organisms_list = []
        organisms_names = set()
        for organism_key, organism_input in organisms_input_list.items():
            try:
                organism = cls._parse_single_organism_input(organism_input=organism_input,
                                                            optimization_cub_index=optimization_cub_index)
            except Exception as e:
                raise ValueError(f"Error in organism input: {organism_key}, re-check your input")
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

        return organisms_list

    @staticmethod
    def _parse_single_organism_input(organism_input: typing.Dict[str, typing.Any],
                                     optimization_cub_index: models.OptimizationCubIndex) -> models.Organism:
        gb_path = organism_input["genome_path"]

        # FIXME - delete
        is_optimized = organism_input["optimized"]
        parsed_organism_file_name = f"{gb_path.strip('.gb')}_{is_optimized}_parsed"
        parsed_organism_file = parsed_organism_file_name + ".json"

        # FIXME - delete
        # if os.path.exists(parsed_organism_file):
        #     with open(parsed_organism_file) as org_file:
        #         organism_data = json.load(org_file)
        #         return models.Organism(name=organism_data["name"],
        #                                cai_profile=organism_data["cai_weights"],
        #                                tai_profile=organism_data["tai_weights"],
        #                                cai_scores=organism_data["cai_scores"],
        #                                tai_scores=organism_data["tai_scores"],
        #                                reference_genes=organism_data["reference_genes"],
        #                                is_optimized=organism_data["is_wanted"],
        #                                optimization_priority=organism_data["optimization_priority"])

        # FIXME - end

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

        cai_weights = None
        cai_scores_dict = None
        tai_weights = None
        tai_scores_dict = None
        reference_genes = None
        if optimization_cub_index.is_codon_adaptation_index:
            cai_weights, reference_genes = calculate_cai_weights_for_input(cds_dict, estimated_expression)
            cai_scores = general_geomean(sequence_lst=cds_dict.values(), weights=cai_weights)
            cai_scores_dict = {gene_names[i]: cai_scores[i] for i in range(len(gene_names))}

        if optimization_cub_index.is_trna_adaptation_index:
            # TODO - move the location of tai calculation to user_input module
            tai_weights = tai_from_tgcnDB(organism_name)
            tai_scores = general_geomean(sequence_lst=cds_dict.values(), weights=tai_weights)
            tai_scores_dict = {gene_names[i]: tai_scores[i] for i in range(len(gene_names))}

        organism_object = models.Organism(name=organism_name,
                                          cai_profile=cai_weights,
                                          tai_profile=tai_weights,
                                          cai_scores=cai_scores_dict,
                                          tai_scores=tai_scores_dict,
                                          reference_genes=reference_genes,
                                          is_optimized=is_optimized,
                                          optimization_priority=optimization_priority)

        # FIXME - delete
        org_summary = organism_object.summary
        org_summary["cai_scores"] = cai_scores_dict
        org_summary["tai_scores"] = tai_scores_dict
        org_summary["cds_dict"] = cds_dict
        org_summary["reference_genes"] = reference_genes
        # with open(parsed_organism_file_name+".fasta", "w") as organism_fasta_file:
        write_fasta(fid=parsed_organism_file_name, list_seq=list(cds_dict.values()), list_name=list(cds_dict.keys()))

        with open(parsed_organism_file, "w") as organism_file:
            json.dump(org_summary, organism_file)
        # FIXME - end

        logger.info(
            F"name={organism_object.name}, cai_std={organism_object.cai_std}, cai_avg={organism_object.cai_avg}")
        return organism_object
