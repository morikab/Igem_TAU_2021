import json

from modules.run_summary import RunSummary
from modules.user_IO.input_functions import *
from modules.shared_functions_and_vars import DEFAULT_ORGANISM_PRIORITY
from modules.shared_functions_and_vars import write_fasta

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
        evaluation_score = models.EvaluationScore(user_input["evaluation_score"]) if \
            user_input.get("evaluation_score") else None
        tuning_parameter = user_input["tuning_param"]
        clusters_count = user_input["clusters_count"]
        output_path = user_input.get("output_path")

        orf_sequence = cls._parse_orf_sequence(user_input)
        logger.info(F"Open reading frame sequence for optimization is: {orf_sequence}")

        organisms_list = cls._parse_organisms_list(organisms_input_list=user_input["organisms"],
                                                   optimization_cub_index=optimization_cub_index)

        user_input = models.UserInput(organisms=organisms_list,
                                      sequence=orf_sequence,
                                      tuning_parameter=tuning_parameter,
                                      optimization_method=optimization_method,
                                      optimization_cub_index=optimization_cub_index,
                                      evaluation_score=evaluation_score,
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
                organism = cls._parse_single_organism_input(
                    organism_input=organism_input,
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

        is_optimized = organism_input["optimized"]
        logger.info(f"Organism is {'optimized' if is_optimized else 'de-optimized'}")

        gb_path = organism_input["genome_path"]

        # FIXME - delete
#         parsed_organism_file_name = f"{gb_path.strip('.gb')}_{is_optimized}_parsed"
#         parsed_organism_file = parsed_organism_file_name + ".json"

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
        try:
            gb_file = SeqIO.read(gb_path, format="gb")
            organism_name = " ".join(gb_file.description.split()[:2])
        except:
            raise ValueError(
                f'Error in genome GenBank file: {gb_path}, make sure you inserted an undamaged .gb file containing '
                f'the full genome sequence and annotations'
            )
        logger.info(f"Information about {organism_name}:")

        cds = extract_gene_data(genbank_path=gb_path)

        exp_csv_type = organism_input['expression_csv_type']
        exp_csv_fid = organism_input['expression_csv']
        estimated_expression = extract_gene_expression(cds=cds,
                                                       expression_csv_fid=exp_csv_fid,
                                                       expression_csv_type=exp_csv_type)

        cds_dict = {cds_record.name_and_function: cds_record.sequence for cds_record in cds}
        gene_names = list(cds_dict.keys())

        cai_weights = None
        cai_scores_dict = None
        tai_weights = None
        tai_scores_dict = None
        reference_genes = None

        if optimization_cub_index.is_codon_adaptation_index:
            reference_genes = get_reference_genes_for_cai(cds_dict, estimated_expression)
            cai = cb.scores.CodonAdaptationIndex(ref_seq=reference_genes.values(), ignore_stop=False)
            cai_weights = cai.weights.to_dict()
            cai_scores_dict = {gene_name: cai.get_score(cds_dict[gene_name]) for gene_name in gene_names}

        if optimization_cub_index.is_trna_adaptation_index:
            tai = calculate_tai_weights(organism_name)
            if tai is not None:
                tai_weights = tai.weights.to_dict()
                tai_scores_dict = {gene_name: tai.get_score(cds_dict[gene_name]) for gene_name in gene_names}

        optimization_priority = organism_input.get("optimization_priority") or DEFAULT_ORGANISM_PRIORITY
        organism_object = models.Organism(name=organism_name,
                                          cai_profile=cai_weights,
                                          tai_profile=tai_weights,
                                          cai_scores=cai_scores_dict,
                                          tai_scores=tai_scores_dict,
                                          reference_genes=list(reference_genes.keys()) if reference_genes else None,
                                          is_optimized=is_optimized,
                                          optimization_priority=optimization_priority)

        # FIXME - delete
#         org_summary = organism_object.summary
#         org_summary["cai_scores"] = cai_scores_dict
#         org_summary["tai_scores"] = tai_scores_dict
#         org_summary["cds_dict"] = cds_dict
#         org_summary["reference_genes"] = reference_genes
#         write_fasta(fid=organism_name, list_seq=list(cds_dict.values()), list_name=list(cds_dict.keys()))

#         with open(parsed_organism_file, "w") as organism_file:
#             json.dump(org_summary, organism_file)
        # FIXME - end

        if optimization_cub_index.is_codon_adaptation_index:
            logger.info(
                F"name={organism_object.name}, cai_std={organism_object.cai_std}, cai_avg={organism_object.cai_avg}")
        if optimization_cub_index.is_trna_adaptation_index:
            logger.info(
                F"name={organism_object.name}, tai_std={organism_object.tai_std}, tai_avg={organism_object.tai_avg}")
        return organism_object
