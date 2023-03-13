import os
import shutil
import time
import traceback
from pathlib import Path
import typing

from logger_factory.logger_factory import LoggerFactory
from modules.run_summary import RunSummary

# Create clean artifacts directory
artifacts_directory = Path(os.path.join(str(Path(__file__).parent.resolve()), "artifacts"))
if artifacts_directory.exists() and artifacts_directory.is_dir():
    shutil.rmtree(artifacts_directory)
artifacts_directory.mkdir(parents=True, exist_ok=True)

from modules import user_IO, ORF, sequence_family
from modules.stats import models as evaluation_models
from modules.stats.evaluation import EvaluationModule
from modules import models

logger = LoggerFactory.get_logger()


def run_modules(user_input_dict: typing.Optional[typing.Dict[str, typing.Any]] = None):
    final_output = {}
    try:
        before_parsing_input = time.time()
        user_input = user_IO.UserInputModule.run_module(user_input_dict)
        after_parsing_input = time.time()
        logger.info(F"Total input processing time: {after_parsing_input-before_parsing_input}")
        # ####################################### Family of sequences #####################################
        # in this part, user input is split into different inputs according to the sequence family theory
        clustered_user_inputs = sequence_family.SequenceFamilyModule.run_module(user_input)
        evaluation_results = []
        for input_cluster in clustered_user_inputs:
            evaluation_result = run_orf_optimization(input_cluster)
            evaluation_results.append(evaluation_result)

        # ###################################### Output Handling ##########################################
        # TODO - think if this is indeed the best location
        RunSummary.save_run_summary(str(artifacts_directory))

        zip_directory = user_input.zip_directory or str(artifacts_directory)
        # TODO - handle multiple results in output generation module
        evaluation_result = evaluation_results[0]
        final_output = user_IO.UserOutputModule.run_module(cds_sequence=evaluation_result.sequence,
                                                           zscore=evaluation_result.optimization_index,
                                                           weakest_score=evaluation_result.weakest_score,
                                                           zip_directory=zip_directory)
    except:
        logger.error("Encountered unknown error when running modules.")
        exception_str = traceback.format_exc()
        final_output = {
            "error_message": exception_str,
        }
    finally:
        logger.info("Final output: %s", final_output)

    return final_output


def choose_orf_optimization_result(
        tai_evaluation_result: typing.Optional[evaluation_models.EvaluationModuleResult],
        cai_evaluation_result: typing.Optional[evaluation_models.EvaluationModuleResult],
) -> evaluation_models.EvaluationModuleResult:
    if tai_evaluation_result is None:
        return cai_evaluation_result
    if cai_evaluation_result is None:
        return tai_evaluation_result
    return cai_evaluation_result if cai_evaluation_result.optimization_index > tai_evaluation_result.optimization_index\
        else tai_evaluation_result


def run_orf_optimization(user_input: models.UserInput) -> evaluation_models.EvaluationModuleResult:
    optimization_cub_score = user_input.optimization_cub_score
    optimization_method = user_input.optimization_method
    tai_evaluation_result = None
    cai_evaluation_result = None

    if optimization_cub_score.is_trna_adaptation_score:
        logger.info("tAI information:")
        trna_adaptation_index_score = optimization_cub_score.trna_adaptation_index
        cds_nt_final_tai = ORF.ORFModule.run_module(user_input=user_input,
                                                    optimization_cub_score=trna_adaptation_index_score,
                                                    optimization_method=optimization_method)
        tai_evaluation_result = EvaluationModule.run_module(final_seq=cds_nt_final_tai,
                                                            user_input=user_input,
                                                            optimization_cub_score=trna_adaptation_index_score)

        logger.info(f"Sequence:\n{cds_nt_final_tai}")
        logger.info(f"Optimized sequences score: {tai_evaluation_result.mean_opt_index}, "
                    f"deoptimized sequence score: {tai_evaluation_result.mean_deopt_index}")
        logger.info(f"Final optimization score: {tai_evaluation_result.optimization_index}")

    if optimization_cub_score.is_codon_adaptation_score:
        logger.info("CAI information:")
        codon_adaptation_index_score = optimization_cub_score.codon_adaptation_index
        cds_nt_final_cai = ORF.ORFModule.run_module(user_input=user_input,
                                                    optimization_cub_score=codon_adaptation_index_score,
                                                    optimization_method=optimization_method)
        cai_evaluation_result = EvaluationModule.run_module(final_seq=cds_nt_final_cai,
                                                            user_input=user_input,
                                                            optimization_cub_score=codon_adaptation_index_score)

        logger.info(f"Sequence:\n{cds_nt_final_cai}")
        logger.info(f"Optimized sequences score: {cai_evaluation_result.mean_opt_index}, "
                    f"deoptimized sequence score: {cai_evaluation_result.mean_deopt_index}")
        logger.info(f"Final optimization score: {cai_evaluation_result.optimization_index}")

    evaluation_result = choose_orf_optimization_result(tai_evaluation_result=tai_evaluation_result,
                                                       cai_evaluation_result=cai_evaluation_result)

    logger.info(f"Sequence:\n{evaluation_result.sequence}")
    logger.info(f"Optimized sequences score: {evaluation_result.mean_opt_index}, "
                f"deoptimized sequence score: {evaluation_result.mean_deopt_index}")
    logger.info(f"Weakest link score: {evaluation_result.weakest_score}")
    logger.info(f"Final optimization score: {evaluation_result.optimization_index}")
    return evaluation_result
