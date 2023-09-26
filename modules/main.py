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
# if artifacts_directory.exists() and artifacts_directory.is_dir():
#     shutil.rmtree(artifacts_directory)
artifacts_directory.mkdir(parents=True, exist_ok=True)

from modules import user_IO, ORF, sequence_family
from modules.evaluation import models as evaluation_models
from modules.evaluation.evaluation import EvaluationModule
from modules import models

logger = LoggerFactory.get_logger()


def run_input_processing(user_input_dict: typing.Optional[typing.Dict[str, typing.Any]] = None) -> models.UserInput:
    run_summary = RunSummary()
    return user_IO.UserInputModule.run_module(user_inp_raw=user_input_dict, run_summary=run_summary)


def run_modules(user_input_dict: typing.Dict[str, typing.Any],
                should_run_output_module: bool = True) -> typing.Dict[str, typing.Any]:
    run_summary = RunSummary()
    final_output = {}
    try:
        before_parsing_input = time.time()
        user_input = user_IO.UserInputModule.run_module(user_inp_raw=user_input_dict, run_summary=run_summary)
        after_parsing_input = time.time()
        logger.info(F"Total input processing time: {after_parsing_input-before_parsing_input}")
        # ####################################### Family of sequences #####################################
        # in this part, user input is split into different inputs according to the sequence family theory
        clustered_user_inputs = sequence_family.SequenceFamilyModule.run_module(user_input)
        evaluation_results = []
        for input_cluster in clustered_user_inputs:
            evaluation_result = run_orf_optimization(user_input=input_cluster, run_summary=run_summary)
            evaluation_results.append(evaluation_result)

        # ###################################### Output Handling ##########################################
        output_path = user_input.output_path or str(artifacts_directory)
        run_summary.save_run_summary(output_path)

        final_output = run_summary.get()

        # TODO - handle multiple results in output generation module
        evaluation_result = evaluation_results[0]
        if should_run_output_module:
            final_output["zip_output_file_path"] = user_IO.UserOutputModule.run_module(
                cds_sequence=evaluation_result.sequence,
                output_path=output_path,
            )
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
        tai_evaluation_results: typing.Optional[typing.Sequence[evaluation_models.EvaluationModuleResult]],
        cai_evaluation_results: typing.Optional[typing.Sequence[evaluation_models.EvaluationModuleResult]],
) -> evaluation_models.EvaluationModuleResult:
    tai_evaluation_results = tai_evaluation_results or []
    cai_evaluation_results = cai_evaluation_results or []
    evaluation_results = tai_evaluation_results + cai_evaluation_results

    best_evaluation_result = None
    for evaluation_result in evaluation_results:
        logger.info(F"Analyzing evaluation result: {evaluation_result.summary}")
        if best_evaluation_result is None:
            best_evaluation_result = evaluation_result
            continue
        # TODO - extend to allow comparison by other scores (make it an argument)
        if evaluation_result.average_distance_score > best_evaluation_result.average_distance_score:
            best_evaluation_result = evaluation_result

    return best_evaluation_result


def run_orf_optimization(user_input: models.UserInput,
                         run_summary: RunSummary) -> evaluation_models.EvaluationModuleResult:
    optimization_cub_index = user_input.optimization_cub_index
    optimization_method = user_input.optimization_method
    tai_evaluation_results = None
    cai_evaluation_results = None

    if optimization_cub_index.is_trna_adaptation_index:
        trna_adaptation_index = optimization_cub_index.trna_adaptation_index
        cds_nt_final_tai = ORF.ORFModule.run_module(user_input=user_input,
                                                    optimization_cub_index=trna_adaptation_index,
                                                    optimization_method=optimization_method,
                                                    run_summary=run_summary)
        tai_evaluation_results = [
            EvaluationModule.run_module(final_sequence=cds_nt_tai,
                                        user_input=user_input,
                                        optimization_cub_index=trna_adaptation_index,
                                        run_summary=run_summary) for cds_nt_tai in cds_nt_final_tai
        ]

    if optimization_cub_index.is_codon_adaptation_index:
        codon_adaptation_index = optimization_cub_index.codon_adaptation_index
        cds_nt_final_cai = ORF.ORFModule.run_module(user_input=user_input,
                                                    optimization_cub_index=codon_adaptation_index,
                                                    optimization_method=optimization_method,
                                                    run_summary=run_summary)

        cai_evaluation_results = [
            EvaluationModule.run_module(final_sequence=cds_nt_cai,
                                        user_input=user_input,
                                        optimization_cub_index=codon_adaptation_index,
                                        run_summary=run_summary) for cds_nt_cai in cds_nt_final_cai
        ]

    evaluation_result = choose_orf_optimization_result(tai_evaluation_results=tai_evaluation_results,
                                                       cai_evaluation_results=cai_evaluation_results)

    logger.info(f"Evaluation result: {evaluation_result.summary}")
    return evaluation_result


def run_orf_module(user_input_dict: typing.Optional[typing.Dict[str, typing.Any]]):
    run_summary = RunSummary()
    user_input = user_IO.UserInputModule.run_module(user_inp_raw=user_input_dict, run_summary=run_summary)

    optimization_cub_index = user_input.optimization_cub_index
    optimization_method = user_input.optimization_method

    if optimization_cub_index.is_trna_adaptation_index:
        logger.info("tAI information:")
        trna_adaptation_index = optimization_cub_index.trna_adaptation_index
        return ORF.ORFModule.run_module(user_input=user_input,
                                        optimization_cub_index=trna_adaptation_index,
                                        optimization_method=optimization_method,
                                        run_summary=run_summary)

    codon_adaptation_index = optimization_cub_index.codon_adaptation_index
    return ORF.ORFModule.run_module(user_input=user_input,
                                    optimization_cub_index=codon_adaptation_index,
                                    optimization_method=optimization_method,
                                    run_summary=run_summary)
