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
# artifacts_directory.mkdir(parents=True, exist_ok=True)

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

        # FIXME - start
        run_summary_content = run_summary.get()
        
        if run_summary_content.get("orf") is None:
            raise ValueError(run_summary_content)
        
        final_output = {
            "initial_optimization_score": run_summary_content["orf"].get("initial_sequence_optimization_score"),
            "final_optimization_score": run_summary_content["orf"].get("final_sequence_optimization_score"),
            "average_distance_score": run_summary_content["evaluation"]["average_distance_score"],
            "weakest_link_score": run_summary_content["evaluation"]["weakest_link_score"],
        }
        
        # run_summary.save_run_summary(output_path)
        # final_output = run_summary.get()
        # FIXME - end 

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
        tai_evaluation_result: typing.Optional[evaluation_models.EvaluationModuleResult],
        cai_evaluation_result: typing.Optional[evaluation_models.EvaluationModuleResult],
) -> evaluation_models.EvaluationModuleResult:
    if tai_evaluation_result is None:
        return cai_evaluation_result
    if cai_evaluation_result is None:
        return tai_evaluation_result
    return cai_evaluation_result if cai_evaluation_result.average_distance_score > tai_evaluation_result.average_distance_score\
        else tai_evaluation_result


def run_orf_optimization(user_input: models.UserInput,
                         run_summary: RunSummary) -> evaluation_models.EvaluationModuleResult:
    optimization_cub_index = user_input.optimization_cub_index
    optimization_method = user_input.optimization_method
    tai_evaluation_result = None
    cai_evaluation_result = None

    if optimization_cub_index.is_trna_adaptation_index:
        logger.info("tAI information:")
        trna_adaptation_index = optimization_cub_index.trna_adaptation_index
        cds_nt_final_tai = ORF.ORFModule.run_module(user_input=user_input,
                                                    optimization_cub_index=trna_adaptation_index,
                                                    optimization_method=optimization_method,
                                                    run_summary=run_summary)
        tai_evaluation_result = EvaluationModule.run_module(final_sequence=cds_nt_final_tai,
                                                            user_input=user_input,
                                                            optimization_cub_index=trna_adaptation_index,
                                                            run_summary=run_summary)

        logger.info(f"Sequence:\n{cds_nt_final_tai}")
        logger.info(f"Final optimization score: {tai_evaluation_result.average_distance_score}")

    if optimization_cub_index.is_codon_adaptation_index:
        logger.info("CAI information:")
        codon_adaptation_index = optimization_cub_index.codon_adaptation_index
        cds_nt_final_cai = ORF.ORFModule.run_module(user_input=user_input,
                                                    optimization_cub_index=codon_adaptation_index,
                                                    optimization_method=optimization_method,
                                                    run_summary=run_summary)
        cai_evaluation_result = EvaluationModule.run_module(final_sequence=cds_nt_final_cai,
                                                            user_input=user_input,
                                                            optimization_cub_index=codon_adaptation_index,
                                                            run_summary=run_summary)

        logger.info(f"Sequence:\n{cds_nt_final_cai}")
        logger.info(f"Final optimization score: {cai_evaluation_result.average_distance_score}")

    evaluation_result = choose_orf_optimization_result(tai_evaluation_result=tai_evaluation_result,
                                                       cai_evaluation_result=cai_evaluation_result)

    logger.info(f"Sequence:\n{evaluation_result.sequence}")
    logger.info(f"Weakest link score: {evaluation_result.weakest_link_score}")
    logger.info(f"Final optimization score: {evaluation_result.average_distance_score}")
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
