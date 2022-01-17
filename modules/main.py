import os
import shutil
import time
import traceback
from pathlib import Path
import typing

from modules.logger_factory import LoggerFactory

# Create clean artifacts directory
artifacts_directory = Path(os.path.join(str(Path(__file__).parent.resolve()), "artifacts"))
if artifacts_directory.exists() and artifacts_directory.is_dir():
    shutil.rmtree(artifacts_directory)
artifacts_directory.mkdir(parents=True, exist_ok=True)

from modules import user_IO, RE, ORF
from modules.stats.evaluation import ZscoreModule
from modules import models

logger = LoggerFactory.create_logger("main")

current_directory = Path(__file__).parent.resolve()
base_path = os.path.join(Path(current_directory).parent.resolve(), "example_data")
default_user_inp_raw = {
    'sequence': os.path.join(base_path, 'mCherry_original.fasta'),
    'selected_promoters': None,
    'tuning_param': 0.75,
    'organisms': {
        'deopt2': {'genome_path': os.path.join(base_path, 'Sulfolobus acidocaldarius.gb'),
                   'optimized': False,
                   'expression_csv': None},
        'blabla': {'genome_path': os.path.join(base_path, 'Mycobacterium tuberculosis.gb'),
                   'optimized': True,
                   'expression_csv': None},
    }
}


def run_modules(user_input_dict: typing.Optional[typing.Dict[str, typing.Any]] = None,
                model_preferences_dict: typing.Optional[typing.Dict[str, str]] = None):
    user_inp_raw = user_input_dict or default_user_inp_raw

    model_preferences = models.ModelPreferences.init_from_dictionary(
        model_preferences_dict
    ) if model_preferences_dict is not None else models.ModelPreferences.init_from_config()

    try:
        before_parsing_input = time.time()
        user_input = user_IO.UserInputModule.run_module(user_inp_raw)
        after_parsing_input = time.time()

        logger.info(F"Total input processing time: {after_parsing_input-before_parsing_input}")

        # ############################################ unit 1 ############################################
        final_cds = None
        optimization_index = None
        weakest_score = None

        if model_preferences.restriction_enzymes or model_preferences.translation:
            final_cds, optimization_index, weakest_score = unit1(user_input, model_preferences)
        ##################################################################################################

        # TODO - get zip_directory from the user

        # TODO - change the user output to remove unneeded parameters
        final_output, zip_file_path = user_IO.UserOutputModule.run_module(cds_sequence=final_cds,
                                                                          zscore=optimization_index,
                                                                          weakest_score=weakest_score,
                                                                          p_name=None,
                                                                          native_prom=None,
                                                                          synth_promoter=None,
                                                                          evalue=None,
                                                                          zip_directory=str(artifacts_directory))
        logger.info("Final output: %s, zip_file_path: %s", final_output, zip_file_path)
    except Exception as e:
        exception_str = traceback.format_exc()
        final_output = {
            'error_message': exception_str,
        }
        zip_file_path = None
        logger.error("Encountered unknown error when running modules. Error message: %s", exception_str)

    return final_output, zip_file_path


def unit1(user_input: models.UserInput, model_preferences: models.ModelPreferences):
    if model_preferences.translation:
        optimization_func = model_preferences.translation_function
        try:
            # both CAI and tAI, select the one with the best optimization index tai optimization
            logger.info('tAI information:')
            cds_nt_final_tai = ORF.ORFModule.run_module(user_input, 'tai', optimization_func)

            if model_preferences.restriction_enzymes:
                cds_nt_final_tai = RE.REModule.run_module(user_input, cds_nt_final_tai)
            tai_mean_opt_index, tai_mean_deopt_index, tai_optimization_index, tai_weakest_score = \
                ZscoreModule.run_module(cds_nt_final_tai, user_input, optimization_type='tai')

            logger.info(f'Sequence:\n{cds_nt_final_tai}')
            logger.info(f'Optimized sequences score: {tai_mean_opt_index}, '
                        f'deoptimized sequence score: {tai_mean_deopt_index}')
            logger.info(f'Final optimization score: {tai_optimization_index}')

            # cai optimization
            logger.info('CAI information:')
            cds_nt_final_cai = ORF.ORFModule.run_module(user_input, 'cai', optimization_func)

            if model_preferences.restriction_enzymes:
                # todo: run both of them together to save time, or split creation of enzyme dict and the actual
                #  optimization (seems like a better solution)
                cds_nt_final_cai = RE.REModule.run_module(user_input, cds_nt_final_cai)
            cai_mean_opt_index, cai_mean_deopt_index, cai_optimization_index, cai_weakest_score = \
                ZscoreModule.run_module(cds_nt_final_cai, user_input, optimization_type='cai')

            logger.info(f'Sequence:\n{cds_nt_final_cai}')
            logger.info(f'Optimized sequences score: {cai_mean_opt_index}, '
                        f'deoptimized sequence score: {cai_mean_deopt_index}')
            logger.info(f'Final optimization score: {cai_optimization_index}')

            if cai_optimization_index > tai_optimization_index:
                logger.info('CAI sequence was selected')
                final_cds = cds_nt_final_cai
                optimization_index = cai_optimization_index
                mean_opt_index = cai_mean_opt_index
                mean_deopt_index = cai_mean_deopt_index
                weakest_score = cai_weakest_score
            else:
                logger.info('tAI sequence was selected')
                final_cds = cds_nt_final_tai
                optimization_index = tai_optimization_index
                mean_opt_index = tai_mean_opt_index
                mean_deopt_index = tai_mean_deopt_index
                weakest_score = tai_weakest_score

        except:
            logger.info('CAI information:')
            final_cds = ORF.ORFModule.run_module(user_input, 'cai', optimization_type=optimization_func)
            if model_preferences.restriction_enzymes:

                final_cds = RE.REModule.run_module(user_input, final_cds)
            mean_opt_index, mean_deopt_index, optimization_index, weakest_score =\
                ZscoreModule.run_module(final_cds, user_input, 'cai')

    else:
        final_cds = RE.REModule.run_module(user_input, user_input.sequence)
        mean_opt_index, mean_deopt_index, optimization_index, weakest_score = \
            ZscoreModule.run_module(final_cds, user_input, 'cai')

    logger.info(f'Sequence:\n{final_cds}')
    logger.info(f'Optimized sequences score: {mean_opt_index}, deoptimized sequence score: {mean_deopt_index}')
    logger.info(f'Weakest link score: {weakest_score}')
    logger.info(f'Final optimization score: {optimization_index}')
    return final_cds, optimization_index, weakest_score


if __name__ == "__main__":
    tic = time.time()
    run_modules()
    toc = time.time()
    modules_run_time = toc - tic
    print('Total modules run time: ', modules_run_time)
