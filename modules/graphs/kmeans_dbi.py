import os
import shutil
import time
from pathlib import Path
import typing
from logger_factory.logger_factory import LoggerFactory
from modules.testing_for_modules import generate_testing_data


# Create clean artifacts directory
artifacts_directory = Path(os.path.join(str(Path(__file__).parent.resolve()), "../artifacts"))
if artifacts_directory.exists() and artifacts_directory.is_dir():
    shutil.rmtree(artifacts_directory)
artifacts_directory.mkdir(parents=True, exist_ok=True)
from modules import user_IO, ORF, sequence_family
from modules.stats.evaluation import ZscoreModule

logger = LoggerFactory.get_logger()


def run_modules(user_input_dict: typing.Optional[typing.Dict[str, typing.Any]] = None):
    user_inp_raw = user_input_dict

    user_input = user_IO.UserInputModule.run_module(user_inp_raw)

    # no clustering
    final_cds = ORF.ORFModule.run_module(user_input, 'cai', optimization_method=user_input.optimization_method)
    avg_opt_index, mean_deopt_index, avg_opt_index, weakest_score = \
        ZscoreModule.run_module(final_cds, user_input, 'cai')

    # with clustering
    opt_indexes = []
    clustered_user_inputs = sequence_family.SequenceFamilyModule.run_module(user_input)
    for input_cluster in clustered_user_inputs:
        final_cds = ORF.ORFModule.run_module(input_cluster, 'cai', optimization_method=user_input.optimization_method)
        c_mean_opt_index, c_mean_deopt_index, c_optimization_index, c_weakest_score = \
            ZscoreModule.run_module(final_cds, input_cluster, 'cai')
        opt_indexes.append(c_optimization_index)

    avg_c_opt_index = sum(opt_indexes)/len(opt_indexes)

    return avg_c_opt_index, avg_opt_index


if __name__ == "__main__":
    current_directory = Path(__file__).parent.resolve()
    base_path = os.path.join(Path(current_directory).parent.resolve(), "example_data")

    no_c_scores = []
    c_scores = []

    for n_organisms in range(10,31,4):
        print('############', n_organisms)
        default_user_inp_raw = generate_testing_data(n_organisms=n_organisms, percent_optimized=0.5)

        tic = time.time()
        avg_c_opt_index, avg_opt_index = run_modules(user_input_dict=default_user_inp_raw)
        toc = time.time()
        modules_run_time = toc - tic
        print('Total modules run time: ', modules_run_time)
        no_c_scores.append(avg_opt_index)
        c_scores.append(avg_c_opt_index)

    print(c_scores)
    print(no_c_scores)

    print('################')
    for i in range(5,31,5):
        print(i)

    print('################')
    for i in no_c_scores:
        print(i)

    print('################')
    for i in c_scores:
        print(i)

