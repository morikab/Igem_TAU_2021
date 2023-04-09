import time

from input_testing_data.generate_input_testing_data_for_modules import \
    generate_testing_data_for_comparing_with_previous_algorithm
from modules.main import run_modules

if __name__ == "__main__":
    tic = time.time()
    default_user_inp_raw = generate_testing_data_for_comparing_with_previous_algorithm(
        optimization_method="zscore_bulk_aa_average",
        optimization_cub_index="CAI",
        clusters_count=1,
        tuning_param=0.5,
        is_ecoli_optimized=False,
    )
    run_modules(default_user_inp_raw)
    toc = time.time()
    modules_run_time = toc - tic
    print(F"Total modules run time: {modules_run_time}")
