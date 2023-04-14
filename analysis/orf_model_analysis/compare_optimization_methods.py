import time

from input_testing_data.generate_input_testing_data_for_modules import \
    generate_testing_data_for_comparing_with_previous_algorithm
from modules.main import run_modules

if __name__ == "__main__":
    for optimization_method in ["single_codon_ratio", "single_codon_diff", "zscore_single_aa_average",
                                "zscore_bulk_aa_average", "zscore_single_aa_weakest_link",
                                "zscore_bulk_aa_weakest_link"]:
        for direction in [True, False]:
            tic = time.time()
            default_user_inp_raw = generate_testing_data_for_comparing_with_previous_algorithm(
                optimization_method=optimization_method,
                optimization_cub_index="CAI",
                clusters_count=1,
                tuning_param=0.5,
                is_ecoli_optimized=direction,
            )
            run_modules(default_user_inp_raw)
            toc = time.time()
            modules_run_time = toc - tic
            print(F"Modules run time for {optimization_method} is: {modules_run_time}")

            break   # FIXME
