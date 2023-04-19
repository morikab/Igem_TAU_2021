import time

from input_testing_data.generate_input_testing_data_for_modules import \
    generate_testing_data_for_comparing_with_previous_algorithm
from modules.main import run_modules


def run_all_methods():
    for optimization_method in ["single_codon_ratio", "single_codon_diff", "zscore_single_aa_average",
                                "zscore_bulk_aa_average", "zscore_single_aa_weakest_link",
                                "zscore_bulk_aa_weakest_link"]:
        for direction in [True, False]:
            run_single_method(optimization_method=optimization_method, is_ecoli_optimized=direction)


def run_single_method(optimization_method: str, is_ecoli_optimized: bool) -> None:
    tic = time.time()
    default_user_inp_raw = generate_testing_data_for_comparing_with_previous_algorithm(
        optimization_method=optimization_method,
        optimization_cub_index="CAI",
        clusters_count=1,
        tuning_param=0.5,
        is_ecoli_optimized=is_ecoli_optimized,
    )
    run_modules(default_user_inp_raw)
    toc = time.time()
    modules_run_time = toc - tic
    print(F"Modules run time for {optimization_method} is: {modules_run_time}")


if __name__ == "__main__":
    # run_all_methods()
    run_single_method(optimization_method="zscore_single_aa_average", is_ecoli_optimized=True)

