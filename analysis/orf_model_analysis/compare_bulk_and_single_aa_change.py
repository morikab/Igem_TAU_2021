import time

from input_testing_data.generate_input_testing_data_for_modules import \
    generate_testing_data_for_comparing_with_previous_algorithm
from modules.main import run_modules

if __name__ == "__main__":
    tic = time.time()
    # TODO - we need to modify the code to adapt for running on genes from the proteome:
    # 1. Create another script for writing the orfs to different fasta files
    # 2. Create a loop of iterating the different fasta files, and invoke the run_modules on each fasta file for the
    # following variations:
    #   (a) bulk + no bulk
    #   (b) ecoli optimized and not optimized 
    default_user_inp_raw = generate_testing_data_for_comparing_with_previous_algorithm(
        optimization_method="single_codon_ratio",
        optimization_cub_index="CAI",
        clusters_count=1,
        tuning_param=0.5,
        is_ecoli_optimized=False,
    )
    run_modules(default_user_inp_raw)
    toc = time.time()
    modules_run_time = toc - tic
    print(F"Total modules run time: {modules_run_time}")
