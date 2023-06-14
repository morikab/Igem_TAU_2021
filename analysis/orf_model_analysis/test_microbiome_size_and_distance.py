import argparse
import json
import os
import pandas as pd
import typing
from pathlib import Path

from Bio import SeqIO


from input_testing_data.generate_input_testing_data_for_modules import generate_testing_data
from input_testing_data.generate_input_testing_data_for_modules import get_organisms_for_testing
from modules.main import run_input_processing
from modules.main import run_modules
from modules.main import run_orf_module


current_directory = Path(__file__).parent.resolve()
base_path = os.path.join(Path(current_directory).parent.resolve(), "example_data")
DEFAULT_SEQUENCE_FILE_PATH = os.path.join(base_path, "mCherry_original.fasta")


def run_for_different_microbiome_size(output_path: str) -> None:
    wanted_hosts = []
    unwanted_hosts = []
    selected_hosts = wanted_hosts + unwanted_hosts
    max_microbiome_size = 30
    added_organisms_per_iteration = 5
    while len(selected_hosts) <= max_microbiome_size:
        added_wanted_hosts, added_unwanted_hosts = get_organisms_for_testing(
            organisms_count=added_organisms_per_iteration,
            excluded_genomes=selected_hosts,
        )
        wanted_hosts.extend(added_wanted_hosts)
        unwanted_hosts.extend(added_unwanted_hosts)

        run_all_methods(orf_sequence_file=DEFAULT_SEQUENCE_FILE_PATH,
                        wanted_hosts=wanted_hosts,
                        unwanted_hosts=unwanted_hosts,
                        output_path=output_path)

        selected_hosts = wanted_hosts + unwanted_hosts


def run_for_sub_microbiome(output_path: str,
                           optimization_method: str,
                           optimization_cub_index: str,
                           fasta_file_path: str,
                           wanted_hosts: typing.Sequence[str] = None,
                           unwanted_hosts: typing.Sequence[str] = None,
                           start_seq: int = 0) -> None:
    microbiome_size = 4
    max_records = 500

    if wanted_hosts is None or unwanted_hosts is None:
        wanted_hosts, unwanted_hosts = get_organisms_for_testing(
            organisms_count=microbiome_size,
            excluded_genomes=[],
        )

    with open(fasta_file_path, "r") as fasta_handle:
        genome_dict = SeqIO.to_dict(SeqIO.parse(fasta_handle, "fasta"))

    endogenous_genes_organism = Path(fasta_file_path).name.strip("fasta").strip(".")
    configuration = f"wanted_{'_'.join(wanted_hosts)}_unwanted_{'_'.join(unwanted_hosts)}"

    results_dict = {}
    count = 0
    for gene_name, gene_sequence in genome_dict.items():
        if count < start_seq or count - start_seq > max_records:
            continue
        count += 1
        gene_sequence = str(gene_sequence.seq)
        if len(gene_sequence) % 3 != 0:
            print(F"Invalid length {len(gene_sequence)} for gene {gene_name}")
            continue
        results_dict[gene_name] = run_single_method_for_orf_sequence(
            optimization_method=optimization_method,
            optimization_cub_index=optimization_cub_index,
            orf_sequence=gene_sequence,
            output_path=os.path.join(output_path, configuration, endogenous_genes_organism, gene_name),
            wanted_hosts=wanted_hosts,
            unwanted_hosts=unwanted_hosts,
        )
    with open(
            os.path.join(output_path,
                         F"{optimization_cub_index}_{optimization_method}_{Path(fasta_file_path).name[:10]}_fasta_results.json"),
            "w") as results_file:
        json.dump(results_dict, results_file)
    with open(
            os.path.join(output_path,
                         F"{optimization_cub_index}_{optimization_method}_{Path(fasta_file_path).name[:10]}_organisms.json"),
            "w") as results_file:
        selected_organisms = {
            "wanted": wanted_hosts,
            "unwanted": unwanted_hosts,
        }
        json.dump(selected_organisms, results_file)


def run_all_methods(wanted_hosts: typing.Sequence[str],
                    unwanted_hosts: typing.Sequence[str],
                    orf_sequence: typing.Optional[str] = None,
                    orf_sequence_file: typing.Optional[str] = None,
                    output_path: typing.Optional[str] = None):
    for optimization_method in [
        # "single_codon_ratio", "single_codon_diff", "single_codon_weakest_link",
        # "zscore_single_aa_ratio",
        #"zscore_bulk_aa_ratio",
        # "zscore_single_aa_diff",
         "zscore_bulk_aa_diff",
        # "zscore_single_aa_weakest_link",
         #"zscore_bulk_aa_weakest_link",
    ]:
        for optimization_cub_index in [
            "CAI",
            # "tAI",
        ]:
            return run_single_method_for_orf_sequence(optimization_method=optimization_method,
                                                      optimization_cub_index=optimization_cub_index,
                                                      orf_sequence=orf_sequence,
                                                      orf_sequence_file=orf_sequence_file,
                                                      output_path=output_path,
                                                      wanted_hosts=wanted_hosts,
                                                      unwanted_hosts=unwanted_hosts)


def run_single_method_for_orf_sequence(optimization_method: str,
                                       optimization_cub_index: str,
                                       wanted_hosts: typing.Sequence[str],
                                       unwanted_hosts: typing.Sequence[str],
                                       orf_sequence: typing.Optional[str] = None,
                                       orf_sequence_file: typing.Optional[str] = None,
                                       output_path: typing.Optional[str] = None):
    default_user_inp_raw = generate_testing_data(
        optimization_method=optimization_method,
        optimization_cub_index=optimization_cub_index,
        wanted_hosts=wanted_hosts,
        unwanted_hosts=unwanted_hosts,
        tuning_param=0.5,
        sequence=orf_sequence,
        sequence_file_path=orf_sequence_file,
        output_path=output_path,
    )
    # return run_orf_module(default_user_inp_raw)

    modules_output = run_modules(default_user_inp_raw)
    if "error_message" in modules_output:
          raise Exception(F"Encountered error while running module: {modules_output['error_message']}")


def analyze_per_microbiome_size(results_directory: str) -> None:
    filename = "run_summary.json"
    prefix = "CAI_zscore_bulk_aa_average"

    results_summary = {}
    for root, dirs, files in os.walk(results_directory):
        for file in files:
            if file == filename:
                directory_name = Path(root).name
                if directory_name.startswith(prefix):
                    microbiome_size = int(directory_name.split(prefix)[1].strip("_").split("_")[0])
                    file_path = os.path.join(root, file)
                    with open(file_path, "r") as summary_file:
                        results_json = json.load(summary_file)
                        results_summary[microbiome_size] = {
                            "average_distance_score": results_json["evaluation"]["average_distance_score"],
                            "weakest_link_score": results_json["evaluation"]["weakest_link_score"],
                        }
    with open(os.path.join(results_directory, "summary.json"), "w") as results_summary_file:
        json.dump(results_summary, results_summary_file)


def group_summary_files(results_directory: str) -> None:
    filename = "summary.json"

    summary_df = pd.DataFrame()
    for root, dirs, files in os.walk(results_directory):
        for file in files:
            if file == filename:
                file_path = os.path.join(root, file)
                with open(file_path, "r") as results_file:
                    summary = json.load(results_file)

                    average_score_data = {F"average_score_{i}": [summary[str(i)]["average_distance_score"]]
                                          for i in range(5, 35, 5)}
                    weakest_link_score_data = {F"weakest_link_score_{i}": [summary[str(i)]["weakest_link_score"]]
                                               for i in range(5, 35, 5)}

                    average_score_data.update(weakest_link_score_data)
                    df = pd.DataFrame(average_score_data)
                summary_df = summary_df.append(df, ignore_index=True)

    summary_df.to_csv(os.path.join(results_directory, "summary.csv"))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Analysis script parser")
    parser.add_argument('-o', '--output', type=str, required=True, help="Output directory name")
    parser.add_argument('-f', '--fasta', type=str, help="Fasta file for orf sequences to run on (for endogenous run)")
    parser.add_argument('-m', '--method', type=str, help="Optimization method")
    parser.add_argument('-i', '--index', type=str, help="Optimization CUB index")
    parser.add_argument('--wanted', type=str, nargs='+', default=[], help="Wanted hosts")
    parser.add_argument('--unwanted', type=str, nargs='+', default=[], help="Unwanted hosts")
    parser.add_argument('-s', '--start', type=int, help="Start index from fasta file")
    args = parser.parse_args()

    run_for_sub_microbiome(
        output_path=args.output,
        optimization_method=args.method,
        optimization_cub_index=args.index,
        fasta_file_path=args.fasta,
        wanted_hosts=args.wanted,
        unwanted_hosts=args.unwanted,
        start_seq=args.start,
    )

    # analyze_per_microbiome_size(
    #     results_directory=rF"C:\projects\Igem_TAU_2021_moran\analysis\orf_model_analysis\results\arabidopsis\{args.output}",
    # )

    # group_summary_files(results_directory=r"C:\projects\Igem_TAU_2021_moran\analysis\orf_model_analysis\results\arabidopsis")

