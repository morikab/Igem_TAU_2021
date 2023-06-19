import json
import pandas as pd
import os
import typing
from collections import defaultdict
from pathlib import Path

from modules.shared_functions_and_vars import nt_to_aa


def get_iterations_scores(zscore_json: typing.Dict) -> typing.Sequence[float]:
    return [x["sequence_score"] for x in zscore_json["orf"]["iterations_summary"]]


def parse_zscore_files_of_bulk_and_single() -> None:
    with open(r"C:\projects\Igem_TAU_2021_moran\analysis\orf_model_analysis\results_human\lcl-NC_000001.11_cds_NP_003386.1_12075\CAI_zscore_single_aa_average_ecoli_opt_True_o9Z7\run_summary.json") as summary_file:
        zscore_single_aa = json.load(summary_file)

    with open(r"C:\projects\Igem_TAU_2021_moran\analysis\orf_model_analysis\results_human\lcl-NC_000001.11_cds_NP_003386.1_12075\CAI_zscore_bulk_aa_average_ecoli_opt_True_dPzw\run_summary.json") as summary_file:
        zscore_bulk_aa = json.load(summary_file)

    # scores
    # for zscore_json in [zscore_single_aa, zscore_bulk_aa]:
    #     print("****** scores_per_iteration ****** ")
    #     for score in get_iterations_scores(zscore_json):
    #         print(score)

    # Selected codons
    for iteration in zscore_single_aa["orf"]["iterations_summary"]:
        print(iteration["selected_codons"][0][1])


def compare_initial_and_final_cds_codons_cai_weights() -> None:
    results_file = r"C:\projects\Igem_TAU_2021_moran\analysis\orf_model_analysis\results\arabidopsis\test_iteration\CAI_zscore_bulk_aa_average_30_Qmuv/run_summary.json"

    with open(results_file, "r") as summary_file:
        summary = json.load(summary_file)

    organisms_list = summary["user_input"]["organisms"]

    initial_sequence = summary["user_input"]["sequence"]
    initial_codons = [initial_sequence[i:i+3] for i in range(0, len(initial_sequence), 3)]

    # initial_scores = {i: [organism["cai_weights"][initial_codons[i]] for organism in organisms_list] for i
    #                   in range(len(initial_codons))}

    final_sequence = summary["evaluation"]["final_sequence"]
    final_codons = [final_sequence[i:i + 3] for i in range(0, len(final_sequence), 3)]

    scores = {}
    for i in range(len(initial_codons)):
        initial_scores = [organism["cai_weights"][initial_codons[i]] for organism in organisms_list]
        final_scores = [organism["cai_weights"][final_codons[i]] for organism in organisms_list]

        initial_scores.extend(final_scores)
        scores[i] = initial_scores

    df = pd.DataFrame(scores)
    df.to_csv("codon_heat_map.csv")


def compare_single_codon_method_scores() -> None:
    root_dir = r"C:\projects\Igem_TAU_2021_moran\analysis\orf_model_analysis\results_human\mcherry"
    filename = "run_summary.json"

    columns = defaultdict(list)

    are_initial_scores_initialized = False
    for root, dirs, files in os.walk(root_dir):
        for file in files:
            if file == filename:
                file_path = os.path.join(root, file)
                with open(file_path, "r") as json_file:
                    summary = json.load(json_file)
                if not are_initial_scores_initialized:
                    for organism in summary["user_input"]["organisms"]:
                        _update_initial_cub_scores(organism=organism, columns=columns, cub_index="CAI")
                        _update_initial_cub_scores(organism=organism, columns=columns, cub_index="tAI")
                    are_initial_scores_initialized = True

                columns["name"].append(Path(root).name)
                columns["cub_index"].append(summary["user_input"]["optimization_cub_index"])

                flattened_scores = {}
                for dictionary in summary["orf_debug"].values():
                    flattened_scores.update(dictionary)

                for codon in nt_to_aa.keys():
                    columns[codon].append(flattened_scores[codon])

    df = pd.DataFrame(columns)
    df.to_csv(os.path.join(root_dir, "single_codon_methods.csv"))


def _update_initial_cub_scores(organism: typing.Dict, columns: typing.Dict, cub_index: str):
    columns["name"].append(organism["name"])
    columns["cub_index"].append(cub_index)
    for codon in nt_to_aa.keys():
        columns[codon].append(organism[F"{cub_index.lower()}_weights"][codon])

from collections import ChainMap
def fix_res(file_path):
    with open(file_path, "r") as res:
        content = json.load(res)

    with open(file_path, "w") as res:
        new_content = dict(ChainMap(*content))
        json.dump(new_content, res)


if __name__ == "__main__":
    fix_res("res.json")
    # compare_single_codon_method_scores()
    # compare_initial_and_final_cds_codons_cai_weights()
    # parse_zscore_files_of_bulk_and_single()
