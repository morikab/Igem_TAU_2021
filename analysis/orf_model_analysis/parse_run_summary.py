import json
import pandas as pd
import typing


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


if __name__ == "__main__":
    compare_initial_and_final_cds_codons_cai_weights()
    # parse_zscore_files_of_bulk_and_single()
