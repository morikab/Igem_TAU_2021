import json
import typing


def get_iterations_scores(zscore_json: typing.Dict) -> typing.Sequence[float]:
    return [x["sequence_score"] for x in zscore_json["orf"]["iterations_summary"]]


def parse_zscore_files() -> None:
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


if __name__ == "__main__":
    parse_zscore_files()
