import json
import os

import openpyxl


aa_to_cell_column = {
    "C": "L",
    "D": "M",
    "S": "N",
    "Q": "O",
    "M": "P",
    "N": "Q",
    "P": "R",
    "K": "S",
    "_": "T",
    "T": "U",
    "F": "V",
    "A": "W",
    "G": "X",
    "I": "Y",
    "L": "Z",
    "H": "AA",
    "R": "AB",
    "W": "AC",
    "V": "AD",
    "E": "AE",
    "Y": "AF",
}


def generate_summary(results_directory: str) -> None:
    filename = "run_summary.json"

    workbook = openpyxl.Workbook()
    worksheet = workbook.active

    initialize_column_headers(worksheet)

    row_offset = 2
    for root, dirs, files in os.walk(results_directory):
        for file in files:
            if file == filename:
                filepath = os.path.join(root, file)
                parse_summary_file(file_path=filepath, worksheet=worksheet, row_offset=row_offset)
                row_offset += 1

    workbook.save("summary.xlsx")


def initialize_column_headers(worksheet) -> None:
    worksheet["A1"] = "optimization_method"

    worksheet["B1"] = "is_ecoli_wanted"
    worksheet["C1"] = "initial_gene_cai_score_ecoli"
    worksheet["D1"] = "final_gene_cai_score_ecoli"
    worksheet["E1"] = "ecoli_dist_score"

    worksheet["F1"] = "is_bacillus_wanted"
    worksheet["G1"] = "initial_gene_cai_score_bacillus"
    worksheet["H1"] = "final_gene_cai_score_bacillus"
    worksheet["I1"] = "bacillus_dist_score"

    worksheet["J1"] = "final_average_distance_score"
    worksheet["K1"] = "final_weakest_link_score"

    for aa, cell in aa_to_cell_column.items():
        worksheet[F"{cell}1"] = aa

    worksheet["AG1"] = "number_of_iterations"
    worksheet["AH1"] = "run_time (seconds)"


def parse_summary_file(file_path: str, worksheet, row_offset: int) -> None:
    with open(file_path, "r") as summary_file:
        summary = json.load(summary_file)

        organisms = summary["evaluation"]["organisms"]
        first_organism = organisms[0]
        if first_organism["name"] == "Escherichia coli":
            e_coli = organisms[0]
            bacillus = organisms[1]
        else:
            e_coli = organisms[1]
            bacillus = organisms[0]

        worksheet[f"A{row_offset}"] = summary["user_input"]["optimization_method"]

        worksheet[f"B{row_offset}"] = e_coli["is_wanted"]
        worksheet[f"C{row_offset}"] = e_coli["cai_initial_score"]
        worksheet[f"D{row_offset}"] = e_coli["cai_final_score"]
        worksheet[f"E{row_offset}"] = e_coli["dist_score"]

        worksheet[f"F{row_offset}"] = bacillus["is_wanted"]
        worksheet[f"G{row_offset}"] = bacillus["cai_initial_score"]
        worksheet[f"H{row_offset}"] = bacillus["cai_final_score"]
        worksheet[f"I{row_offset}"] = bacillus["dist_score"]

        worksheet[f"J{row_offset}"] = summary["evaluation"]["average_distance_score"]
        worksheet[f"K{row_offset}"] = summary["evaluation"]["weakest_link_score"]

        orf_summary = summary["orf"]
        aa_to_optimal_codon = orf_summary["aa_to_optimal_codon"]
        for aa, cell in aa_to_cell_column.items():
            worksheet[F"{cell}{row_offset}"] = aa_to_optimal_codon.get(aa) or ""

        worksheet[f"AG{row_offset}"] = orf_summary.get("iterations_count") or 1
        worksheet[f"AH{row_offset}"] = orf_summary["run_time"]


if __name__ == "__main__":
    # generate_summary(results_directory="results")
    generate_summary(results_directory=r"results\CAI_zscore_single_aa_average_ecoli_opt_True_0T45")
