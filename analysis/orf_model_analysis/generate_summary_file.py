import argparse
import json
import os
import typing
from pathlib import Path

import openpyxl


aa_list = ["C", "D", "S", "Q", "M", "N", "P", "K", "_", "T", "F", "A", "G", "I", "L", "H", "R", "W", "V", "E", "Y"]


def accumulate_summary_files(results_directory: str, requested_optimization_method: str) -> None:
    filename = "summary.xlsx"

    with open("description_to_gene_mapping.json", "r") as genes_file:
        description_to_gene_mapping = json.load(genes_file)

    columns_per_organism_count = 4
    general_columns_count = 4 + len(aa_list)
    for root, dirs, files in os.walk(results_directory):
        for file in files:
            if file == filename:
                file_path = os.path.join(root, file)
                sequence_name = Path(root).name.replace("-", "|")
                gene_name = description_to_gene_mapping.get(sequence_name)
                workbook = openpyxl.load_workbook(filename=file_path)
                worksheet = workbook.active

                original_headers = [cell.value for cell in worksheet[1]]
                headers = ["Sequence", "Gene"] + original_headers
                optimization_method_to_scores = {}
                for row in worksheet.iter_rows(min_row=2, values_only=True):
                    if row[0] != requested_optimization_method:
                        continue
                    # optimization_method = models.OptimizationMethod(value=row[1])

                    organisms_optimization_index = "_".join(
                        F"{original_headers[i]}_{str(row[i])}" for i in
                        range(1, len(row) - general_columns_count, columns_per_organism_count)
                    )
                    key = "-".join([row[0], organisms_optimization_index])

                    # Needed to remove duplicates from the summary file
                    optimization_method_to_scores[key] = (sequence_name, gene_name) + row

                for index, row in optimization_method_to_scores.items():
                    summary_file_path = F"{index}-summary.xlsx"
                    if not os.path.exists(summary_file_path):
                        summary_workbook = openpyxl.Workbook()
                        summary_worksheet = summary_workbook.active
                        summary_worksheet.append(headers)
                        summary_worksheet.append(row)
                        summary_workbook.save(summary_file_path)
                    else:
                        summary_workbook = openpyxl.load_workbook(summary_file_path)
                        summary_worksheet = summary_workbook.active
                        summary_worksheet.append(row)
                        summary_workbook.save(summary_file_path)


def generate_summary(results_directory: str) -> None:
    filename = "run_summary.json"

    workbook = openpyxl.Workbook()
    worksheet = workbook.active

    for root, dirs, files in os.walk(results_directory):
        for file in files:
            if file == filename:
                file_path = os.path.join(root, file)
                update_from_summary(file_path=file_path, worksheet=worksheet)

    workbook.save(os.path.join(results_directory, "summary.xlsx"))


def add_cell_with_value(worksheet, row: int, column: int, value: typing.Any) -> None:
    worksheet.insert_cols(column)
    worksheet.cell(row=row, column=column, value=value)


def initialize_column_headers(summary: typing.Dict, worksheet) -> None:
    add_cell_with_value(worksheet=worksheet, row=1, column=worksheet.max_column, value="optimization_method")

    organisms = summary["user_input"]["organisms"]
    # Display unwanted organisms first
    organisms.sort(key=lambda x: x.get("is_wanted"))
    for organism in organisms:
        organism_name = organism["name"]
        formatted_organism_name = "_".join(organism_name.split(" ")).lower()
        add_cell_with_value(worksheet=worksheet,
                            row=1,
                            column=worksheet.max_column,
                            value=F"{formatted_organism_name}_is_wanted")
        add_cell_with_value(worksheet=worksheet,
                            row=1,
                            column=worksheet.max_column,
                            value=F"{formatted_organism_name}_initial_gene_cai_score")
        add_cell_with_value(worksheet=worksheet,
                            row=1,
                            column=worksheet.max_column,
                            value=F"{formatted_organism_name}_final_gene_cai_score")
        add_cell_with_value(worksheet=worksheet,
                            row=1,
                            column=worksheet.max_column,
                            value=F"{formatted_organism_name}_dist_score")

    add_cell_with_value(worksheet=worksheet, row=1, column=worksheet.max_column, value="final_average_distance_score")
    add_cell_with_value(worksheet=worksheet, row=1, column=worksheet.max_column, value="final_weakest_link_score")

    for aa in aa_list:
        add_cell_with_value(worksheet=worksheet, row=1, column=worksheet.max_column, value=aa)

    add_cell_with_value(worksheet=worksheet, row=1, column=worksheet.max_column, value="number_of_iterations")
    add_cell_with_value(worksheet=worksheet, row=1, column=worksheet.max_column, value="run_time (seconds)")


def update_from_summary(file_path: str, worksheet) -> None:
    with open(file_path, "r") as summary_file:
        summary = json.load(summary_file)

    if worksheet.max_row == 1:
        initialize_column_headers(summary=summary, worksheet=worksheet)

    optimization_method = summary["user_input"]["optimization_method"]
    summary_row = [optimization_method]

    organisms = summary["evaluation"]["organisms"]
    organisms.sort(key=lambda x: x.get("is_wanted"))
    for organism in organisms:
        summary_row.append(organism["is_wanted"])
        summary_row.append(organism["cai_initial_score"])
        summary_row.append(organism["cai_final_score"])
        summary_row.append(organism["dist_score"])

    summary_row.append(summary["evaluation"]["average_distance_score"])
    summary_row.append(summary["evaluation"]["weakest_link_score"])

    orf_summary = summary["orf"]
    aa_to_optimal_codon = orf_summary["aa_to_optimal_codon"]
    for aa in aa_list:
        summary_row.append(aa_to_optimal_codon.get(aa) or "")

    summary_row.append(orf_summary.get("iterations_count") or 1)
    summary_row.append(orf_summary["run_time"])

    worksheet.append(summary_row)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Summary script parser")
    parser.add_argument('-m', '--method', type=str, required=True, help="Optimization method to collect results for.")

    args = parser.parse_args()
    requested_optimization_method = args.method
    # 1
    # generate_summary(results_directory=r"results\CAI_zscore_single_aa_average_ecoli_opt_True_0T45")

    # 2
    root_dir = r"C:\projects\Igem_TAU_2021_moran\analysis\orf_model_analysis\results_human"
    # for file_name in os.listdir(root_dir):
    #     generate_summary(results_directory=os.path.join(root_dir, file_name))

    # 3
    accumulate_summary_files(results_directory=root_dir, requested_optimization_method=requested_optimization_method)
