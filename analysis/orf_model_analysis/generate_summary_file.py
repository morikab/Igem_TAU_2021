import json
import os
import typing

import openpyxl


aa_list = ["C", "D", "S", "Q", "M", "N", "P", "K", "_", "T", "F", "A", "G", "I", "L", "H", "R", "W", "V", "E", "Y"]


def generate_summary(results_directory: str) -> None:
    filename = "run_summary.json"

    workbook = openpyxl.Workbook()
    worksheet = workbook.active

    for root, dirs, files in os.walk(results_directory):
        for file in files:
            if file == filename:
                file_path = os.path.join(root, file)
                update_from_summary(file_path=file_path, worksheet=worksheet)

    workbook.save("summary.xlsx")


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
    # generate_summary(results_directory="results")
    generate_summary(results_directory=r"results\CAI_zscore_single_aa_average_ecoli_opt_True_0T45")
