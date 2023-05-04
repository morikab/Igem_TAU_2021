import argparse
import os
import time
import typing
from pathlib import Path

from Bio import SeqIO

from input_testing_data.generate_input_testing_data_for_modules import \
    generate_testing_data_for_ecoli_and_bacillus
from modules.main import run_modules


current_directory = Path(__file__).parent.resolve()
base_path = os.path.join(Path(current_directory).parent.resolve(), "example_data")
DEFAULT_SEQUENCE_FILE_PATH = os.path.join(base_path, "mCherry_original.fasta")


def run_from_fasta_file(file_path: str, start_record: str, max_records_count: int) -> None:
    with open(file_path, "r") as fasta_handle:
        genome_dict = SeqIO.to_dict(SeqIO.parse(fasta_handle, "fasta"))

    # nc_records = []
    # for key in genome_dict.keys():
    #     splitted = key.split("_")
    #     splitted = [x.strip("lcl|") for x in splitted]
    #     # Consider only reference sequences that were manually validated.
    #     if "NC" in splitted and "XP" not in splitted:
    #         nc_records.append(key)
    #
    # with open("records", "w") as records_file:
    #     for record in nc_records:
    #         records_file.write(record + "\n")

    with open("all_records") as records_file:
        records = records_file.readlines()

    is_record_found = False
    count = 0

    for record in records:
        if record == start_record:
            is_record_found = True
        if not is_record_found:
            continue
        if count >= max_records_count:
            break

        count += 1
        value = genome_dict[record]
        run_all_methods(orf_sequence=str(value.seq),
                        output_path=os.path.join("results_human", record.replace('|', '-')))

    # for key, value in genome_dict.items():
    #     if key == start_record:
    #         is_record_found = True
    #     if not is_record_found:
    #         continue
    #     if count >= max_records_count:
    #         break
    #
    #     count += 1
    #     run_all_methods(orf_sequence=str(value.seq),
    #                     output_path=os.path.join("results_human", key.replace('|', '-')))


def run_all_methods(orf_sequence: typing.Optional[str] = None,
                    orf_sequence_file: typing.Optional[str] = None,
                    output_path: typing.Optional[str] = None) -> None:
    for optimization_method in [
        "single_codon_ratio", # "single_codon_diff",
        # "zscore_single_aa_average", "zscore_bulk_aa_average",
        # "zscore_single_aa_weakest_link", "zscore_bulk_aa_weakest_link",
    ]:
        for direction in [True, False]:
            run_single_method_for_orf_sequence(optimization_method=optimization_method,
                                               is_ecoli_optimized=direction,
                                               orf_sequence=orf_sequence,
                                               orf_sequence_file=orf_sequence_file,
                                               output_path=output_path)


def run_single_method_for_orf_sequence(optimization_method: str,
                                       is_ecoli_optimized: bool,
                                       orf_sequence: typing.Optional[str] = None,
                                       orf_sequence_file: typing.Optional[str] = None,
                                       output_path: typing.Optional[str] = None) -> None:
    tic = time.time()
    default_user_inp_raw = generate_testing_data_for_ecoli_and_bacillus(
        optimization_method=optimization_method,
        optimization_cub_index="CAI",
        clusters_count=1,
        tuning_param=0.5,
        is_ecoli_optimized=is_ecoli_optimized,
        sequence=orf_sequence,
        sequence_file_path=orf_sequence_file,
        output_path=output_path,
    )
    run_modules(default_user_inp_raw)
    toc = time.time()
    modules_run_time = toc - tic
    print(F"Modules run time for {optimization_method} is: {modules_run_time}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Analysis script parser")
    parser.add_argument('-s', '--start', type=str, required=True, help="Fasta record description to start running from")
    parser.add_argument('-n', '--number', type=int, help="Number of records to parse from the given start record")

    args = parser.parse_args()
    # run_all_methods(orf_sequence=DEFAULT_SEQUENCE)

    # run_single_method_for_orf_sequence(optimization_method="zscore_single_aa_average",
    #                                    is_ecoli_optimized=True,
    #                                    orf_sequence_file=DEFAULT_SEQUENCE_FILE_PATH)

    # 145289 records in dict
    # Reference - https://www.ncbi.nlm.nih.gov/data-hub/genome/GCF_000001405.40/
    run_from_fasta_file(
        file_path=r"C:\Users\Kama\Documents\Moran\biomedical-engineering\microbiome-optimization\articles\ORF\ncbi_homo_sapiens_dataset\ncbi_dataset\data\GCF_000001405.40\cds_from_genomic.fna",
        start_record=args.start,
        max_records_count=args.number or 500,
    )

    # 121766 records in dict
    # run_from_fasta_file(
    #     file_path=r"C:\Users\Kama\Documents\Moran\biomedical-engineering\microbiome-optimization\articles\ORF\Homo_sapiens.GRCh38.cds.all.fa\Homo_sapiens.GRCh38.cds.all.fa")