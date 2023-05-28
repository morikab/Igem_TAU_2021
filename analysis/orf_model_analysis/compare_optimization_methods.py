import argparse
import json
import os
import re
import typing
from collections import defaultdict
from pathlib import Path

from Bio import SeqIO

from input_testing_data.generate_input_testing_data_for_modules import \
    generate_testing_data_for_ecoli_and_bacillus
from modules.main import run_modules


current_directory = Path(__file__).parent.resolve()
base_path = os.path.join(Path(current_directory).parent.resolve(), "example_data")
DEFAULT_SEQUENCE_FILE_PATH = os.path.join(base_path, "mCherry_original.fasta")


def run_from_fasta_file(fasta_file_path: str,
                        records_file_path: str,
                        start_record: str,
                        max_records_count: int) -> None:
    with open(fasta_file_path, "r") as fasta_handle:
        genome_dict = SeqIO.to_dict(SeqIO.parse(fasta_handle, "fasta"))

    with open(records_file_path, "r") as records_file:
        records = json.load(records_file)

    is_record_found = False
    count = 0
    for record in records.values():
        if record == start_record:
            is_record_found = True
        if not is_record_found:
            continue
        if count >= max_records_count:
            break

        count += 1
        value = genome_dict[record]
        run_all_methods(orf_sequence=str(value.seq),
                        output_path=record.replace('|', '-'))


def run_all_methods(orf_sequence: typing.Optional[str] = None,
                    orf_sequence_file: typing.Optional[str] = None,
                    output_path: typing.Optional[str] = None) -> None:
    # TODO - add a new optimization method here..
    for optimization_method in [
        # "single_codon_ratio", "single_codon_diff",
        "zscore_single_aa_average", # "zscore_bulk_aa_average",
        # "zscore_single_aa_weakest_link", "zscore_bulk_aa_weakest_link",
    ]:
        for direction in [True, False]:
            run_single_method_for_orf_sequence(optimization_method=optimization_method,
                                               is_ecoli_optimized=direction,
                                               orf_sequence=orf_sequence,
                                               orf_sequence_file=orf_sequence_file,
                                               output_path=output_path)


def extract_sequences_for_analysis(fasta_file_path: str) -> None:
    with open(fasta_file_path, "r") as fasta_handle:
        genome_dict = SeqIO.to_dict(SeqIO.parse(fasta_handle, "fasta"))

    missing_genes = []
    gene_mapping = defaultdict(list)
    description_to_gene_mapping = {}

    for record, value in genome_dict.items():
        if "NC" not in record or "XP" in record:
            continue
        parameters = re.findall("gene=.*]", value.description)
        if not parameters:
            missing_genes.append(record)
            continue
        gene_parameter = parameters[0].split("]")[0]
        gene_name = gene_parameter.strip("gene=")

        if len(value.seq) % 3 != 0:
            print(F"Skip sequence {record} for {gene_name} because of length not divisible by 3.")
            continue

        gene_mapping[gene_name].append(record)
        description_to_gene_mapping[record] = gene_name

    with open("gene_mapping.json", "w") as genes_file:
        json.dump(gene_mapping, genes_file)
    with open("description_to_gene_mapping.json", "w") as description_to_genes_file:
        json.dump(description_to_gene_mapping, description_to_genes_file)
    with open("missing_genes.txt", "w") as missing_genes_file:
        for gene in missing_genes:
            missing_genes_file.write(gene)

    gene_to_longest_sequence = {
        key: max(value, key=lambda x: len(genome_dict[x].seq)) for key, value in gene_mapping.items()
    }
    with open("gene_to_longest_sequence.json", "w") as genes_file:
        json.dump(gene_to_longest_sequence, genes_file)


def run_single_method_for_orf_sequence(optimization_method: str,
                                       is_ecoli_optimized: bool,
                                       orf_sequence: typing.Optional[str] = None,
                                       orf_sequence_file: typing.Optional[str] = None,
                                       output_path: typing.Optional[str] = None) -> None:
    default_user_inp_raw = generate_testing_data_for_ecoli_and_bacillus(
        optimization_method=optimization_method,
        optimization_cub_index="CAI",
        clusters_count=1,
        tuning_param=0.5,
        is_ecoli_optimized=is_ecoli_optimized,
        sequence=orf_sequence,
        sequence_file_path=orf_sequence_file,
        output_path=os.path.join("results_human", output_path),
    )
    run_modules(default_user_inp_raw)


def compare_gene_mappings() -> None:
    with open("gene_to_longest_sequence_all.json") as gene_mapping_file:
        previous_mapping = json.load(gene_mapping_file)

    with open("gene_to_longest_sequence.json") as gene_mapping_file:
        current_mapping = json.load(gene_mapping_file)

    changed_genes = []
    uncovered_genes = []
    for gene in previous_mapping.keys():
        if current_mapping.get(gene) is None:
            uncovered_genes.append(gene)
            continue

        if previous_mapping[gene] != current_mapping[gene]:
            changed_genes.append(current_mapping[gene])

    with open("changed_genes.txt", "w") as changed_genes_file:
        for gene in changed_genes:
            changed_genes_file.write(current_mapping[gene]+"\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Analysis script parser")
    parser.add_argument('-s', '--start', type=str, required=True, help="Fasta record description to start running from")
    parser.add_argument('-n', '--number', type=int, help="Number of records to parse from the given start record")

    args = parser.parse_args()

    # Reference - https://www.ncbi.nlm.nih.gov/data-hub/genome/GCF_000001405.40/

    # extract_sequences_for_analysis(
    #     fasta_file_path=r"C:\Users\Kama\Documents\Moran\biomedical-engineering\microbiome-optimization\articles\ORF\ncbi_homo_sapiens_dataset\ncbi_dataset\data\GCF_000001405.40\cds_from_genomic.fna")

    run_from_fasta_file(
        fasta_file_path=r"C:\Users\Kama\Documents\Moran\biomedical-engineering\microbiome-optimization\articles\ORF\ncbi_homo_sapiens_dataset\ncbi_dataset\data\GCF_000001405.40\cds_from_genomic.fna",
        records_file_path="gene_to_longest_sequence.json",
        start_record=args.start,
        max_records_count=args.number or 500,
    )
