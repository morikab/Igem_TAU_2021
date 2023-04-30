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


def run_from_fasta_file(file_path: str) -> None:
    # with open(file_path, "r") as fasta_handle:
    #     genome_dict = SeqIO.to_dict(SeqIO.parse(fasta_handle, "fasta"))

    with open(r"C:\Users\Kama\Documents\Moran\biomedical-engineering\microbiome-optimization\articles\ORF\ncbi_homo_sapiens_dataset\ncbi_dataset\data\GCF_000001405.40\genomic.gff", "r") as gff_handle:
        gff_dict = SeqIO.to_dict(SeqIO.parse(gff_handle, "gff3", target_seqs=genome_dict))

    len(gff_dict)


def run_all_methods(orf_sequence: typing.Optional[str] = None, orf_sequence_file: typing.Optional[str] = None) -> None:
    for optimization_method in ["single_codon_ratio", "single_codon_diff", "zscore_single_aa_average",
                                "zscore_bulk_aa_average", "zscore_single_aa_weakest_link",
                                "zscore_bulk_aa_weakest_link"]:
        for direction in [True, False]:
            run_single_method_for_orf_sequence(optimization_method=optimization_method,
                                               is_ecoli_optimized=direction,
                                               orf_sequence=orf_sequence,
                                               orf_sequence_file=orf_sequence_file)


def run_single_method_for_orf_sequence(optimization_method: str,
                                       is_ecoli_optimized: bool,
                                       orf_sequence: typing.Optional[str] = None,
                                       orf_sequence_file: typing.Optional[str] = None) -> None:
    tic = time.time()
    default_user_inp_raw = generate_testing_data_for_ecoli_and_bacillus(
        optimization_method=optimization_method,
        optimization_cub_index="CAI",
        clusters_count=1,
        tuning_param=0.5,
        is_ecoli_optimized=is_ecoli_optimized,
        sequence=orf_sequence,
        sequence_file_path=orf_sequence_file,
    )
    run_modules(default_user_inp_raw)
    toc = time.time()
    modules_run_time = toc - tic
    print(F"Modules run time for {optimization_method} is: {modules_run_time}")


if __name__ == "__main__":
    # run_all_methods(orf_sequence=DEFAULT_SEQUENCE)

    # run_single_method_for_orf_sequence(optimization_method="zscore_single_aa_average",
    #                                    is_ecoli_optimized=True,
    #                                    orf_sequence_file=DEFAULT_SEQUENCE_FILE_PATH)

    run_from_fasta_file(file_path=r"C:\Users\Kama\Documents\Moran\biomedical-engineering\microbiome-optimization\articles\ORF\ncbi_homo_sapiens_dataset\ncbi_dataset\data\GCF_000001405.40\GCF_000001405.40_GRCh38.p14_genomic.fna")
