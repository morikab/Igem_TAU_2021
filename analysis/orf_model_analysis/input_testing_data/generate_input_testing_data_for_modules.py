import os
import random
import string
from pathlib import Path
from os import listdir
from os.path import isfile, join


current_directory = Path(__file__).parent.resolve()
base_path = os.path.join(Path(current_directory).parent.parent.resolve(), "example_data")
default_genome_path = os.path.join(base_path, 'arabidopsis_microbiome')

DEFAULT_ORGANISM_PRIORITY = 50
DEFAULT_CLUSTERS_COUNT = 2
DEFAULT_TUNING_PARAM = 0.5
DEFAULT_SEQUENCE_FILE_PATH = os.path.join(base_path, "mCherry_original.fasta")


def generate_random_string(length: int) -> str:
    letters_and_digits = string.ascii_letters + string.digits
    return "".join(random.choice(letters_and_digits) for _ in range(length))


def generate_testing_data(n_organisms=15,
                          percent_optimized=0.5,
                          clusters_count=DEFAULT_CLUSTERS_COUNT,
                          tuning_param=DEFAULT_TUNING_PARAM,
                          genome_path=default_genome_path):

    inp_dict = {
            'sequence': os.path.join(base_path, 'mCherry_original.fasta'),
            'tuning_param': tuning_param,
            'organisms': {},
            'clusters_count': clusters_count,
    }
    genome_list = [f for f in listdir(genome_path) if isfile(join(genome_path, f))]
    if n_organisms > len(genome_list):
        raise ValueError('not enough genomes in data. select less genomes')

    opt_genomes = random.sample(genome_list,
                                round(n_organisms*percent_optimized))

    for opt_genome in opt_genomes:
        inp_dict['organisms'][opt_genome[:-3]] = {
            'genome_path': os.path.join(genome_path, opt_genome),
            'optimized': True,
            'expression_csv': None,
            'optimization_priority': DEFAULT_ORGANISM_PRIORITY,
        }

    deopt_genomes = random.sample([i for i in genome_list if i not in opt_genomes],
                                   round(n_organisms*(1-percent_optimized)))

    for deopt_genome in deopt_genomes:
        inp_dict['organisms'][deopt_genome[:-2]] = {
            'genome_path': os.path.join(genome_path, deopt_genome),
            'optimized': False,
            'expression_csv': None,
            'optimization_priority': DEFAULT_ORGANISM_PRIORITY,
        }

    return inp_dict


def generate_testing_data_for_ecoli_and_bacillus(
        optimization_method: str,
        optimization_cub_index: str,
        clusters_count: int = DEFAULT_CLUSTERS_COUNT,
        tuning_param: float = DEFAULT_TUNING_PARAM,
        is_ecoli_optimized: bool = False,
        sequence_file_path: str = None,
        sequence: str = None,
        output_path: str = None,
):
    assert (sequence is not None or sequence_file_path is not None), \
        "Should provide either a sequence or a sequence file path"

    output_path = output_path or "results"
    output_directory = os.path.join(output_path, F"{optimization_cub_index}_{optimization_method}_ecoli_opt_"
                                                 F"{is_ecoli_optimized}_{generate_random_string(4)}")
    Path(output_directory).mkdir(parents=True, exist_ok=True)

    if is_ecoli_optimized:
        opt_genome = "Escherichia coli.gb"
        opt_mrna_levels = "ecoli_mrna_level.csv"
        deopt_genome = "Bacillus subtilis.gb"
        deopt_mrna_levels = "bacillus_mrna_level.csv"
    else:
        opt_genome = "Bacillus subtilis.gb"
        opt_mrna_levels = "bacillus_mrna_level.csv"
        deopt_genome = "Escherichia coli.gb"
        deopt_mrna_levels = "ecoli_mrna_level.csv"

    inp_dict = {
        "sequence_file_path": sequence_file_path,
        "sequence": sequence,
        "tuning_param": tuning_param,
        "organisms": {},
        "clusters_count": clusters_count,
        "optimization_method": optimization_method,
        "optimization_cub_index": optimization_cub_index,
        "output_path": output_directory,
    }

    inp_dict['organisms'][opt_genome[:-3]] = {
        "genome_path": os.path.join(base_path, opt_genome),
        "optimized": True,
        "expression_csv": os.path.join(base_path, opt_mrna_levels),
        "optimization_priority": DEFAULT_ORGANISM_PRIORITY,
    }

    inp_dict["organisms"][deopt_genome[:-2]] = {
        "genome_path": os.path.join(base_path, deopt_genome),
        "optimized": False,
        "expression_csv": os.path.join(base_path, deopt_mrna_levels),
        "optimization_priority": DEFAULT_ORGANISM_PRIORITY,
    }

    return inp_dict
