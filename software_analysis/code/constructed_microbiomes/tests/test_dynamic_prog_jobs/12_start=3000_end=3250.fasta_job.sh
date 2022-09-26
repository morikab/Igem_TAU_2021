#!/bin/sh 
 cd /tamir1/liyamlevi/projects/communique/Igem_TAU_2021/software_analysis/code/constructed_microbiomes/tests
 module load Python3.9Plus
python -c "from pairwise_align_function.py import align_seq_to_fasta; pairwise_align_function.align_seq_to_fasta('../../../data/tested_results/KDVY_example_metagenome/KDVY_first_entry_KDVY01000001.1.fasta', '../../../data/tested_results/genomes_16s_sliced/12_start=3000_end=3250.fasta')"
