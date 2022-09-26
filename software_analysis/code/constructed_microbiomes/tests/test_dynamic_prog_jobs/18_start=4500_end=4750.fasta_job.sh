#!/bin/sh 
 cd /tamir1/liyamlevi/projects/communique/Igem_TAU_2021/software_analysis/code/constructed_microbiomes/tests
 module load Python3.9Plus
python -c 'from pairwise_alignment_test import align_seq_to_fasta;  pairwise_alignment_test.align_seq_to_fasta('../../../data/tested_results/KDVY_example_metagenome/KDVY_first_entry_KDVY01000001.1.fasta', '../../../data/tested_results/genomes_16s_sliced/18_start=4500_end=4750.fasta')'