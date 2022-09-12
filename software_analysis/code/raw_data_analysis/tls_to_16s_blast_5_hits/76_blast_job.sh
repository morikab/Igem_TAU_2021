#!/bin/sh 
 cd /tamir1/liyamlevi/projects/communique/Igem_TAU_2021/software_analysis/code/raw_data_analysis
/tamir1/liyamlevi/tools/ncbi-blast-2.11.0+/bin/blastn -db /tamir1/liyamlevi/projects/communique/Igem_TAU_2021/software_analysis/data/processed_genomes/filtered_16s_blastdb/filtered_16s_blastdb -query ../../data/genbank_tls/tls.KBPM.1.fsa_nt -out ../../data/genbank_tls/tls.KBPM.1_5_hits.csv - -max_target_seqs 5 -outfmt 10 -num_threads 1
