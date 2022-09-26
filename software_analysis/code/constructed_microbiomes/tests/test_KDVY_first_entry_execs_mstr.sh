#!/bin/sh
cd /tamir1/liyamlevi/projects/communique/Igem_TAU_2021/software_analysis/code/constructed_microbiomes/tests/test_local_align_KDVY_first_entry
chmod 777 ./*
for file in ./*; do  qsub -q TullerNano -r y  -e 0_blast_job_error.txt -o 0_blast_job_output.txt -l cput=01:00:00,pmem=1gb,mem=1gb,pvmem=1gb,vmem=1gb "$file"  echo "$file"; done
