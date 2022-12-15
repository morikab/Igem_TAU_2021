#!/bin/sh 
 cd /tamir1/liyamlevi/projects/communique/Igem_TAU_2021/software_analysis/code/raw_data_analysis/tls_to_16s_blast_5_hits
qsub -q TullerNano -r y  -e 0_blast_job_error.txt -o 0_blast_job_output.txt -l cput=08:00:00,pmem=1gb,mem=1gb,pvmem=1gb,vmem=1gb 0_blast_job.sh
qsub -q TullerNano -r y  -e 1_blast_job_error.txt -o 1_blast_job_output.txt -l cput=08:00:00,pmem=1gb,mem=1gb,pvmem=1gb,vmem=1gb 1_blast_job.sh
qsub -q TullerNano -r y  -e 2_blast_job_error.txt -o 2_blast_job_output.txt -l cput=08:00:00,pmem=1gb,mem=1gb,pvmem=1gb,vmem=1gb 2_blast_job.sh
qsub -q TullerNano -r y  -e 3_blast_job_error.txt -o 3_blast_job_output.txt -l cput=08:00:00,pmem=1gb,mem=1gb,pvmem=1gb,vmem=1gb 3_blast_job.sh
qsub -q TullerNano -r y  -e 4_blast_job_error.txt -o 4_blast_job_output.txt -l cput=08:00:00,pmem=1gb,mem=1gb,pvmem=1gb,vmem=1gb 4_blast_job.sh
qsub -q TullerNano -r y  -e 5_blast_job_error.txt -o 5_blast_job_output.txt -l cput=08:00:00,pmem=1gb,mem=1gb,pvmem=1gb,vmem=1gb 5_blast_job.sh
