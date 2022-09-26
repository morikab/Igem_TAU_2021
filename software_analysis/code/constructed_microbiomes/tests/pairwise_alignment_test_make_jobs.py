import os
from split_genome_16s_fasta import output_fid
from Bio import pairwise2, SeqIO


# python -c 'from cross_tls_with_genome_blast_job.py import filename_to_sent_job; filename_to_sent_job.filename_to_sent_job('vvv')'
def write_job(lines, job_fid):
    f = open(job_fid, 'w')
    f.write(
        '#!/bin/sh \n'
        ' cd /tamir1/liyamlevi/projects/communique/Igem_TAU_2021/software_analysis/code/constructed_microbiomes/tests\n'
        ' module load Python3.9Plus\n'
    )
    for line in lines:
        f.write(line + '\n')
    f.close()



def create_alignment_jobs(seq_fasta, split_genomes_dir=output_fid):
    genomes_files = [os.path.join(split_genomes_dir, i)
                     for i in os.listdir(split_genomes_dir) if '.fasta' in i ]
    for genomes_fid in genomes_files:
        job_name = 'test_dynamic_prog_jobs/' + genomes_fid.split('/')[-1] + '_job.sh'
        command = 'python -c "from pairwise_align_function.py import align_seq_to_fasta; ' + f'pairwise_align_function.align_seq_to_fasta{seq_fasta, genomes_fid}' + '"'
        write_job([command], job_fid = job_name )


metagenome = {name: str(i.seq) for name, i in
              SeqIO.index("../../../data/tested_results/KDVY_example_metagenome/tls.KDVY.1.fsa_nt",
                          "fasta").items()}

seq = metagenome[list(metagenome.keys())[0]]
seq_fasta = "../../../data/tested_results/KDVY_example_metagenome/KDVY_first_entry_KDVY01000001.1.fasta"
create_alignment_jobs(seq_fasta=seq_fasta, split_genomes_dir=output_fid)

