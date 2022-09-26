from Bio import pairwise2, SeqIO
import os
import pandas as pd
from split_genome_16s_fasta import output_fid


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



def align_seq_to_fasta(seq_fid, genomes_fid):
    # print(genomes_fid)
    seq = [str(i.seq) for name, i in
                 SeqIO.index(seq_fid,
                             "fasta").items()][0]

    genomes = {name: str(i.seq) for name, i in
                 SeqIO.index(genomes_fid,
                             "fasta").items()}
    scores = {name:pairwise2.align.localxx(seq, rrna_seq)[0].score
              for name, rrna_seq in genomes.items()}
    scores_df = pd.DataFrame.from_dict(scores)
    print(scores_df)
    scores_df.to_csv(genomes_fid[:6] + '_dynamic_programing.csv')



def create_alignment_jobs(seq_fasta, split_genomes_dir=output_fid):
    genomes_files = [os.path.join(split_genomes_dir, i)
                     for i in os.listdir(split_genomes_dir) if '.fasta' in i ]
    print(genomes_files)
    for genomes_fid in genomes_files:
        job_name = 'test_dynamic_prog_jobs/' + genomes_fid.split('/')[-1] + '_job.sh'
        print(job_name)
        command = f"python -c 'from pairwise_alignment_test import align_seq_to_fasta;  pairwise_alignment_test.align_seq_to_fasta{seq_fasta, genomes_fid}'"
        write_job([command], job_fid = job_name )
        print(command)


metagenome = {name: str(i.seq) for name, i in
              SeqIO.index("../../../data/tested_results/KDVY_example_metagenome/tls.KDVY.1.fsa_nt",
                          "fasta").items()}

seq = metagenome[list(metagenome.keys())[0]]
seq_fasta = "../../../data/tested_results/KDVY_example_metagenome/KDVY_first_entry_KDVY01000001.1.fasta"
create_alignment_jobs(seq_fasta=seq_fasta, split_genomes_dir=output_fid)

