import os
from split_genome_16s_fasta import output_fid
from Bio import pairwise2, SeqIO


# python -c 'from cross_tls_with_genome_blast_job.py import filename_to_sent_job; filename_to_sent_job.filename_to_sent_job('vvv')'
def write_job(seq_fasta, genomes_fid, job_fid):
    f = open(job_fid, 'w')

    out_dir = '../../../../data/tested_results/KDVY_example_metagenome/sliced_alignment_first_entry/'
    f.write(
        '#!/powerapps/share/centos7/miniconda/miniconda3-4.7.12-environmentally/envs/Python3.9Plus/bin/python\n'
        'import sys\n'
        "sys.path.insert(1, '../')\n"
        'from pairwise_align_function import align_seq_to_fasta\n'
        f"METAGENOME_SEQ = '../{seq_fasta}'\n"
        f"GENOME_SLICE = '../{genomes_fid}'\n"
        f"OUT_DIR = '../{out_dir}'\n"
        'align_seq_to_fasta(METAGENOME_SEQ, GENOME_SLICE, OUT_DIR)\n'
    )
    f.close()



def create_alignment_jobs(seq_fasta, split_genomes_dir=output_fid):
    genomes_files = [os.path.join(split_genomes_dir, i)
                     for i in os.listdir(split_genomes_dir) if '.fasta' in i ]
    for genomes_fid in genomes_files:
        job_name = 'test_local_align_KDVY_first_entry/' + genomes_fid.split('/')[-1] + '_exec.py'
        write_job(seq_fasta, genomes_fid, job_fid = job_name)




metagenome = {name: str(i.seq) for name, i in
              SeqIO.index("../../../data/tested_results/KDVY_example_metagenome/tls.KDVY.1.fsa_nt",
                          "fasta").items()}

seq = metagenome[list(metagenome.keys())[0]]
seq_fasta = "../../../data/tested_results/KDVY_example_metagenome/KDVY_first_entry_KDVY01000001.1.fasta"
create_alignment_jobs(seq_fasta=seq_fasta, split_genomes_dir=output_fid)

