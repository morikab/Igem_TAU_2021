import os
from split_genome_16s_fasta import output_fid
from Bio import pairwise2, SeqIO


# python -c 'from cross_tls_with_genome_blast_job.py import filename_to_sent_job; filename_to_sent_job.filename_to_sent_job('vvv')'
def write_job(seq_fasta, genomes_fid, job_fid):
    exec_fname = job_fid + '_exec.py'
    exec_file = open(exec_fname, 'w')
    print('./' + exec_fname.split('/')[-1], ' >>  ', exec_fname.split('/')[-1][:-2], 'txt')
    out_dir = '../../../data/tested_results/KDVY_example_metagenome/sliced_alignment_first_entry/'
    exec_file.write(
        '#!/powerapps/share/centos7/miniconda/miniconda3-4.7.12-environmentally/envs/Python3.9Plus/bin/python\n'
        'import sys\n'
        "sys.path.insert(1, '../')\n"
        'from pairwise_align_function import align_seq_to_fasta\n'
        f"METAGENOME_SEQ = '../{seq_fasta}'\n"
        f"GENOME_SLICE = '../{genomes_fid}'\n"
        f"OUT_DIR = '../{out_dir}'\n"
        'align_seq_to_fasta(METAGENOME_SEQ, GENOME_SLICE, OUT_DIR)\n'
    )

    exec_file.close()

    sh_file = open(job_fid + '_job.sh', 'w')
    sh_file.write(
        '# !/bin/sh\n'
        'cd /tamir1/liyamlevi/projects/communique/Igem_TAU_2021/software_analysis/code/constructed_microbiomes/tests/test_local_align_KDVY_first_entry\n'
        'module load Python3.9Plus\n'
        f'python {exec_fname}\n'
    )
    sh_file.close()

def create_alignment_jobs(seq_fasta, split_genomes_dir=output_fid):
    genomes_files = [os.path.join(split_genomes_dir, i)
                     for i in os.listdir(split_genomes_dir) if '.fasta' in i ]

    job_names  = []
    for genomes_fid in genomes_files:
        job_name = 'test_local_align_KDVY_first_entry/' + genomes_fid.split('/')[-1][:-6]
        write_job(seq_fasta, genomes_fid, job_fid = job_name)
        job_names.append(job_name)


    return job_names

def make_mstr_job(mstr_fid, job_names):
    f = open(mstr_fid, 'w')

    f.write(
        '# !/bin/sh\n '
        ' cd /tamir1/liyamlevi/projects/communique/Igem_TAU_2021/software_analysis/code/constructed_microbiomes/tests/test_local_align_KDVY_first_entry/\n'
        'chmod 777 ./*\n'
    )

    for idx, job in enumerate(job_names):
        job = job.split('/')[-1] + '_job.sh'
        f.write(f'qsub -q TullerNano -r y  -e {idx}_job_error.txt -o {idx}_output.txt -l cput=01:00:00,pmem=1gb,mem=1gb,pvmem=1gb,vmem=1gb {job}\n')

    f.close()





seq_fasta = "../../../data/tested_results/KDVY_example_metagenome/KDVY_first_entry_KDVY01000001.1.fasta"
job_names = create_alignment_jobs(seq_fasta=seq_fasta, split_genomes_dir=output_fid)
print(job_names)
make_mstr_job('local_alignment_fist_entry_mstr.sh', job_names)
