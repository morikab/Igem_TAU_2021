# -*- coding: utf-8 -*-
import os
import subprocess
import pandas as pd
from pathlib import Path
from os import listdir
from os.path import isfile, join

###makes blast jobs for tls (with the 16s sequences as a db
### will run jobs only if there isn't already an existing empty error file for the job

def fastafile_to_blastcsv(fname:str):
    return fname[:-7] +'_' + n_hits + '_hits.csv'

def blastcsv_to_fastafile(fname:str):
    return fname[:-11] + '_' + n_hits + '_hits.csv'

def blastn_run(tls_inp):
    blastn_loc = '/tamir1/liyamlevi/tools/ncbi-blast-2.11.0+/bin/blastn'
    db_loc = '/tamir1/liyamlevi/projects/communique/Igem_TAU_2021/software_analysis/data/processed_genomes/filtered_16s_blastdb/filtered_16s_blastdb'
    other_preferences = ' -max_target_seqs ' + n_hits + ' -outfmt 10 -num_threads 1 -perc_identity 95'
    tls_output = fastafile_to_blastcsv(tls_inp)
    command = blastn_loc + ' -db ' + db_loc + ' -query ' + tls_inp + ' -out ' + tls_output + other_preferences
    print(tls_output)
    print(command)
    return tls_output, command


def run_all_tls(tls_blast_path = '../../data/genbank_tls/'):
    tls_files = [f for f in listdir(tls_blast_path) if isfile(join(tls_blast_path, f))]
    fasta_loc_list = [tls_blast_path + f for f in tls_files if 'fsa_nt' in f]

    commands = []
    outputs = []
    for tls_fid in fasta_loc_list:
        # if Path(tls_fid[:-6]+'csv').is_file(): #instead, skips the run only if the error file exists
        #     print(tls_fid)
        #     continue
        output, command  = blastn_run(tls_fid)
        commands.append(command)
        outputs.append(output)
        # print(output)


    # tls_metadata['blast_command'] = commands
    # tls_metadata['blast_csv'] = outputs
    # tls_metadata.to_csv(metadata_fid[:-4]+ '_with_blast.csv')
    return commands

# run_all_tls(metadata_fid = '../../data/processed_tls/tls_assembly_metadata.csv')


def write_job(lines, job_fid):
    f = open(job_fid, 'w')
    f.write(
        '#!/bin/sh \n cd /tamir1/liyamlevi/projects/communique/Igem_TAU_2021/software_analysis/code/raw_data_analysis\n')
    for line in lines:
        f.write(line + '\n')
    f.close()


def filename_to_sent_job(sh_file, cput = '10:00:00', mem = 3):
    send_prefix = 'qsub -q TullerNano -r y '
    send_suffix = f' -l cput={cput},pmem={mem}gb,mem={mem}gb,pvmem={mem}gb,vmem={mem}gb '
    error_file = sh_file[:-3] + '_error.txt'
    print(error_file)
    output_file = sh_file[:-3] + '_output.txt'
    line = send_prefix + ' -e ' + error_file + ' -o ' + output_file + send_suffix +sh_file
    return line, error_file

def write_mstr_file(job_files, n_hits, output_dir, cput = '10:00:00', mem = 3):

    master_commands = []
    for sh_file in job_files:
        line, error_file = filename_to_sent_job(sh_file, cput, mem)
        error_file = os.path.join(output_dir, error_file)
        if os.path.exists(error_file):
            print(error_file)
            if os.stat(error_file).st_size == 0 : #check that I didn't already successfully run it
                print(error_file, '****')
                continue
        else:
            master_commands.append(line)
    f = open(os.path.join(output_dir, 'mstr_job.sh'), 'w')
    f.write(
        '#!/bin/sh \n cd /tamir1/liyamlevi/projects/communique/Igem_TAU_2021/software_analysis/code/raw_data_analysis/tls_to_16s_blast_' + n_hits + '_hits\n')
    for line in master_commands:
        f.write(line + '\n')
    f.close()

if __name__ == "__main__":

    n_hits= '200'
    print('Start')
    command_list = run_all_tls('../../data/genbank_tls/')
    print(len(command_list))
    job_files = []
    for idx, command in enumerate(command_list):
        filename = str(idx) + '_blast_job.sh'
        job_files.append(filename)
        write_job([command], 'tls_to_16s_blast_' + n_hits + '_hits/' + filename)
    write_mstr_file(job_files, output_dir = f'tls_to_16s_blast_{n_hits}_hits',
                    n_hits = n_hits, cput='10:00:00', mem=3)
    print("The end")






