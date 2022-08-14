# -*- coding: utf-8 -*-
import os
import subprocess
import pandas as pd
from pathlib import Path

destination_dir = '../../data/refseq_genomes/'

def run_cmd(cmd, verbose=False, *args, **kwargs):
    process = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        shell=True
    )
    std_out, std_err = process.communicate()
    if verbose:
        print(std_out.strip(), std_err)



#/tamir1/liyamlevi/tools/ncbi-blast-2.11.0+/bin/blastn -db processed_genomes/blast_db/16s_genomes -query genbank_tls/tls.KDSB.1.fsa_nt -out processed_genomes/out_file3.csv -subject_besthit -outfmt 10 -max_target_seqs 1 -num_threads 4

def blastn_run(tls_inp):
    blastn_loc = '/tamir1/liyamlevi/tools/ncbi-blast-2.11.0+/bin/blastn'
    db_loc = '../../data/processed_genomes/blast_db/16s_genomes'
    other_preferences = ' -subject_besthit -outfmt 10 -max_target_seqs 1 -num_threads 1'
    tls_output = tls_inp[:-6]+'csv'
    command = blastn_loc + ' -db ' + db_loc + ' -query ' + tls_inp + ' -out ' + tls_output + other_preferences
    # run_cmd(command)
    # print(command)
    return tls_output, command


def run_all_tls(metadata_fid ):
    tls_metadata = pd.read_csv(metadata_fid)
    fasta_loc_list = tls_metadata['fasta'].to_list()

    commands = []
    outputs = []
    for tls_fid in fasta_loc_list:
        if Path(tls_fid[:-6]+'csv').is_file():
            continue
        output, command  = blastn_run(tls_fid)
        commands.append(command)
        outputs.append(output)
        print(output)

    # tls_metadata['blast_command'] = commands
    # tls_metadata['blast_csv'] = outputs
    # tls_metadata.to_csv(metadata_fid[:-4]+ '_with_blast.csv')
    return commands, tls_metadata

# run_all_tls(metadata_fid = '../../data/processed_tls/tls_assembly_metadata.csv')


def write_job(lines, job_fid):
    f = open(job_fid, 'w')
    f.write(
        '#!/bin/sh \n cd /tamir1/liyamlevi/projects/communique/Igem_TAU_2021/software_analysis/code/raw_data_analysis\n')
    for line in lines:
        f.write(line + '\n')
    f.close()


def filename_to_sent_job(sh_file):
    send_prefix = 'qsub -q TullerNano -r y '
    send_suffix = ' -l cput=01:00:00,pmem=1gb,mem=1gb,pvmem=1gb,vmem=1gb '
    error_file = sh_file[:-3] + '_error.txt'
    output_file = sh_file[:-3] + '_output.txt'
    line = send_prefix + ' -e ' + error_file + ' -o ' + output_file + send_suffix +sh_file
    print(line)
    return line


if __name__ == "__main__":
    print('Start')
    command_list, df = run_all_tls(metadata_fid = '../../data/processed_tls/tls_assembly_metadata.csv')
    job_files = []
    for idx, command in enumerate(command_list):
        filename = str(idx) + '_blast_job.sh'
        job_files.append(filename)
        write_job([command], 'jobs/' + filename)


    master_commands = [filename_to_sent_job(sh_file) for sh_file in job_files]
    f = open('jobs/mstr_job.sh', 'w')
    f.write(
        '#!/bin/sh \n')
    for line in master_commands:
        f.write(line + '\n')
    f.close()
        # proc = subprocess.Popen(command, shell=True,
        #                         stdin=subprocess.PIPE, stdout=subprocess.PIPE)
        # out, err = proc.communicate()
        # print("{} : {}".format(command, out.decode()))
    print("The end")






