# -*- coding: utf-8 -*-
import os
import subprocess
import pandas as pd

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
    other_preferences = ' -subject_besthit -outfmt 10 -max_target_seqs 1 -num_threads 4'
    tls_output = tls_inp[:-6]+'csv'
    command = blastn_loc + ' -db ' + db_loc + ' -query ' + tls_inp + ' -out ' + tls_output + other_preferences
    run_cmd(command)
    print(command)
    return tls_output, command


def run_all_tls(metadata_fid ):
    tls_metadata = pd.read_csv(metadata_fid)
    fasta_loc_list = tls_metadata['fasta'].to_list()
    print(fasta_loc_list)

    commands = []
    outputs = []
    for tls_fid in fasta_loc_list:
        command, output = blastn_run(tls_fid)
        commands.append(command)
        outputs.append(output)

    tls_metadata['blast_command'] = commands
    tls_metadata['blast_csv'] = outputs
    tls_metadata.to_csv(metadata_fid[:-4]+ '_with_blast.csv')

run_all_tls(metadata_fid = '../../data/processed_tls/tls_assembly_metadata.csv')





