# -*- coding: utf-8 -*-
import os
import subprocess
import pandas as pd
from pathlib import Path
import json
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




def blastn_run(tls_csv, genomes_df):
    tls_df = pd.read_csv(tls_csv)
    found_genomes = [ i[1] for i in tls_df. iloc[:, 0:2].values.tolist()]
    unique_found_genomes = sorted(set(found_genomes))
    genomes_df[tls_csv] = [i in unique_found_genomes for i in range(len(genomes_df))]
    cai_weights = genomes_df.iloc[unique_found_genomes, list(range(4,68))].transpose().to_dict()
    return cai_weights, genomes_df


def check_all_blast_res(genomes_df, tls_new_metadata_df):
    tls_new_metadata_df = tls_new_metadata_df.iloc[:, 1:]
    print(tls_new_metadata_df)

    blast_results_dict = {}
    for idx, blast_csv in enumerate(tls_new_metadata_df['blast_csv'].to_list()):
        # if idx>2:
        #     break
        tls_dict = tls_new_metadata_df.iloc[idx, :].to_dict()
        cai_weights, genomes_df = blastn_run(tls_csv=blast_csv, genomes_df=genomes_df)
        tls_dict['cai'] = cai_weights
        entry_name = blast_csv.split('.')[-3]
        blast_results_dict[entry_name] = tls_dict

    return blast_results_dict, genomes_df


def save_data(blast_results_dict, genomes_df, out_fid):
    genomes_df.to_csv(out_fid + 'tls_genome_matches.csv')
    with open(out_fid + 'tls_genome_matches.json', "w") as outfile:
        json.dump(blast_results_dict, outfile)

if __name__ == "__main__":
    print('Start')
    genomes_df = pd.read_csv('../../data/processed_genomes/cai_and_16s_for_genomes.csv')
    tls_new_metadata_df = pd.read_csv('../../data/processed_tls/tls_assembly_metadata_with_blast.csv', index_col=None)
    out_fid = '../../data/tls_genome_match/'
    blast_results_dict, genomes_df = check_all_blast_res(genomes_df, tls_new_metadata_df)
    save_data(blast_results_dict, genomes_df, out_fid)
        # print(blast_results_dict)







