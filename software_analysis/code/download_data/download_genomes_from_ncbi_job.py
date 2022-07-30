# -*- coding: utf-8 -*-
import os
import subprocess

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


def download_files():
    # with open("rna_files.txt", "r") as rna_files:
    #     file_names = rna_files.readlines()
    #     for file_full_name in file_names:
    #         f_name = file_full_name.strip().split('/')[-1]
    #         rna_destination_file = destination_dir+'ncbi_genome_rna'
    #         downloaded_files = os.listdir(rna_destination_file)
    #         if f_name not in downloaded_files:
    #             if f_name[:-3] not in downloaded_files:
    #                 command = 'wget -P ' + rna_destination_file +  ' ' + file_full_name.strip() +' --no-check-certificate'
    #                 run_cmd(command)
    #                 print(command)
    #             command = 'gzip -d ' + rna_destination_file + '/'+  f_name
    #             run_cmd(command)
    #             print(command)
    #             continue

    with open("cds_files.txt", "r") as cds_file:
        file_names = cds_file.readlines()
        for file_full_name in file_names:
            f_name = file_full_name.strip().split('/')[-1]
            cds_destination_file = destination_dir+'ncbi_genome_cds'
            downloaded_files = os.listdir(cds_destination_file)
            if f_name not in downloaded_files:
                if f_name[:-3] not in downloaded_files:
                    command = 'wget -P ' + cds_destination_file +  ' ' + file_full_name.strip() +' --no-check-certificate'
                    run_cmd(command)
                    print(command)
                command = 'gzip -d ' + cds_destination_file + '/' + f_name
                run_cmd(command)
                print(command)
                continue


if __name__ == "__main__":
    print("Start")
    download_files()
    print("The end")





