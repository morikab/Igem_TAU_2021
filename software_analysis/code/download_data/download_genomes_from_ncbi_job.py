# -*- coding: utf-8 -*-
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
    # with open("cds_files.txt", "r") as cds_file:
    #     file_names = cds_file.readlines()
    #     for file_name in file_names:
    #         cds_destination_file = destination_dir+'ncbi_genome_cds'
    #         command = 'wget -P ' + cds_destination_file +  ' ' + file_name.strip() +' --no-check-certificate'
    #         run_cmd(command)
    #         print(command)
    #     command = 'gzip -d ' + cds_destination_file + '/*'
    #     print(command)
    #     run_cmd(command)

    with open("rna_files.txt", "r") as rna_files:
        file_names = rna_files.readlines()
        for file_name in file_names:
            rna_destination_file = destination_dir+'ncbi_genome_rna'
            command = 'wget -P ' + rna_destination_file +  ' ' + file_name.strip() +' --no-check-certificate'
            run_cmd(command)
            print(command)
        command = 'gzip -d ' +  rna_destination_file +'/*'
        print(command)
        run_cmd(command)

if __name__ == "__main__":
    print("Start")
    download_files()
    print("The end")





