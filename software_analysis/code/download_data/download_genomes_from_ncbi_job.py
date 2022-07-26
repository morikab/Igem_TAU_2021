# -*- coding: utf-8 -*-
import subprocess

destination_dir = '../../data/'

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
    with open("cds_files.txt", "r") as cds_file:
        file_names = cds_file.readlines()
        for file_name in file_names:
            cds_destination_file = destination_dir+'ncbi_genome_cds'
            run_cmd("wget" +  file_name.strip() + " ‐‐directory-prefix " + cds_destination_file)
            print("wget" +  file_name.strip() + " ‐‐directory-prefix " + cds_destination_file)
            run_cmd("gzip -d" +  cds_destination_file)

    with open("rna_files.txt", "r") as rna_files:
        file_names = rna_files.readlines()
        for file_name in file_names:
            rna_destination_file = destination_dir+'ncbi_genome_rna'
            run_cmd("wget " +  file_name.strip() + " ‐‐directory-prefix " + rna_destination_file)
            print("wget " +  file_name.strip() + " ‐‐directory-prefix " + rna_destination_file)
            run_cmd("gzip -d" +  rna_destination_file)



if __name__ == "__main__":
    print("Start")
    download_files()
    print("The end")

