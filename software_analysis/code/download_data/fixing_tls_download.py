# -*- coding: utf-8 -*-
import codecs
import os
import subprocess

destination_dir = '../../data/genbank_tls'
index_file = codecs.open("tls_index.html", 'r')


def txt_file_to_lst (fid):
    txt_file = open(fid, "r")
    file_content = txt_file.read()
    print("The file content are: ", file_content)

    content_list = file_content.split(",")
    txt_file.close()
    print("The list is: ", content_list)
    return content_list


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
    Lines = index_file.readlines()


    count = 0
    run_cmd('cd '+ destination_dir )
    for line in Lines:
        count += 1
        line = line.strip()
        if "fsa_nt.gz" in line or ".mstr.gbff.gz" in line:
            f_name = line.split('"')[1]
            print(os.listdir(destination_dir))
            command = 'wget -P ' + destination_dir +  ' https://ftp.ncbi.nlm.nih.gov/genbank/tls/K/' +f_name +' --no-check-certificate'

            run_cmd(command)
            print(command)
            # destination_file = destination_dir + '/' + f_name
            command = 'gzip -d ' + destination_dir +'/' +f_name
            print(command)
            run_cmd(command)




# if __name__ == "__main__":
#     print("Start!")
#     download_files()
#     print("The end")

print(len(txt_file_to_lst('all_files.txt')))

