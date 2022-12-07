from cross_tls_with_genome_blast_job import *



if __name__ == "__main__":
    print('Start')
    command_list = run_all_tls('../../data/genbank_tls/')
    print(len(command_list))
    job_files = []
    for idx, command in enumerate(command_list):
        filename = str(idx) + '_blast_job.sh'
        job_files.append(filename)
        write_job([command], 'tls_to_16s_blast_5_hits/' + filename)


    master_commands = [filename_to_sent_job(sh_file, cput='15:00:00') for sh_file in job_files]
    f = open('tls_to_16s_blast_5_hits/rerun_killed_jobs_mstr.sh', 'w')
    f.write(
        '#!/bin/sh \n cd /tamir1/liyamlevi/projects/communique/Igem_TAU_2021/software_analysis/code/raw_data_analysis/tls_to_16s_blast_5_hits\n')
    for line in master_commands:
        f.write(line + '\n')
    f.close()
        # proc = subprocess.Popen(command, shell=True,
        #                         stdin=subprocess.PIPE, stdout=subprocess.PIPE)
        # out, err = proc.communicate()
        # print("{} : {}".format(command, out.decode()))
    print("The end")