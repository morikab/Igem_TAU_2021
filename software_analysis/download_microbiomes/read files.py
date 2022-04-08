import codecs
index_file=codecs.open("index.html", 'r')
Lines = index_file.readlines()



commands = open("commands.sh", 'w') ### writes the files to the directory is was activated from

count = 0
for line in Lines:
    count += 1
    line = line.strip()
    print(line)
    if "fsa_nt.gz" in line or ".mstr.gbff.gz" in line:
        f_name = line.split('"')[1]
        commands.write(f'wget https://ftp.ncbi.nlm.nih.gov/genbank/tls/K/{f_name} --no-check-certificate &\n')

commands.close()