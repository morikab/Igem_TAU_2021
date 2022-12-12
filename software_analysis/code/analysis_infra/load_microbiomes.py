from os import listdir
from os.path import isfile, join
import json
import matplotlib.pyplot as plt

class obj:
    def __init__(self, dict1):
        self.__dict__.update(dict1)




if __name__ == "__main__":
    tls_match_dir = '../../data/tls_genome_match/'
    tls_files = [join(tls_match_dir, f) for f in listdir(tls_match_dir) if isfile(join(tls_match_dir, f))]

    for tls_file in tls_files:
        print(tls_file)
        tls = json.load(open(tls_file), object_hook=obj)


