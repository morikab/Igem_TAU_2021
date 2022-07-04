from RE.multi_org_functions import *
from RE.re_functions import insert_site_CDS, translate
import json
from Bio.SeqIO import read
from modules.testing_for_modules import *
from shared_functions_and_vars import *
import pandas as pd
import statistics
import matplotlib.pyplot as plt
import math

base_path = os.path.join(os.path.dirname(__file__))

f = open(os.path.join(base_path, 'REbase.txt'))
content = f.readlines()
REbase = [x.strip() for x in content]


def REbase_org(org):  # DONE
    n_unknown_in_site = 0
    '''
    parsing the REbase and creating a dictionary with the relevant sites from
    enzymes originated in org
    :param org: scientific name of investigated organism
    :param cds_aa: aa sequence of the cds
    :return: a dict in the following format:
                enzyme_dict[enzyme_name] = {
                'ambiguous_site':ambiguous_site, #contaiinng characters from the ambiguous_code dictionary
                'nt_to_aa': nt_to_aa} #non ambigiosus sites translated (in all 3 RFs)
    '''
    org = org.split(' ')
    org = org[0]# todo: this uses only the fam name and not specie
    print(org)
    RE_org_ind = [i for i, k in enumerate(REbase)
                  if org in k]  ## index of the 3rd section of the "block"
    print(RE_org_ind)
    enzyme_dict = {}  ## REbase dict.
    for k in RE_org_ind:
        org_list = []
        site = REbase[k + 2][3:]
        if '?' in site or ',' in site:
            n_unknown_in_site +=1
            continue
        enzyme_name = REbase[k - 2][3:]
        org_list.append(enzyme_name)

    print(f'number of sites with unknown char {n_unknown_in_site}')
    return enzyme_dict


n_unknown = 0
RE_org_name = [k[3:] for i, k in enumerate(REbase) if '<3>' in k]
RE_org_site = [k[3:] for i, k in enumerate(REbase) if '<5>' in k]

org_dict = {}
for index, org  in enumerate(RE_org_name):
    if 'Unidentified bacterium' in org:
        n_unknown +=1
        continue
    org = ' '.join(str.split(org, ' ')[:2])
    try:
        org_dict[org].append(RE_org_site[index])
    except:
        org_dict[org] = [RE_org_site[index]]

RE_org_list = unique(RE_org_name)
print('total number of organisms: ', len(RE_org_list))
print('total number of organisms unknown organisms: ', len([i for i in RE_org_list if 'Unidentified bacterium' in i]))

print(org_dict)

org_n_dict  = {k:len(v) for k,v in org_dict.items()}
log_org_n_dict  = {k:math.log10(len(v)) for k,v in org_dict.items()}

print(f'total number of sites {sum(org_n_dict.values())}')
print(f'avg number of sites {sum(org_n_dict.values())/len(org_n_dict)}')
print(f'std number of sites {statistics.stdev(org_n_dict.values())}')

more_than100 = {k:i for  k,i in org_n_dict.items() if i>50}
print(more_than100)
#
# plt.hist(org_n_dict.values(), bins= 100)
# plt.show()
#
# plt.hist(log_org_n_dict.values(), bins= 10)
# plt.show()