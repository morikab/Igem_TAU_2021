from RE.multi_org_functions import *
from RE.re_functions import insert_site_CDS, translate
import json
from Bio.SeqIO import read
from modules.RE import REModule

with open('C:\\Users\\labuser\\Documents\\Liyam\\tamir projects\\igem article\\Igem_TAU_2021\\model_analysis\\data_for_analysis\\org_name_to_dict.json') as f:
    org_dict = json.load(f)

gene = read('C:\\Users\\labuser\\Documents\\Liyam\\tamir projects\\igem article\\Igem_TAU_2021\\model_analysis\\zorA anti-phage defense.fasta', 'fasta')
cds_nt = str(gene.seq)

optimized_org_names = list(org_dict.keys())[:5]

deoptimized_org_names = list(org_dict.keys())[6:13]

cds_aa = translate(cds_nt)

optimized_RE_dict = multi_organism_RE_dict(optimized_org_names, cds_aa)
deoptimized_RE_dict = multi_organism_RE_dict(deoptimized_org_names, cds_aa)

re_positions,  add_cds_nt= insert_site_CDS(deoptimized_RE_dict, cds_nt)
final_cds_nt = multi_org_remove_site(optimized_RE_dict, add_cds_nt)


opt_sites = multi_org_final_found_sites(optimized_RE_dict, cds_nt)
print('opt sites', opt_sites)
deopt_sites = multi_org_final_found_sites(deoptimized_RE_dict, final_cds_nt)
print('deopt sites', deopt_sites)
