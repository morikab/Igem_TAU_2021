from RE.multi_org_functions import *
from user_input import user_inp

print('\n\n\n###############################')
print('# RESTRICTION ENZYME ANALYSIS #')
print('###############################')

cds_nt = user_inp['sequence']
cds_aa = translate(cds_nt)

#finding the optimixed and deoptimized RE dictionaries
optimized_org_names, deoptimized_org_names = parse_inp1(user_inp)
print('\nOptimized organisms:')
optimized_RE_dict = multi_organism_RE_dict(optimized_org_names, cds_aa)
print('\nDeptimized organisms:')
deoptimized_RE_dict = multi_organism_RE_dict(deoptimized_org_names, cds_aa)

# inserting sites for all optimized organisms
# todo: use multi_organism_RE_dict after completing and finish this- do not call it just cds_nt
add_cds_nt = multi_org_insert_site(deoptimized_RE_dict, cds_nt)

# removing sites for all deoptimized organism
final_cds_nt = multi_org_remove_site(optimized_RE_dict, add_cds_nt)

# checking result:

print(f'Initial sequence before translation and restriction enzyme optimization: '
      f'{cds_nt}\n')
total_sequence_analysis(optimized_RE_dict, deoptimized_RE_dict, final_cds_nt)

print(f'Final sequence after translation and restriction enzyme optimization: '
      f'{final_cds_nt}\n')
total_sequence_analysis(optimized_RE_dict, deoptimized_RE_dict, final_cds_nt)

