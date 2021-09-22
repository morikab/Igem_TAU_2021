import glob
from modules.promoters.globals_and_shared_methods import *
from modules.promoters.intersect_motifs_2_org_final import *

######################################################
###### Motif filtering for multi-organism model ######
######################################################

"""
Checks if an intergenic STREME output file is from an optimized organism

@param fname: name of a STREME output file

@return: True if the file is from an optimized organism, False otherwise
"""
def is_optimized(fname):
    name = fname.split(os.sep)[-2]
    org_name = '_'.join(name.split('_')[:2])
    if org_name in organism_dict['opt']:
        return True
    return False
    

"""
Finds all xml STREME output files from intergenic runs

@return: a list of file names
"""
def find_all_inter_files():
    all_files = glob.glob(os.path.join(start, 'streme_outputs', "**", "*.xml"), recursive=True)
    inter_files = [f for f in all_files if 'inter' in f]
    return inter_files


"""
Finds all xml STREME output files from selective runs

@return: a list of file names
"""
def find_all_selective_files():
    all_files = glob.glob(os.path.join(start, 'streme_outputs', "**", "*.xml"), recursive=True)
    selective_files = [f for f in all_files if 'inter' not in f]
    return selective_files


"""
Creates a new xml STREME file with all intergenic motifs in all runs

@return: a pointer to an ElementTree object containing the data in xml format
"""
def unionize_motifs():
    inter_files = [f for f in find_all_inter_files() if is_optimized(f)]

    base_file = inter_files[0]
    base_tree = et.parse(base_file)
    base_root = base_tree.getroot()
    base_element = base_root.find('.//motifs')
    
    for xml_file in inter_files[1:]:
        tree = et.parse(xml_file)
        root = tree.getroot()
        element = root.find('.//motifs')
        for motif in element.findall('motif'):
            base_element.append(motif)

    base_tree.write(os.path.join(start, 'unionized_motifs.xml'))
    return base_tree


"""
Finds all intergenic motifs with low correlation to selective ones

@param C_set: dictionary of the unionized intergenic motifs with pssms
@param file_list: a list of files to calculate their motifs' correlation to C_set
@threshold: used to decide which correlation scores are too low

@return: a list of motif indices to delete from the xml file
"""
def get_motifs_to_delete(C_set, file_list, threshold: float):
    to_delete = set()
    for file in file_list:
        if is_optimized(file):
            S_x = extract_pssm_from_xml(file)
            corr_df, pval_df = compare_pssm_sets(C_set, S_x)
            corr_df['max'] = corr_df.max(axis=1)
            low_idx = list(np.where(corr_df['max'] < threshold)[0])
            to_delete.update(low_idx) #mark all motifs with low correlation

    return to_delete


"""
Creates a final STREME file with only motifs that are both selective and intergenic in all optimized organisms

@param base_tree: an ElementTree object containing all intergenic motifs
@param D1: threshold for intergenic correlation calculation
@param D2: threshold for selective correlation calculation

@return: the name of the new motif file created
"""
def multi_filter_motifs(base_tree, D1, D2):
    selective_files = find_all_selective_files()
    inter_files = find_all_inter_files()
    C_set = extract_pssm_from_xml(os.path.join(start, 'unionized_motifs.xml'))

    base_root = base_tree.getroot()
    base_element = base_root.find('.//motifs')

    #first check: each set S_x where x in A, has at least one motif with correlation > D1 to the C_set motif
    to_delete1 = get_motifs_to_delete(C_set, inter_files, D1)

    #second check: each set S_xy where x in A and y in B, has at least one motif with correlation > D2 to the C_set motif
    to_delete2 = get_motifs_to_delete(C_set, selective_files, D2)
    
    to_delete = set.union(to_delete1, to_delete2)
    
    for i, motif in enumerate(base_element.findall('motif')):
        if i in to_delete:
            base_element.remove(motif)

    fname = 'final_motif_set.xml'
    fname = os.path.join(start, fname)
    base_tree.write(fname)
    
    return fname


"""
Calls the xml pipeline to create final motif set for MAST

@param D1: threshold for intergenic correlation calculation
@param D2: threshold for selective correlation calculation

@return: the name of the new motif file created
""" 
def create_final_motif_xml(D1: float, D2: float):
    tree = unionize_motifs()
    motif_file = multi_filter_motifs(tree, D1, D2)
    
    return motif_file

