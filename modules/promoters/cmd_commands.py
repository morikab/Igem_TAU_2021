import glob
import xml.etree.ElementTree as et
from Bio import SeqIO
from modules.shared_functions_and_vars import write_fasta
from modules.promoters.intersect_motifs_2_org_final import *
import numpy as np
#from mast_analysis import extract_motifs_from_xml


global dna, opt_path, deopt_path, end, organism_dict
dna = "ACGT"
opt_path = 'promoters_not_for_user/opt_files'
deopt_path = 'promoters_not_for_user/deopt_files'
end = '.fasta'
organism_dict = {'opt': [], 'deopt': []}
"""
organism_dict = {
    'opt': ['Escherichia_coli', 'Mycobacterium_tuberculosis', 'Pantoea_ananatis', 'Azospirillum_brasilense'],
    'deopt': ['Bacillus_subtilis', 'Sulfolobus_acidocaldarius']
    }
"""
"""
This method uses STREME to find motifs that are enriched in one set of promoters and aren't present in the other

@param primary_path: path to a .fasta file with promoter sequences, in which you want to find enriched motifs
@param control_path: path to a .fasta file with promoter sequences, in which you don't want the motifs to occur
@param kmer: length of subsequences that will be created from the input for the first step of the algorithm (default is 3)
@param output_path: name of the folder where output files will be stored
@param minw: minimal width of motifs
@param maxw: maximal width of motifs
@param pvt: p-value threshold for when to stop searching (default is 0.05)
@param ab: alphabet of the sequences (default is dna)
"""
def streme(primary_path, control_path="", kmer=3, output_path=".", minw=8, maxw=15, pvt=0.05, ab="dna"):
    if len(control_path) > 0:
        control_path = "--n " + control_path + " "
    primary_path = "--p " + primary_path + " "
    kmer = "--order " + str(kmer - 1) + " "
    output_path = "--oc " + output_path + " "
    minw = "--minw " + str(minw) + " "
    maxw = "--maxw " + str(maxw) + " "
    pvt = "--pvt " + str(pvt)
    ab = "--" + ab + " "
    command = "streme --verbosity 1 " + output_path + ab + primary_path + control_path + kmer + minw + maxw + pvt
    os.system(command)


"""
This method uses MAST to rank the presence of motifs in a sequence set

@param motif_path: path to a .txt or .html file with motif data - output of STREME in find_motifs()
@param seq_path: path to a .fasta file with promoter sequences
@param output_path: name of the folder where output files will be stored
@param ev: e-value threshold for when to stop ranking (default is 10)
"""
def mast(motif_path, seq_path, output_path=".", ev=10):
    command = "mast -oc " + output_path + " -nostatus -norc -remcorr -ev " + str(ev) + " " + motif_path + " " + seq_path
    os.system(command)


"""
This method creates fasta files from input dictionary
@param data_dict: a dictionary containing all extracted data for each input organism
"""
def create_files_for_meme(data_dict):
    os.mkdir(opt_path)
    os.mkdir(deopt_path)
    orgs_dict = data_dict['organisms']
    for org in orgs_dict.keys():
        org_data = orgs_dict[org]
        org = org.replace(' ', '_')
        if org_data['optimized']:
            organism_dict['opt'].append(org)
            org_path = os.path.join(opt_path, org)
        else:
            organism_dict['deopt'].append(org)
            org_path = os.path.join(deopt_path, org)
        os.mkdir(org_path)
        file_path = os.path.join(org_path, org)
        #inter and 100 not necessary for control organisms???
        write_fasta(file_path + '_100_200', list(org_data['200bp_promoters'].values()), list(org_data['200bp_promoters'].keys()))
        write_fasta(file_path + '_inter',   list(org_data['intergenic'].values()),      list(org_data['intergenic'].keys()))
        write_fasta(file_path + '_33_200',  list(org_data['third_most_HE'].values()),   list(org_data['third_most_HE'].keys()))

    promoter_dict = data_dict['selected_prom']
    new_p_dict = dict()
    for p in promoter_dict.keys():
        new_p_key = p.replace(' ', '_')
        new_p_dict[new_p_key] = promoter_dict[p]
    write_fasta('promoters', list(new_p_dict.values()), list(new_p_dict.keys()))


"""
This method executes one run of STREME for a specific set of files
@param name1: name of the primary file
@param start1: path to the directory containing the file name1
@param name2: name of the control file
@param start2: path to the directory containing the file name2
@param out_path: (optional) directory in which to store output files from STREME
"""
def one_streme(name1, start1, name2, start2, out_path=""):
    k = 3
    minw = 6
    maxw = 20

    path1 = os.path.join(start1, name1 + end)
    path2 = os.path.join(start2, name2 + end)
    out_dir = '_'.join([name1, name2])
    if len(out_path) > 0:
        out_dir = os.path.join(out_path, out_dir)
    streme(primary_path=path1, control_path=path2, kmer=k, output_path=out_dir, minw=minw, maxw=maxw)


"""
This method runs STREME several times to create all necessary files for subsequent motif ranking and filtering
"""
def run_streme():
    out_path = 'promoters_not_for_user/streme_outputs'
    os.mkdir(out_path)
    for opt_org in glob.glob(os.path.join(opt_path, "*", "")):

        #first run: opt organism vs. intergenic
        org_name = opt_org.split(os.sep)[-2]
        dir_name = os.path.join(opt_path, org_name)
        name1 = '_'.join([org_name, '100', '200'])
        name2 = '_'.join([org_name, 'inter'])
        one_streme(name1, dir_name, name2, dir_name, out_path)

        #second run: opt organism HE vs. all deopt organisms HE
        for deopt_org in glob.glob(os.path.join(deopt_path, "*", "")):
             de_org_name = deopt_org.split(os.sep)[-2]
             de_dir_name = os.path.join(deopt_path, de_org_name)
             name1 = '_'.join([org_name, '33', '200'])
             name2 = '_'.join([de_org_name, '33', '200'])
             one_streme(name1, dir_name, name2, de_dir_name, out_path)


"""
This method runs MAST (not much different than mast(), created for symmetry with run_streme() )
Saves results in the same folder (doesn't create a new folder)
@param motif_path: path to the file containing motifs (already filtered for most selective and also most transcription inducing)
@param promoter_path200: path to the file containing promoter sequences of length 200 (endogenic by default, otherwise chosen by user)
@param promoter_path400: path to the file containing promoter sequences of length 400 (only if promoter_path200 is not chosen by user)
"""
def run_mast(motif_path, promoter_path):
    motif_name = motif_path.split(os.sep)[-1].split('.')[0] #modified motif file- not in original folder anymore (get only name of file without ext or path)
    promoter_name = promoter_path.split(os.sep)[-1].split('.')[0] #get only file name without ext or path
    mast(motif_path, promoter_path, output_path='_'.join(['motif', motif_name, 'seq', promoter_name]))


"""
This method filters motifs outputted by STREME according to a given dictionary of motifs.
The motifs which stay in the file are those that were found to be both selective and transcription inducing
@param motif_d: a dictionary of motifs where the key is the motif id (as defined by STREME) and the value is its correlation coefficient
@param streme_xml: path to a STREME output file in XML format from the STREME run that differentiates between primary and control organisms

@return name of modified xml file, located in the working directory
"""
def filter_motifs(motif_d, streme_xml):
    chosen = list(motif_d.keys())
    name = streme_xml.split(os.sep)[-2] #get parameters-written in the name of the folder that was outputted from streme (file name is simply 'streme'..)
    tree = et.parse(streme_xml)
    root = tree.getroot()
    for element in root.findall('.//motifs'):
        for motif in element.findall('motif'):
            if motif.get('id') not in chosen:
                element.remove(motif)

    tree.write(name + '_modified.xml')
    return name + '_modified.xml'


"""
This method modifies promoters to align better with motifs that match them.
NOTE: currently only modifies sites where motif match was with the original sequence and not the reverse complement
@param p_file: path to a fasta file of promoters (must be same promoters as the MAST run)
@param mast_xml: MAST output XML file

@return a dictionary of promoters where the keys are the promoter names and the values are dictionaries with promoter information
"""
def modify_promoter(p_file, mast_xml):
    promoters = SeqIO.to_dict(SeqIO.parse(p_file, 'fasta'))
    #print(list(promoters.keys())[0:10])
    new_promoters = dict()
    tree = et.parse(mast_xml)
    root = tree.getroot()

    motifs = list(root.findall('.//motif'))

    for seq in root.findall('.//sequence'): #for each promoter that was found to have significant matches to motifs
        strand = seq.find('score').get('strand')
        if strand == 'forward':
            p_name = seq.get('name').strip().replace(' ', '_')
            p_comment = seq.get('comment').strip().replace(' ', '_')
            p_name = '_'.join([p_name, p_comment])
            #print(p_name)
            p_len = seq.get('length')
            promoter = str(promoters[p_name].seq).upper()
            org_promoter = promoter
            evalue = float(seq.find('score').get('evalue'))
            for seg in seq.findall('seg'):
                for hit in seg.findall('hit'): #for each motif that matched the promoter
                    start = int(hit.get('pos'))
                    motif_idx = int(hit.get('idx'))
                    motif_seq = motifs[motif_idx].get('id').split('-')[1].upper()
                    match = hit.get('match')
                    mismatches = [m.start() for m in re.finditer(' ', match)]

                    for mm in mismatches: #for each mismatch in the local alignment of the motif and the promoter
                        pos = start + mm - 1
                        letter = find_letter(motifs[motif_idx], motif_seq, mm)
                        promoter = promoter[:pos] + letter + promoter[pos + 1:]

        new_promoters[p_name] = {'original_seq': org_promoter, 'evalue': evalue, 'modified_seq': promoter}

    return new_promoters


"""
This method finds the letter in the motif that will replace the current letter in the promoter, to match the motif better
@param motif: an Element Tree element encoding the relevant motif
@param motif_seq: the sequence of the motif
@param index: the index of the letter in the motif

@return the letter with highest PSSM score
"""
def find_letter(motif, motif_seq, index):
    if motif_seq[index] in dna:
        return motif_seq[index]
    pos = list(motif.findall('pos'))[index]
    pssm = [float(pos.get(letter)) for letter in dna]
    best_letter = dna[pssm.index(max(pssm))]
    return best_letter


"""
This method writes the results of the promoter optimization into a text file
@param promoter_dict: a dictionary where promoter names are keys and a dictionary containing promoter info is the value; output of modify_promoter
@param p_num: number of promoters to write into the file. These will be the ones with best E-value scores.
"""
def write_results(promoter_dict, p_num):
    text = ''
    for p_name, p_info in list(promoter_dict.items())[:p_num]:
        text += 'promoter name:\n' + p_name.replace('_', ' ') + '\n\n'
        org_seq = p_info['original_seq']
        text += 'original sequence:\n' + org_seq + '\n\n'
        new_seq = p_info['modified_seq']
        text += 'modified sequence:\n' + new_seq + '\n\n'
        evalue = p_info['evalue']
        text += 'E-value:\n' + str(evalue) + '\n\n'

    with open('chosen_promoter.txt', 'w') as f:
        f.write(text)

######################################################
###### Motif filtering for multi-organism model ######
######################################################

"""
This method checks if an intergenic STREME output file is from
an optimized organism
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
This method gets all xml STREME output files from intergenic runs

@return: a list of file names
"""
def find_all_inter_files():
    all_files = glob.glob(os.path.join('streme_outputs', "**", "*.xml"), recursive=True)
    inter_files = [f for f in all_files if 'inter' in f]
    return inter_files


"""
This method gets all xml STREME output files from selective runs

@return: a list of file names
"""
def find_all_selective_files():
    all_files = glob.glob(os.path.join('streme_outputs', "**", "*.xml"), recursive=True)
    selective_files = [f for f in all_files if 'inter' not in f]
    return selective_files


"""
This method creates a new xml STREME file with all intergenic motifs in all runs

@return: a pointer to an ElementTree object containing the data in xml format
"""
def unionize_motifs():
    ### all outputs or only opt organisms??? currently all #####################
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

    base_tree.write('unionized_motifs.xml')
    return base_tree


"""
This method finds all intergenic motifs with low correlation to selective ones
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
This method creates a final STREME file with only motifs that are both
selective and intergenic in all optimized organisms
@param base_tree: an ElementTree object containing all intergenic motifs
@param D1: threshold for intergenic correlation calculation
@param D2: threshold for selective correlation calculation
"""
def multi_filter_motifs(base_tree, D1, D2):
    selective_files = find_all_selective_files()
    inter_files = find_all_inter_files()
    C_set = extract_pssm_from_xml('unionized_motifs.xml')

    base_root = base_tree.getroot()
    base_element = base_root.find('.//motifs')

    #first check: each set S_x where x in A, has at least one motif with correlation > D1 to the C_set motif
    to_delete1 = get_motifs_to_delete(C_set, inter_files, D1)

    #second check: each set S_xy where x in A and y in B, has at least one motif with correlation > D2 to the C_set motif
    to_delete2 = get_motifs_to_delete(C_set, selective_files, D2)

    to_delete = set.union(to_delete1, to_delete2)
    print(to_delete)
    for i, motif in enumerate(base_element.findall('motif')):
        if i in to_delete:
            base_element.remove(motif)

    base_tree.write('final_motif_set.xml')


"""
This method calls the xml pipeline to create final motif set
@param D1: threshold for intergenic correlation calculation
@param D2: threshold for selective correlation calculation
"""
def create_final_motif_xml(D1: float, D2: float):
    tree = unionize_motifs()
    multi_filter_motifs(tree, D1, D2)




base_path = os.path.join(os.path.dirname(__file__), 'example_data')
user_inp_raw = {
    'sequence': os.path.join(base_path, 'mCherry_original.fasta'),
    'selected_promoters': None,
    'organisms': {
                   'opt1': {'genome_path': os.path.join(base_path, 'Escherichia coli.gb'),
                            'optimized': True,
                            'expression_csv': None},

                    'deopt1': {'genome_path': os.path.join(base_path, 'Bacillus subtilis.gb'),
                               'optimized': False,
                               'expression_csv': None},

                    'deopt2': {'genome_path': os.path.join(base_path, 'Sulfolobus acidocaldarius.gb'),
                              'optimized': False,
                              'expression_csv': None},

                    'opt2': {'genome_path': os.path.join(base_path, 'Mycobacterium tuberculosis.gb'),
                             'optimized': True,
                             'expression_csv': None},

                    'opt3': {'genome_path': os.path.join(base_path, 'Pantoea ananatis.gb'),
                             'optimized': True,
                             'expression_csv': None},

                    'opt4': {'genome_path': os.path.join(base_path, 'Azospirillum brasilense.gb'),
                             'optimized': True,
                             'expression_csv': None}
    }
}

#input_dict = user_IO.UserInputModule.run_module(user_inp_raw) #keys: sequence, selected_prom, organisms


