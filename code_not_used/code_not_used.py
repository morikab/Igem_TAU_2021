import os
import xml.etree.ElementTree as et


"""
Filters motifs outputted by STREME according to a given dictionary of motifs.
The motifs that are kept in the file are those that were found to be both selective and transcription inducing

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
Writes the results of the promoter optimization into a text file

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
        
