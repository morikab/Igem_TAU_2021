import os
import glob
import xml.etree.ElementTree as et
from Bio import SeqIO
import re
from user_input import parse_input, write_fasta
import csv

global dna, opt_path, deopt_path, end
dna = "ACGT"
opt_path = 'opt_files'
deopt_path = 'deopt_files'
end = '.fasta'

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
    kmer = "--kmer " + str(kmer) + " "
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
    for org in data_dict.keys():
        org_data = data_dict[org]
        if org_data['optimized']:
            org_path = os.path.join(opt_path, org)
        else:
            org_path = os.path.join(deopt_path, org)
        os.mkdir(org_path)
        file_path = os.path.join(org_path, org)
        ###inter and 100 not necessary for control organisms???
        write_fasta(file_path + '_100_200', list(org_data['200bp_promoters'].values()), list(org_data['200bp_promoters'].keys()))
        write_fasta(file_path + '_inter',   list(org_data['intergenic'].values()),      list(org_data['intergenic'].keys()))
        write_fasta(file_path + '_33_200',  list(org_data['third_most_HE'].values()),   list(org_data['third_most_HE'].keys()))


"""
This method runs STREME several times to create all necessary files for subsequent motif ranking and filtering
@param data_dict: a dictionary containing all input information for each organism
"""
def run_streme(data_dict):
    #create files for STREME from input dictionary
    create_files_for_meme(data_dict) ### maybe call it outside???
    out_path = 'streme_outputs'
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
This method runs MAST (not much different than mast(), created for symmetry with run_streme() )
Saves results in the same folder (doesn't create a new folder)
@param motif_path: path to the file containing motifs (already filtered for most selective and also most transcription inducing)
@param promoter_path200: path to the file containing promoter sequences of length 200 (endogenic by default, otherwise chosen by user)
@param promoter_path400: path to the file containing promoter sequences of length 400 (only if promoter_path200 is not chosen by user)
"""
def run_mast(motif_path, promoter_path200, promoter_path400=None):
    motif_name = motif_path.split(os.sep)[-1].split('.')[0] #modified motif file- not in original folder anymore (get only name of file without ext or path)
    promoter_name200 = promoter_path200.split(os.sep)[-1].split('.')[0] #get only file name without ext or path
    mast(motif_path, promoter_path200, output_path='_'.join(['motif', motif_name, 'seq', promoter_name200]))
    if promoter_path400 is not None:
        promoter_name400 = promoter_path400.split(os.sep)[-1]
        mast(motif_path, promoter_path400, output_path='_'.join(['motif', motif_name, 'seq', promoter_name400]))
    

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

@return a dictionary of promoters where the keys are the promoter names and the values are the modified sequences
"""
def modify_promoter(p_file, mast_xml):
    promoters = SeqIO.to_dict(SeqIO.parse(p_file, 'fasta'))
    new_promoters = dict()
    tree = et.parse(mast_xml)
    root = tree.getroot()

    motifs = list(root.findall('.//motif'))

    for seq in root.findall('.//sequence'): #for each promoter that was found to have significant matches to motifs
        p_name = seq.get('name')
        p_len = seq.get('length')
        promoter = str(promoters[p_name].seq).upper()
        
        for seg in seq.findall('seg'):
            for hit in seg.findall('hit'): #for each motif that matched the promoter
                if hit.get('rc') == 'n': ### check if need to always run without RC???
                    start = int(hit.get('pos'))
                    motif_idx = int(hit.get('idx'))
                    motif_seq = motifs[motif_idx].get('id').split('-')[1].upper()
                    match = hit.get('match')
                    mismatches = [m.start() for m in re.finditer(' ', match)]

                    for mm in mismatches: #for each mismatch in the local alignment of the motif and the promoter
                        pos = start + mm - 1
                        letter = find_letter(motifs[motif_idx], motif_seq, mm)
                        promoter = promoter[:pos] + letter + promoter[pos + 1:]

        new_promoters[p_name] = promoter

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
base_path = os.path.join(os.path.dirname(__file__), 'example_data')
user_inp1_raw = {
     'bacillus': {'genome_path': os.path.join(base_path, 'Bacillus subtilis.gb'),
              'optimized': True},
     'ecoli': {'genome_path': os.path.join(base_path, 'Escherichia coli.gb'),
               'optimized': False}}
user_inp1 = parse_input(user_inp1_raw)
print('finished dict')
#create_files_for_meme(user_inp1)
run_streme(user_inp1)
print('finished files')
"""
