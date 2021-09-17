import os
import re
import xml.etree.ElementTree as et
from Bio import SeqIO
from globals_and_shared_methods import *

"""
Modifies best promoter to align better with motifs that match them.
Only modifies sites where motif match was with the original sequence and not the reverse complement

@param p_file: path to a fasta file of promoters (must be same promoters as the MAST run)
@param mast_folder: path to MAST output directory

@return: name, original sequence, modified sequence and E-value of promoter with best E-value (first in xml file)
"""
def modify_promoter(p_file, mast_folder):
    promoters = SeqIO.to_dict(SeqIO.parse(p_file, 'fasta'))

    mast_file = os.path.join(mast_folder, 'mast.xml')

    tree = et.parse(mast_file)
    root = tree.getroot()
    
    motifs = list(root.findall('.//motif'))

    seq = root.find('.//sequence')
    strand = seq.find('score').get('strand')
    if strand == 'forward':
        p_name = seq.get('name')
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

    #replace '_' with ' ' (only after '|' sign)
    if '|' in p_name:
        temp = p_name.split('|')
        for i in range(1, len(temp)):
            temp[i] = temp[i].replace('_', ' ')
        p_name = "|".join(temp)
        
    return p_name, org_promoter, promoter, evalue
                

"""
Finds the letter in the motif that will replace the current letter in the promoter, to match the motif better
This is calculated based on PSSM

@param motif: an Element Tree element encoding the relevant motif
@param motif_seq: the sequence of the motif
@param index: the index of the letter in the motif

@return: the letter with best PSSM score
"""
def find_letter(motif, motif_seq, index):
    if motif_seq[index] in dna:
        return motif_seq[index]
    pos = list(motif.findall('pos'))[index]
    pssm = [float(pos.get(letter)) for letter in dna]
    best_letter = dna[pssm.index(max(pssm))]
    return best_letter
