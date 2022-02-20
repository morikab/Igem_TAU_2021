nt_to_aa = {
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
    'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
    'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
    'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_',
    'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W',
}


synonymous_codons = {
    "C": ["TGT", "TGC"],
    "D": ["GAT", "GAC"],
    "S": ["TCT", "TCG", "TCA", "TCC", "AGC", "AGT"],
    "Q": ["CAA", "CAG"],
    "M": ["ATG"],
    "N": ["AAC", "AAT"],
    "P": ["CCT", "CCG", "CCA", "CCC"],
    "K": ["AAG", "AAA"],
    "_": ["TAG", "TGA", "TAA"],
    "T": ["ACC", "ACA", "ACG", "ACT"],
    "F": ["TTT", "TTC"],
    "A": ["GCA", "GCC", "GCG", "GCT"],
    "G": ["GGT", "GGG", "GGA", "GGC"],
    "I": ["ATC", "ATA", "ATT"],
    "L": ["TTA", "TTG", "CTC", "CTT", "CTG", "CTA"],
    "H": ["CAT", "CAC"],
    "R": ["CGA", "CGC", "CGG", "CGT", "AGG", "AGA"],
    "W": ["TGG"],
    "V": ["GTA", "GTC", "GTG", "GTT"],
    "E": ["GAG", "GAA"],
    "Y": ["TAT", "TAC"],
}

DEFAULT_ORGANISM_PRIORITY = 5


def write_fasta(fid, list_seq, list_name):
    ofile = open(fid + '.fasta', "w+")
    for i in range(len(list_seq)):
        ofile.write(">" + str(list_name[i]) + "\n" + str(list_seq[i]) + "\n")
    ofile.close()


# write ideas for the promoter model
def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    reverse_complement = "".join(complement.get(base, base) for base in reversed(seq))
    return reverse_complement


def translate(seq, table=nt_to_aa):  # DONE
    protein = ""
    if len(seq) % 3 == 0:
        for i in range(0, len(seq), 3):
            codon = seq[i:i + 3]
            protein += table[codon]
        return protein
    else:
        return ValueError('len(seq)%3 !=0')


def unique(list1):
    return sorted(set(list1))
