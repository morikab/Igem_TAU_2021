from Bio import SeqIO

#todo: for all shared functions and variables (or things that should be defined as external data that we might want to change
# later, like the sij) add them to this file and import into correct module ( look at RE!!!)

def fasta_to_dict(fasta_fid):
    fasta_dict = {record.description:str(record.seq) for record in SeqIO.parse(fasta_fid, 'fasta') }
    return fasta_dict

def write_fasta(fid, list_seq, list_name):
    ofile = open(fid + '.fasta', "w+")
    for i in range(len(list_seq)):
        ofile.write(">" + list_name[i] + "\n" + list_seq[i] + "\n")
    ofile.close()

# write ideas for the promoter model
def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    reverse_complement = "".join(complement.get(base, base) for base in reversed(seq))
    return reverse_complement