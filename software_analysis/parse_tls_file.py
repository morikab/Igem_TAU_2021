from Bio import SeqIO

def parse_tls(fid):
    tls_seqs = {}

    with open(fid, "rt") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            tls_seqs[record.description] = record.seq
            print(record.description)



parse_tls('download_microbiomes/id.vdb_wgsnc.0301.2019.KAAB.1.fsa_nt')