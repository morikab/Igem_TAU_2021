from Bio import SeqIO

# module load Python3.9Plus
# python -c 'from cross_tls_with_genome_blast_job.py import filename_to_sent_job; filename_to_sent_job.filename_to_sent_job('vvv')'

output_fid = '../../../data/tested_results/genomes_16s_sliced/'
filtered_genome_fid = '../../../data/processed_genomes/filtered/cai_and_16s_for_genomes_filtered.fasta'
genomes = {name: str(i.seq) for name, i in
              SeqIO.index(filtered_genome_fid,
                          "fasta").items()}

n_in_split = 250
split_idxs_raw = [i for i in range(len(genomes)) if i%n_in_split ==0]
split_idxs = [(split_idxs_raw[i], split_idxs_raw[i+1]) for i in range(len(split_idxs_raw)-1)]
print(split_idxs_raw[-1])
split_idxs.append((split_idxs_raw[-1], len(genomes)))

def write_fasta(fid, list_seq, list_name):
    ofile = open(fid + '.fasta', "w+")
    for i in range(len(list_seq)):
        ofile.write(">" + str(list_name[i]) + "\n" + str(list_seq[i]) + "\n")
    ofile.close()


fasta_files = {}
for idx, split in enumerate(split_idxs):
    start, end = split
    fasta_name = f'{idx}_start_{start}_end_{end}'
    fasta_dict = dict(list(genomes.items())[start:end])
    fasta_files[fasta_name] = fasta_dict
    write_fasta(fid = output_fid + fasta_name,
                list_seq = list(fasta_dict.values()),
                list_name = list(fasta_dict.keys()))

    print(fasta_name)

