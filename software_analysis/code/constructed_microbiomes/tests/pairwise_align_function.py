from Bio import pairwise2, SeqIO
import pandas as pd
import json


def align_seq_to_fasta(seq_fid, genomes_fid, output_dir):
    # print(genomes_fid)

    genome_slice_name = genomes_fid.split('/')[-1][:-6]
    csv_fname =  output_dir + genome_slice_name+ '_local_align.csv'
    json_fname = output_dir + genome_slice_name + '_local_align.json'

    seq = [str(i.seq) for name, i in
                 SeqIO.index(seq_fid,
                             "fasta").items()][0]

    genomes = {name: str(i.seq) for name, i in
                 SeqIO.index(genomes_fid,
                             "fasta").items()}
    scores = {name:pairwise2.align.localxx(seq, rrna_seq)[0].score
              for name, rrna_seq in genomes.items()}
    scores_df = pd.DataFrame(scores, index=['alignment_score']).T

    print(scores_df)


    scores_df.to_csv(csv_fname)
    with open(json_fname, 'w') as fp:
        json.dump(scores, fp)

