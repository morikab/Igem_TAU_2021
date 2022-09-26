from Bio import pairwise2, SeqIO
import pandas as pd


def align_seq_to_fasta(seq_fid, genomes_fid):
    # print(genomes_fid)
    seq = [str(i.seq) for name, i in
                 SeqIO.index(seq_fid,
                             "fasta").items()][0]

    genomes = {name: str(i.seq) for name, i in
                 SeqIO.index(genomes_fid,
                             "fasta").items()}
    scores = {name:pairwise2.align.localxx(seq, rrna_seq)[0].score
              for name, rrna_seq in genomes.items()}
    scores_df = pd.DataFrame.from_dict(scores)
    print(scores_df)
    scores_df.to_csv(genomes_fid[:6] + '_dynamic_programing.csv')
