from Bio import pairwise2, SeqIO
import matplotlib.pyplot as plt
import time
import heapq


metagenome = {name: str(i.seq) for name, i in
              SeqIO.index("../../../data/tested_results/KDVY_example_metagenome/tls.KDVY.1.fsa_nt",
                          "fasta").items()}

scores = []
for name1, seq1 in metagenome.items():
    amplicon_scores = []
    for name2, seq2 in metagenome.items():
        if name2 == name1:
            continue
        alignments = pairwise2.align.localxx(seq2, seq1)[0].score
        v = pairwise2.align.localxx(seq2, seq1)
        j = sum([1 for i in range(len(v[0].seqB)) if v[0].seqB[i] == v[0].seqA[i]])*100/len(v[0].seqA)
        amplicon_scores.append(j)
    scores.append(max(amplicon_scores))
    print(len(scores))

print('two closest sequences: ', max(scores))
plt.hist(scores)
plt.title('KDVY difference between tls sequences hist')
plt.xlabel('percent identity to closes sequence')
plt.show()