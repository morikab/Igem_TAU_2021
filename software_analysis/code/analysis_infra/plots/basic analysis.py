from os import listdir
from os.path import isfile, join
import json
import matplotlib.pyplot as plt
import pandas as pd
import math

class obj:
    def __init__(self, dict1):
        self.__dict__.update(dict1)
    def __getattr__(self, attr):
        return self[attr]


def avg(lst):
    return sum(lst)/len(lst)

if __name__ == "__main__":
    tls_match_dir = '../../../data/tls_genome_match/'
    tls_files = [join(tls_match_dir, f) for f in listdir(tls_match_dir) if isfile(join(tls_match_dir, f))
                 and 'json'in f]

    blast_res_df = pd.DataFrame(
        columns=['n_hits','avg_pident', 'avg_match_len', 'avg_reads_per_hit'])

    titles = []
    n_reads_plt = plt
    hits_pident_plt = plt
    for idx, tls_file in enumerate(tls_files):
        print(idx, tls_file)
        tls = json.load(open(tls_file), object_hook=obj)
        titles.append(tls.title)
        n_reads_per_hit = [math.log(genome.scores.counts) for genome in tls.match_data.__dict__.values()]
        hits_pident = [genome.scores.pident for genome in tls.match_data.__dict__.values()]
        avg_reads_per_hit = avg(n_reads_per_hit)
        avg_pident = avg(hits_pident)
        blast_res_df.loc[idx+1] = [len(tls.match_data.__dict__), avg_pident, tls.avg_match_len, avg_reads_per_hit]

        #hits_pident_plt.hist(hits_pident, bins = 100, alpha = 0.3)
        n_reads_plt.hist(n_reads_per_hit, bins = 100, alpha = 0.3)

    blast_res_df['titles'] = titles

    # hits_pident_plt.xlabel('highest pident for hit')
    # hits_pident_plt.title('pident hist')
    # plt.show()
    # hits_pident_plt.savefig('pident_hist.png')

    n_reads_plt.xlabel('ln(n_reads_per_hit)')
    n_reads_plt.title('n_reads_per_hit hist')
    plt.show()
    n_reads_plt.savefig('ln_n_reads_per_hit_hist.png')

    corr_df = blast_res_df.corr(method='spearman')
    corr_df.to_csv('spearman_correlations_for_res.df')
    blast_res_df.to_csv('blast_job_summary')
    print(corr_df)

    plt.scatter(blast_res_df['n_hits'], blast_res_df['avg_pident'], n)
    plt.xlabel('n_hits')
    plt.ylabel('avg_pident')
    plt.show()
    plt.savefig('fig1.png')

    plt.scatter(blast_res_df['n_hits'], blast_res_df['avg_match_len'])
    plt.xlabel('n_hits')
    plt.ylabel('avg_match_len')
    plt.show()
    plt.savefig('fig2.png')

    plt.scatter(blast_res_df['avg_pident'], blast_res_df['avg_reads_per_hit'])
    plt.xlabel('avg_pident')
    plt.ylabel('avg_reads_per_hit')
    plt.show()
    plt.savefig('fig3.png')

    plt.scatter(blast_res_df['avg_match_len'], blast_res_df['avg_reads_per_hit'])
    plt.xlabel('avg_match_len')
    plt.ylabel('avg_reads_per_hit')
    plt.show()
    plt.savefig('fig4.png')



v=1
