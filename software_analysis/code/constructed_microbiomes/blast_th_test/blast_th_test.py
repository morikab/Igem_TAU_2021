import pandas as pd
import matplotlib.pyplot as plt
import math

def for_df(lst1):
    lst2 = []
    i=0
    for idx,key in enumerate(lst1):
        print(key)
        if lst1[idx-1] == key:
            i+=1
        else:
            i=0
        print(i)
        lst2.append(i)
    # print(lst1[:20])
    # print(lst2[:20])
    return lst2

def generate_plot(ssqids, pidents):
    read_idx = 0
    check_pair = False
    second_best_score_diff= []
    second_best_hit_place = []
    for idx, pair in enumerate(zip(ssqids, pidents)):
        read, score = pair
        if ssqids[idx-1] != read:
            print(read)
            check_pair=True
            read_idx = 0
        else:
            read_idx +=1
            if score<pidents[idx-1] and check_pair:
                second_best_score_diff.append(pidents[idx-1]-score)
                second_best_hit_place.append(read_idx)
                print('***', read_idx,pidents[idx-1], score)
                check_pair = False

    plt.scatter(second_best_hit_place, second_best_score_diff,s =0.5)
    plt.xlabel('number of hits with identical best score')
    plt.ylabel('diff between best and second best score')
    plt.show()

    plt.hist(second_best_hit_place, bins = 500)
    plt.title('number of hits with identical best score')
    plt.show()

    plt.hist(second_best_score_diff, bins = 500)
    plt.title('diff between best and second best score')
    plt.show()




df = pd.read_csv('tls_500_hits.csv')
blast_col = ['sseqid', 'qseqid', 'pident', 'length',
             'mismatch', 'gapopen', 'qstart', 'qend',
             'sstart', 'send', 'evalue', 'bitscore']


df.columns = blast_col
df.sort_values(['sseqid', 'evalue'], inplace=True, ascending=False)

# #n_genomes per read plot
# df_counts = df['sseqid'].value_counts()
# plt.hist(df_counts, bins=100)
# plt.title('Number of genomes per read')
# plt.xlabel('Number of mapped genomes')
# plt.ylabel('number of reads with corresponding number of mapped genomes')
# plt.show()


generate_plot(df['sseqid'].to_list(), [math.log(i) for i in df['evalue'].to_list()])




