import json
import numpy as np
from sklearn.metrics import silhouette_score, davies_bouldin_score
from sklearn.cluster import AgglomerativeClustering, KMeans

# TODO: try different clustering methods, and cluster eval methods- DBI, Silhouette- find min val
# TODO: write efficiently



def dict_to_cluster_np_array(genome_inp_dict):
    clustering_mat=[]
    for cub_data in genome_inp_dict.values():
        cub_profile = cub_data['cai_profile']
        clustering_mat.append(list(cub_profile.values()))
    clustering_mat=np.array(clustering_mat)
    return clustering_mat

#sklearn
def find_best_clustering(clustering_mat, max_clus_num = 10):
    # dist_metric = 'euclidean'
    scores = []
    clusters = []
    for n_clus in range(2, min(np.size(clustering_mat), max_clus_num)):

        ##### clustering options ######
        # clustering = AgglomerativeClustering(n_clusters=n_clus, affinity= dist_metric).fit(clustering_mat)
        clustering = KMeans(n_clusters=n_clus).fit(clustering_mat)

        labels = clustering.labels_

        ##### cluster eval indexes ####
        # score = silhouette_score(clustering_mat, labels, metric= dist_metric)
        score = davies_bouldin_score(clustering_mat, labels)

        scores.append(score)
        clusters.append(labels)

    best_score = min(scores)
    best_clusturing = clusters[scores.index(best_score)]
    return best_clusturing, best_score

f = open('analysis_data.json')
genome_inp_dict = json.load(f)
clustering_mat = dict_to_cluster_np_array(genome_inp_dict)
best_clusturing, best_score = find_best_clustering(clustering_mat, max_clus_num = 10)
print(best_clusturing, best_score)