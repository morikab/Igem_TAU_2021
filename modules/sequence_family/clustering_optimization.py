import json
import numpy as np
from modules import models
from sklearn.metrics import silhouette_score, davies_bouldin_score
from sklearn.cluster import AgglomerativeClustering, KMeans
from modules.logger_factory import LoggerFactory


logger = LoggerFactory.create_logger("user_input") #todo: is there anything here we want ot add to the logger?


# TODO: try different clustering methods, and cluster eval methods- DBI, Silhouette- find min val
# TODO: write efficiently



def dict_to_cluster_np_array(user_input: models.UserInput): #todo: make this work on tai as well
    clustering_mat=[]
    opt_org_list = []
    for organism in user_input.organisms:
        if organism.is_optimized:
            cai_profile = organism.cai_profile
            clustering_mat.append(list(cai_profile.values()))
            opt_org_list.append(organism.name)
    clustering_mat=np.array(clustering_mat)
    return clustering_mat, opt_org_list

#sklearn
def find_best_clustering(clustering_mat, max_clus_num, c_index = 'dbi', c_method = 'kmeans' ):
    dist_metric = 'euclidean'
    scores = []
    clusters = []
    for n_clus in range(2, min(np.size(clustering_mat), max_clus_num)):

        ##### clustering options ######
        if c_method == 'kmeans':
            clustering = KMeans(n_clusters=n_clus).fit(clustering_mat)
        else:
            clustering = AgglomerativeClustering(n_clusters=n_clus, affinity= dist_metric).fit(clustering_mat)

        labels = clustering.labels_

        ##### cluster eval indexes ####
        if c_index == 'dbi':
            score = davies_bouldin_score(clustering_mat, labels)
        else:
            score = silhouette_score(clustering_mat, labels, metric= dist_metric)

        scores.append(score)
        clusters.append(labels)

    best_score = min(scores)
    best_clusturing = clusters[scores.index(best_score)]
    return best_clusturing, best_score


def return_list_of_sub_microbiomes(best_clusturing:list, user_input:models.UserInput):
    opt_org_list = [org.name for org in user_input.organisms if org.is_optimized]
    deopt_org_list =  [org.name for org in user_input.organisms if not org.is_optimized]
    c_assignment_dict = dict(zip(opt_org_list, best_clusturing))
    opt_org_clusters = list({n: [k for k in c_assignment_dict.keys() if c_assignment_dict[k] == n]
         for n in set(c_assignment_dict.values())}.values())

    inp_obj_list = []
    opt_and_deopt = []
    for c_opt_org_list in opt_org_clusters:
        opt_and_deopt = c_opt_org_list+deopt_org_list
        new_user_input = models.UserInput(organisms=[],
                                sequence=user_input.sequence,
                                tuning_parameter=user_input.tuning_parameter)

        new_user_input.organisms =  [user_input.organisms[i] for i in range(len(opt_org_list+deopt_org_list))
                                    if user_input.organisms[i].name in opt_and_deopt]
        inp_obj_list.append(new_user_input)
    return inp_obj_list
