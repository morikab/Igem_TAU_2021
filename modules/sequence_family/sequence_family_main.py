from modules.logger_factory import LoggerFactory
from modules import models
from modules.sequence_family.clustering_optimization import *


class SequenceFamilyModule(object):
    @staticmethod
    def run_module(user_input: models.UserInput ,c_index = 'dbi', c_method = 'kmeans'):
        clustering_mat, opt_org_list= dict_to_cluster_np_array(user_input)
        #todo: decide minimal val to devide into 2 sequences
        best_clusturing, best_score = find_best_clustering(clustering_mat,
                                                           max_clus_num=4,
                                                           c_index=c_index,
                                                           c_method=c_method)
        object_list = return_list_of_sub_microbiomes(best_clusturing, user_input)
        return object_list
