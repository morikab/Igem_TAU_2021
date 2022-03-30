from modules.logger_factory import LoggerFactory
from modules import models
from modules.sequence_family.clustering_optimization import *


class SequenceFamilyModule(object):
    @staticmethod
    def run_module(user_input: models.UserInput):
        clustering_mat, opt_org_list = dict_to_cluster_np_array(user_input)
        #todo: decide minimal val to devide into 2 sequences
        best_clustering = create_n_clusters(clustering_mat, n_clus=user_input.clusters_count)
        object_list = return_list_of_sub_microbiomes(best_clustering, user_input)
        return object_list
