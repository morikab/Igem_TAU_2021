from collections import namedtuple

from Bio import SeqIO
from Bio.Seq import Seq

from Igem_TAU_2021.ORF.CAI import load_genome, CAI
from Igem_TAU_2021.ORF.TAI import TAI

import os
############################
#
#
# Definitions of Organism and Gene classes
#
#
############################

Feature = namedtuple(typename='Feature', field_names='index_name weights ratio')

class Organism(object):
    # todo: implement creating tGCN from genome (.gb), so that no tGCN_path argument is needed
    def __init__(self, name, genome_path, features_to_generate, mrna_levels_path=None, mrna_level_threshold=0.2):

        """

        :param name:
        :param genome_path:
        :param mrna_levels_path:
        :param features_to_generate: tuple (feature_name, ratio)
        """

        self.name = name

        self.genome_path = genome_path
        self.cds_path = os.path.dirname(self.genome_path)+'/'+name+'_cds'
        self.mrna_levels_path = mrna_levels_path

        if self.mrna_levels_path is None:
            genes_HE_path = None
        else:
            genes_HE_path = os.path.dirname(self.mrna_levels_path)+'/'+name+'_highly_expressed_genes'

        load_genome(self.genome_path, self.cds_path, self.mrna_levels_path, mrna_level_threshold,
                    genes_HE_path=genes_HE_path)

        self.features_to_generate = features_to_generate
        self.features = []

        if self.features_to_generate is not None:
            for index, (feature_name, feature_ratio) in enumerate(features_to_generate):
                if feature_name == 'CAI':
                    self.add_cai_index(ratio=feature_ratio)
                elif feature_name == 'TAI':
                    self.add_tai_index(ratio=feature_ratio)

        else:
            self.add_tai_index(ratio=0.5)
            self.add_cai_index(ratio=0.5)


    def add_tai_index(self, ratio):

        tai_index = TAI(self.genome_path, self.cds_path)
        self.features.append(Feature(index_name='TAI', weights=tai_index, ratio=ratio))


    def add_cai_index(self, ratio):

        cai_index = CAI(self.cds_path).index
        self.features.append(Feature(index_name='CAI', weights=cai_index, ratio=ratio))


# --------------------------------------------------------------------------------

class Gene(Seq):

    def __init__(self, dna_seq_path):

        """
        :param dna_seq_path: .txt file - nucleotides sequence
        """

        with open(dna_seq_path, 'r') as f:
            dna_data = f.read()

        super().__init__(dna_data)
        self.protein_seq = self.translate()