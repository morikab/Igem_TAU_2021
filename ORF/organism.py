from collections import namedtuple

from Bio.Seq import Seq

from ORF.TAI import TAI

############################
#
#
# Definitions of Organism and Gene classes
#
#
############################

Feature = namedtuple(typename='Feature', field_names='index_name weights ratio')

class Organism(object):

    def __init__(self, name, tai_weights, cai_weights):

        """

        :param name:
        :param genome_path:
        :param mrna_levels_path:
        :param features_to_generate: tuple (feature_name, ratio)
        """

        self.name = name
        print('=========== Defining ' + self.name + ' ===============')


        self.tai_weights = tai_weights
        self.cai_weights = cai_weights
        self.features = []

        self.features.append(Feature(index_name='TAI', weights=self.tai_weights, ratio=0.5))
        self.features.append(Feature(index_name='CAI', weights=self.cai_weights, ratio=0.5))




# --------------------------------------------------------------------------------

class Gene(Seq):

    def __init__(self, dna_seq_record):

        """
        :param dna_seq_path: fasta file - nucleotides sequence
        """

        dna_data = str(dna_seq_record)

        super().__init__(dna_data)
        self.protein_seq = self.translate()