from collections import namedtuple

from Bio.Seq import Seq

from Igem_TAU_2021.ORF.TAI import TAI

############################
#
#
# Definitions of Organism and Gene classes
#
#
############################

Feature = namedtuple(typename='Feature', field_names='index_name weights ratio')

class Organism(object):

    def __init__(self, name, gene_cds, tgcn, cai_weights, features_to_generate=None):

        """

        :param name:
        :param genome_path:
        :param mrna_levels_path:
        :param features_to_generate: tuple (feature_name, ratio)
        """

        self.name = name
        print('=========== Defining ' + self.name + ' ===============')

        self.gene_cds = gene_cds
        self.tgcn = tgcn
        self.cai_weights = cai_weights
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

        tai_weights = TAI(self.gene_cds, self.tgcn).index
        self.features.append(Feature(index_name='TAI', weights=tai_weights, ratio=ratio))


    def add_cai_index(self, ratio):

        self.features.append(Feature(index_name='CAI', weights=self.cai_weights, ratio=ratio))


# --------------------------------------------------------------------------------

class Gene(Seq):

    def __init__(self, dna_seq_record):

        """
        :param dna_seq_path: fasta file - nucleotides sequence
        """

        dna_data = str(dna_seq_record.seq)

        super().__init__(dna_data)
        self.protein_seq = self.translate()