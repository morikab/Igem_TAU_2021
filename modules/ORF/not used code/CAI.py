from Bio.SeqUtils import CodonUsage as cu
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd
import pickle
import os

# ------------------------------------------------

def load_genome(genome_path, output_path, mrna_levels_path, mrna_level_threshold, genes_HE_path):
    """

    :param genome_path: original genbank file path
    :param output_path: output path (fasta format)
    :param mrna_levels_path: .xlsx file: first column - gene name, second column - mrna level. If None, take all genes
    :param mrna_level_threshold:
    :param genes_HE_path:
    :return:
    The function extracts all the coding sequences from the supplied genome
    and writes them to a file (output path) in the "fasta" format

    """

    cds_records = []
    genes_HE = extract_highly_expressed_genes(mrna_levels_path, mrna_level_threshold, genes_HE_path)
    with open(genome_path) as input_handle:
        for record in SeqIO.parse(input_handle, "genbank"):
            for feature in record.features:
                if feature.type == "CDS":
                    if mrna_levels_path is not None:
                        if feature.qualifiers['gene'][0] not in genes_HE:
                            continue
                    if 'protein_id' in feature.qualifiers.keys() and 'gene' in feature.qualifiers.keys() and 'product' in feature.qualifiers.keys():
                        cds_record = SeqRecord(
                            feature.extract(record.seq),
                            id=feature.qualifiers['protein_id'][0],
                            name=feature.qualifiers['gene'][0],
                            description=feature.qualifiers['product'][0]
                        )

                        cds_records.append(cds_record)

    with open(output_path, "w") as output_handle:
        SeqIO.write(cds_records, output_handle, "fasta")


# ------------------------------------------------

def extract_highly_expressed_genes(mrna_levels_path, mrna_level_threshold, genes_HE_path):
    """

    :param mrna_levels_path: .xlsx file: first column - gene name, second column - mrna level.
    :param mrna_level_threshold:
    :param genes_HE_path:
    :return:
    """

    mrna = pd.read_excel(mrna_levels_path)
    mrna = mrna[pd.to_numeric(mrna['mRNA_level'], errors='coerce').notnull()]
    mrna = mrna.sort_values('mRNA_level', ascending=True)
    genes_HE = list(mrna['gene'][0:int(len(mrna['mRNA_level']) * mrna_level_threshold)])

    if genes_HE_path is not None:

        with open(genes_HE_path, 'wb') as handle:
            pickle.dump(genes_HE, handle)

    return genes_HE


# ------------------------------------------------

def generate_cai_index(cds_path):
    cai_object = cu.CodonAdaptationIndex()
    cai_object.generate_index(cds_path)
    cai_index = cai_object.index
    return cai_index

# ------------------------------------------------

class CAI(cu.CodonAdaptationIndex):

    def __init__(self, cds_path):
        super().__init__()
        self.generate_index(cds_path)

    def __call__(self, sequence):
        """

        :param sequence: sequence as a Seq object
        :return:
        """
        return self.cai_for_gene(str(sequence))

