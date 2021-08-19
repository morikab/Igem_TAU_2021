from Bio.SeqUtils import CodonUsage as cu
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# ------------------------------------------------

def load_genome(genome_path, output_path, genes_HE):
    """

    :param genome_path: original genbank file path
    :param output_path: output path (fasta format)
    :param genes_HE: 20% Highly Expressed Genes, optional
    :return:
    The function extracts all the coding sequences from the supplied genome
    and writes them to a file (output path) in the "fasta" format

    """

    cds_records = []
    with open(genome_path) as input_handle:
        for record in SeqIO.parse(input_handle, "genbank"):
            for feature in record.features:
                if feature.type == "CDS":
                    if genes_HE is not None:
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

