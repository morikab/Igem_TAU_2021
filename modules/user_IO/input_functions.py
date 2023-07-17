import os

import operator
import typing

from Bio import SeqIO
import pandas as pd

from logger_factory.logger_factory import LoggerFactory
from modules.configuration import Configuration
from modules.ORF.calculating_cai import relative_adaptiveness
from modules.ORF.TAI import TAI
from modules.shared_functions_and_vars import *


logger = LoggerFactory.get_logger()
config = Configuration.get_config()


# RE model
def find_org_name(gb_file):
    org_name = gb_file.description
    return ' '.join(org_name.split()[:2])


def tai_from_tgcnDB(org_name):
    base_path = os.path.join(os.path.dirname(__file__))
    tgcn_df = pd.read_csv(os.path.join(base_path, 'filtered_tgcn.csv'), index_col=0)
    all_org_tgcn_dict = tgcn_df.T.to_dict('list')

    if org_name in all_org_tgcn_dict.keys():
        tgcn_dict = dict(zip(tgcn_df.columns, all_org_tgcn_dict[org_name]))
        tai_weights = TAI(tgcn_dict).index
        logger.info(f'tGCN values were found, tAI profile was calculated')
    else:
        tai_weights = {}
        logger.info(f'tGCN values were not found, tAI profile was not calculated')
    return tai_weights


def extract_gene_data(genbank_path: str, expression_csv_fid=None):
    genome = str(SeqIO.read(genbank_path, format='gb').seq)
    cds_seqs = []
    gene_names = []
    functions = []
    starts = []
    ends = []
    strands = []
    estimated_expression = {}

    if expression_csv_fid is not None:
        try:
            expression_df = pd.read_csv(expression_csv_fid)

            gene_name_to_mrna_level = {}
            for idx, pair in enumerate(zip(expression_df.gene.to_list(), expression_df.mRNA_level.to_list())):
                measured_gene_name, expression_level = pair
                try:
                    gene_name_to_mrna_level[measured_gene_name.lower()] = float(expression_level)
                except:
                    continue
            mrna_levels = list(gene_name_to_mrna_level.values())
            mrna_names = list(gene_name_to_mrna_level.keys())
        except:
            expression_csv_fid = None
            logger.info('Expression data file is corrupt. \nMake sure that: ')
            logger.info('1. File is in csv format')
            logger.info('2. Gene names fit their NCBI naming ')
            logger.info('3. Column with the gene names is labeled "gene"  ')
            logger.info('4. Column with the gene expression levels is labeled "mRNA_level" ')

    with open(genbank_path) as input_handle:
        for record in SeqIO.parse(input_handle, "genbank"):
            for feature in record.features:
                if feature.type == "CDS":
                    if "gene" in feature.qualifiers.keys():
                        name = feature.qualifiers['gene'][0]
                    elif "locus_tag" in feature.qualifiers.keys():
                        name = feature.qualifiers["locus_tag"][0]
                    else:
                        continue
                    if feature.location is not None and name not in gene_names:
                        cds = genome[feature.location.start: feature.location.end]
                        function = " ".join(feature.qualifiers["product"])
                        if feature.location.strand == -1:
                            cds = reverse_complement(cds)

                        if len(cds) % 3 != 0:
                            continue

                        gene_names.append(name)
                        cds_seqs.append(cds)
                        functions.append(function)
                        starts.append(feature.location.start)
                        ends.append(feature.location.end)
                        strands.append(feature.location.strand)

                        if expression_csv_fid is not None:
                            try:
                                mrna_level = [mrna_levels[index] for index in range(len(mrna_names)) if mrna_names[index] == name.lower()][0]
                                estimated_expression[name + '|'+ function] = mrna_level
                            except:
                                continue

    entry_num = len(gene_names)
    name_and_function = [gene_names[i] + '|' + functions[i] for i in range(entry_num)]
    # prom200_dict = extract_prom(starts, ends, strands, name_and_function, prom_length=200, genome=genome)  # fix!!
    cds_dict = {name_and_function[i]: cds_seqs[i] for i in range(entry_num)}

    return cds_dict, estimated_expression


def calculate_cai_weights_for_input(
        cds_dict: typing.Dict[str, typing.Any],
        estimated_expression_dict: typing.Dict[str, float],
) -> typing.Tuple[typing.Dict[str, float], typing.Sequence[str]]:
    """
    calculates the cai weights - if estimated_expression dictionary has more than 3 times the number of ribosomal genes,
    30% most highly expressed genes will be used as reference set.
    in any other case, ribosomal genes will be used
    """
    ribosomal_proteins_count_threshold = config["INPUT"]["RIBOSOMAL_PROTEINS_COUNT_THRESHOLD"]
    ribosomal_proteins = {description: cds for description, cds in cds_dict.items() if "ribosom" in description}
    ribosomal_proteins_count = len(ribosomal_proteins)
    logger.info(F"Found {ribosomal_proteins_count} ribosomal proteins in input genome.")

    reference_genes = []

    if len(estimated_expression_dict) < max(ribosomal_proteins_count, ribosomal_proteins_count_threshold) * 3:
        logger.info("Estimated expression dictionary does not have enough expression levels. CAI will be calculated "
                    "from a reference set of ribosomal proteins or the entire genome.")

        if ribosomal_proteins_count < ribosomal_proteins_count_threshold:
            logger.warning(F"Less than {ribosomal_proteins_count_threshold} ribosomal genes were found. This "
                           F"annotation is likely to be of low quality and therefore results may be less accurate."
                           F"CAI will be calculated from the entire genome as a reference set.")
            cai_weights = relative_adaptiveness(sequences=list(cds_dict.values()))
        else:
            logger.info("CAI will be calculated from a reference set of ribosomal proteins.")
            cai_weights = relative_adaptiveness(ribosomal_proteins.values())
            reference_genes = list(ribosomal_proteins.keys())
    else:
        logger.info("CAI will be calculated from a reference set of estimated expression dictionary.")
        logger.info(F"Expression levels were found for {len(estimated_expression_dict)}")
        estimated_expression_threshold = config["INPUT"]["EXPRESSION_PERCENTAGE_THRESHOLD"]
        sorted_estimated_expression = dict(
            sorted(estimated_expression_dict.items(), key=operator.itemgetter(1), reverse=True)
        )
        highly_expressed_genes_count = round(len(sorted_estimated_expression) * estimated_expression_threshold)
        logger.info(F"Calculate CAI weights from a reference set of {highly_expressed_genes_count} highly expressed "
                    F"genes from estimated expression dictionary.")
        highly_expressed_names = list(sorted_estimated_expression.keys())[:highly_expressed_genes_count]
        highly_expressed_cds_seqs = [cds for description, cds in cds_dict.items() if description in
                                     highly_expressed_names]
        cai_weights = relative_adaptiveness(sequences=highly_expressed_cds_seqs)

    return cai_weights, reference_genes
