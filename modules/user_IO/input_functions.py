from modules.ORF.calculating_cai import relative_adaptiveness
from modules.ORF.TAI import TAI
from modules.shared_functions_and_vars import *
import pandas as pd
from modules.logger_factory import LoggerFactory
import operator
import os

logger = LoggerFactory.create_logger("user_input")


# RE model
def find_org_name(gb_file):
    org_name = gb_file.description
    return ' '.join(org_name.split()[:2])


def tai_from_tgcnDB(org_name):
    base_path = os.path.join(os.path.dirname(__file__))
    tgcn_df = pd.read_csv(os.path.join(base_path, 'filtered_tgcn.csv'), index_col=0)
    all_org_tgcn_dict = tgcn_df.T.to_dict('list')

    if org_name in all_org_tgcn_dict.keys():
        tgcn_dict = dict(zip(tgcn_df.columns, all_org_tgcn_dict[org_name] ))
        tai_weights = TAI(tgcn_dict).index
        logger.info(f'tGCN values were found, tAI profile was calculated')
    else:
        tai_weights = {}
        logger.info(f'tGCN values were not found, tAI profile was not calculated')
    return tai_weights

def extract_gene_data(genbank_path, expression_csv_fid = None):
    """
    regorgnize relevant genebank data
    :param genbank_path: cds file path
    :return: a dictionary, where the
    # prom_dic: gene name to prom
    # cds_dict: gene name to cds
    # intergenic_dict: idx is a placing along the genome, and the value is the intergenic sequence
        """

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
                    if 'gene' in feature.qualifiers.keys():
                        name = feature.qualifiers['gene'][0]
                    elif 'locus_tag' in feature.qualifiers.keys():
                        name = feature.qualifiers['locus_tag'][0]
                    else:
                        continue
                    if feature.location is not None \
                            and name not in gene_names:
                        cds = genome[feature.location.start: feature.location.end]
                        function = ' '.join(feature.qualifiers['product'])
                        if feature.location.strand == -1:
                            cds = reverse_complement(cds)

                        if len(cds)%3 !=0:
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
                                estimated_expression[name + '|'+function] = mrna_level
                            except:
                                continue

    entry_num = len(gene_names)
    name_and_function = [gene_names[i] + '|' + functions[i] for i in range(entry_num)]
    # prom200_dict = extract_prom(starts, ends, strands, name_and_function, prom_length=200, genome=genome)  # fix!!
    cds_dict = {name_and_function[i]: cds_seqs[i] for i in range(entry_num)}

    # return prom200_dict, cds_dict, intergenic_dict, estimated_expression
    return cds_dict, estimated_expression


def calculate_cai_weights_for_input(cds_dict, estimated_expression, expression_csv_fid):
    """
    calculates the cai weights- if estimated_expression dictionary has more than 3 times the number of ribosomal genes,
    30% most highly expressed genes will be used as reference set.
    in any other case, ribosomal genes will be used
    """
    ribosomal_proteins = [cds for description, cds in cds_dict.items() if 'ribosom' in description]

    if len(ribosomal_proteins) < 10:
        logger.info('WARNING: less than 10 ribosomal genes were found, this annotation is likely to be of low quality.\n'
                    'All genes will be used as a reference set for CAI calculation, results may be less accurate.')
        cai_weights = relative_adaptiveness(sequences=list(cds_dict.values()))
    else:
        if len(estimated_expression) <len(ribosomal_proteins)*3: #if we found less than 50 expression levels for genes (or no expression csv supplied), the CAI will be used as an estimation
            cai_weights = relative_adaptiveness(ribosomal_proteins)
            if expression_csv_fid is not None:
                logger.info(
                    f'Not enough genes have supplied expression levels, are the gene names the same as the NCBI genbank convention?')
            logger.info('CAI will be calculated from a reference set of ribosomal proteins and used as estimated expression')

        else:
            sorted_estimated_expression = dict(
                sorted(estimated_expression.items(), key=operator.itemgetter(1), reverse=True))
            highly_expressed_names = list(sorted_estimated_expression.keys())[:round(len(sorted_estimated_expression)* 0.3 )]
            highly_expressed_cds_seqs = [cds for description, cds in cds_dict.items() if description in highly_expressed_names]
            cai_weights = relative_adaptiveness(sequences=highly_expressed_cds_seqs)
            logger.info(f'Expression levels were found for {len(estimated_expression)}')
    return  cai_weights
