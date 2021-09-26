from modules.ORF.calculating_cai import relative_adaptiveness, general_geomean
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


# promoter model
def extract_intergenic(cds_start, cds_stop, cds_strand, prom_length, genome, len_th=30):
    """
    extract intergenic sequences for streme intergenic motifs
    :param cds_start:
    :param cds_stop:
    :param cds_strand: -1 if the sequence is on the reverse strand, 1 for forward
    :param prom_length: length of the sequences considered as promoters
    :param genome: string of the whole 5'->3' genome on the fwd strand
    :param len_th:shortest sequence to add to the file
    :return: list of intergenic sequences
    """
    for idx, vals in enumerate(zip(cds_start, cds_stop, cds_strand)):
        start, end, strand = vals
        if strand == 1:
            genome = genome[:start - prom_length] + '-' * (end - start + prom_length) + genome[end:]
        else:  # strand ==-1
            genome = genome[:start] + '-' * (end - start + prom_length) + genome[end + prom_length:]
    intergenic_list = genome.split('-')
    return {k: i for k, i in enumerate(intergenic_list) if len(i) > len_th}


def extract_prom(cds_start, cds_stop, cds_strand, cds_names, prom_length, genome):
    """
    extracts prom_length bases before the cds on the correct strand
    :param cds_start:
    :param cds_stop:
    :param cds_strand: -1 for reverse strand, 1 for forward
    :param cds_names: list of gene names
    :param prom_length: number of bases before cds to take
    :param genome: string of the whole 5'->3' genome on the fwd strand
    :return: dict of promoters {gene_name:promoter}
    """
    prom_dict = {}
    genome_fwd = genome
    genome_rev = genome
    for idx, vals in enumerate(zip(cds_start, cds_stop, cds_strand)):
        start, end, strand = vals
        if strand == 1:
            genome_fwd = genome[:start] + '-' * (end - start) + genome[end:]
        else:
            genome_rev = genome[:start] + '-' * (end - start) + genome[end:]
    for idx, vals in enumerate(zip(cds_start, cds_strand, cds_names)):
        start, strand, name = vals
        if strand == 1:
            prom = genome_fwd[start - prom_length:start]
        else:
            prom = reverse_complement(genome_rev[start:start + prom_length])
        if prom != '-' * prom_length and len(prom.replace('-', '')) > 0:
            prom_dict[name] = prom.replace('-', '')
    return prom_dict



def extract_highly_expressed_promoters(expression_estimation, prom_dict, percent_used =1/3):
    exp_list = list(expression_estimation.values())
    exp_list.sort(reverse=True)
    exp_th = exp_list[round(percent_used*len(exp_list))]
    print(len(prom_dict))
    print(len(expression_estimation))
    # print(prom_dict['B6A19_RS00005|AAA family ATPase'])
    print(expression_estimation['B6A19_RS00005|AAA family ATPase'])
    highly_exp_prom = {gene_name:prom_dict[gene_name]
                       for gene_name in expression_estimation.keys()
                       if expression_estimation[gene_name]>exp_th
                       }

    return highly_exp_prom


def extract_gene_data(genbank_path, expression_csv_fid=None):
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
            mrna_names =  list(gene_name_to_mrna_level.keys())
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
    prom200_dict = extract_prom(starts, ends, strands, name_and_function, prom_length=200, genome=genome)  # fix!!
    cds_dict = {name_and_function[i]: cds_seqs[i] for i in range(entry_num)}
    intergenic_dict = extract_intergenic(starts, ends, strands, prom_length=2000, genome=genome, len_th=30)

    return prom200_dict, cds_dict, intergenic_dict, estimated_expression



def calculate_cai_weights_for_input (cds_dict, estimated_expression, expression_csv_fid):
    """
    calculates the cai weights- if estimated_expression dictionary has more than 3 times the number of ribosomal genes,
    30% most highly expressed genes will be used as reference set.
    in any other case, ribosomal genes will be used
    """
    ribosomal_proteins = [cds for description, cds in cds_dict.items() if 'ribosom' in description]
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


