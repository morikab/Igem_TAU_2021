from Bio import SeqIO
from Igem_TAU_2021.ORF.calculating_cai import CAI

#functions for input_main code

def fasta_to_dict(fasta_fid):
    fasta_dict = {record.description:str(record.seq) for record in SeqIO.parse(fasta_fid, 'fasta') }
    return fasta_dict

def write_fasta(fid, list_seq, list_name):
    ofile = open(fid + '.fasta', "w+")
    for i in range(len(list_seq)):
        ofile.write(">" + list_name[i] + "\n" + list_seq[i] + "\n")
    ofile.close()

# write ideas for the promoter model
def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    reverse_complement = "".join(complement.get(base, base) for base in reversed(seq))
    return reverse_complement


# RE model
def find_org_name(gb_file):
    org_name = gb_file.description
    return ' '.join(org_name.split()[:2])


# ORI model
def find_tgcn(gb_path):
    #todo: fix this-
    # each genbank file has a different way to represent the tRNAs in it-
    # some do not specify the anticodon (just the amino acid)
    # some have it in the notes and some have an "anticodon" sub-feature to use.
    tgcn_dict = {}
    for record in SeqIO.parse(gb_path, "genbank"):
        for feature in record.features:
            if str(feature.type).lower() == "trna":
                if 'anticodon' in feature.qualifiers.keys():
                    anticodon = feature.qualifiers['anticodon'][0].split(':')[-1][:3]
                elif 'note' in feature.qualifiers.keys():
                    note = ' '.join(feature.qualifiers['note'])
                    anticodon = note[note.find("(") + 1:note.find(")")]
                else:
                    continue
                if anticodon in tgcn_dict.keys():
                    tgcn_dict[anticodon] += 1
                else:
                    tgcn_dict[anticodon] = 1
    return tgcn_dict


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



def extract_highly_expressed_gene_names(expression_estimation, percent_used =1/3):
    exp_list = list(expression_estimation.values())
    exp_list.sort(reverse=True)
    exp_th = exp_list[round(percent_used*len(exp_list))]
    list_of_highly_exp_genes = [gene_name for gene_name in expression_estimation.keys() if expression_estimation[gene_name]>exp_th]
    return list_of_highly_exp_genes


def extract_gene_data(genbank_path):
    """
    regorgnize relevant genebank data
    :param genbank_path: cds file path
    :return: a dictionary, where the
    # prom_dic: gene name to prom
    # cds_dict: gene name to cds
    # intergenic_dict: idx is a placing along the genome, and the value is the intergenic sequence
    # cai_dict: cai score for each cds
        """
    genome = str(SeqIO.read(genbank_path, format='gb').seq)
    cds_seqs = []
    gene_names = []
    functions = []
    starts = []
    ends = []
    strands = []
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

    entry_num = len(gene_names)
    name_and_function = [gene_names[i] + '|' + functions[i] for i in range(entry_num)]
    prom200_dict = extract_prom(starts, ends, strands, name_and_function, prom_length=200, genome=genome)  # fix!!
    prom400_dict = extract_prom(starts, ends, strands, name_and_function, prom_length=400, genome=genome)  # fix!!
    cds_dict = {name_and_function[i]: cds_seqs[i] for i in range(entry_num)}
    intergenic_dict = extract_intergenic(starts, ends, strands, prom_length=400, genome=genome, len_th=30)


    cai_reference_set = [cds for description, cds in cds_dict.items() if 'ribosom' in description]
    cai_scores, cai_weights = CAI(cds_seqs, reference=cai_reference_set)

    ribosomal_proteins = [cds for description, cds in cds_dict.items() if 'ribosom' in description]
    if len(estimated_expression) <len(ribosomal_proteins): #if we found less than 50 expression levels for genes (or no expression csv supplied), the CAI will be used as an estimation
        cai_scores, cai_weights = CAI(cds_seqs, reference=ribosomal_proteins)
        cai_weights = CAI()

        estimated_expression = cai_scores
        if expression_csv_fid is not None:
            logger.info(
                f'Not enough genes have supplied expression levels, are the gene names the same as the NCBI genbank convention?')
        else:
            logger.info('CAI will be used as an estimated expression measurement')

    else:
        sorted_estimated_expression = dict(
            sorted(estimated_expression.items(), key=operator.itemgetter(1), reverse=True))
        highly_expressed_names = list(sorted_estimated_expression.keys())[:round(len(sorted_estimated_expression)* 0.3 )]
        highly_expressed_cds_seqs = [cds for description, cds in cds_dict.items() if description in highly_expressed_names]
        cai_scores, cai_weights = CAI(cds_seqs, reference=highly_expressed_cds_seqs)
        logger.info(f'Expression levels were found for {len(estimated_expression)}')


    cai_dict = {name_and_function[i]: cai_scores[i] for i in range(entry_num)}
    return prom200_dict, prom400_dict, cds_dict, intergenic_dict, cai_dict, cai_weights

