from Bio import SeqIO
import pandas as pd
import os


def write_fasta(fid, list_seq, list_name):
    """
    writes fasta file
    :param fid: output path
    :param list_seq: list of sequences to write
    :param list_name: list of fasta headers
    :return:
    """
    ofile = open(fid + '.fasta', "w+")
    for i in range(len(list_seq)):
        ofile.write(">" + str(list_name[i]) + "\n" + str(list_seq[i]) + "\n")
    ofile.close()

def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    reverse_complement = "".join(complement.get(base, base) for base in reversed(seq))
    return reverse_complement


def make_genomic_dataframe(genbank_path, exp_level_csv):
    """
    regorgnize relevant genebank data
    :param genbank_path: cds file path
    :param exp_level_csv: expression levels of all genes- first col is locus tag (not used), seconf is the gene name, third is the exp level
    :return: The function returns a dataframe with the specified columns for all CDSs:
    ['gene_name', 'mrna_level', 'function', 'cds_start', 'cds_stop', 'cds_strand', 'cds_seq']
    """
    genome = SeqIO.read(genbank_path, format='gb')
    cds_data = pd.DataFrame(columns=['gene_name', 'mrna_level', 'function', 'cds_start', 'cds_stop', 'cds_strand', 'cds_seq'])
    mrna_levels = pd.read_csv(exp_level_csv,low_memory=False)
    mrna_levels.drop(inplace=True, index = [index for index in range(len(mrna_levels.mRNA_level.to_list())) if mrna_levels.mRNA_level.to_list()[index] =='N.D.'])
    mrna_names = [i.lower() for i in mrna_levels.gene.to_list()]
    mrna_levels = [float(i) for i in mrna_levels.mRNA_level.to_list()]
    i=0
    with open(genbank_path) as input_handle:
            for record in SeqIO.parse(input_handle, "genbank"):
                for feature in record.features:
                    if feature.type == "CDS" and 'gene' in feature.qualifiers.keys():
                        name= feature.qualifiers['gene'][0].lower()
                        if feature.location is not None \
                            and name in mrna_names \
                            and name not in cds_data.gene_name.to_list():
                            mrna_level = [mrna_levels[index] for index in range(len(mrna_names)) if mrna_names[index]==name][0]
                            cds = genome[feature.location.start: feature.location.end]
                            function = ' '.join(feature.qualifiers['product'])
                            if feature.location.strand ==-1:
                                cds = reverse_complement(cds)
                            gene_info = [name, mrna_level, function,feature.location.start, feature.location.end, feature.location.strand, cds]
                            cds_data.loc[i]=gene_info
                            i+=1
    cds_data.sort_values('mrna_level', inplace=True, ascending=False)
    return cds_data

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
        if strand ==1:
            genome = genome[:start-prom_length] + '-'*(end-start+prom_length) + genome[end:]
        else: #strand ==-1
            genome = genome[:start] + '-' * (end - start + prom_length) + genome[end+prom_length:]
    intergenic_list = genome.split('-')
    return [i for i in intergenic_list if len(i)>len_th]


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
        if strand ==1:
            genome_fwd = genome[:start ] + '-' * (end - start) + genome[end:]
        else:
            genome_rev = genome[:start] + '-' * (end - start) + genome[end:]
    for idx, vals in enumerate(zip(cds_start, cds_strand, cds_names)):
        start, strand, name = vals
        if strand ==1:
            prom = genome_fwd[start-prom_length:start]
        else:
            prom = reverse_complement(genome_rev[start:start+prom_length])
        if prom != '-'*prom_length and len(prom.replace('-',''))>0:
            prom_dict[name] = prom.replace('-','')
    return prom_dict


### code
org_groups = ['ecoli', 'bacillus']
for org in org_groups:

    # input fids
    expression_level_csv_fid = org + '_mrna_level.csv'
    genbank_path = org + '.gb'
    #----------------------------------

    # output fids
    output_path = os.path.join('outputs', org)
    gene_df_csv = os.path.join(output_path, 'cds_data_' + org + '.csv')
    prom200_fasta = os.path.join(output_path, '200bp_promoters_' + org )
    prom200_highest50_fasta =  os.path.join(output_path, '200bp_highest50%_promoters_' + org )
    prom400_fasta = os.path.join(output_path, '400bp_promoters_' + org )
    intergenic_fasta = os.path.join(output_path, 'intergenic_' + org )
    #-----------------------------------

    # expression csv
    cds_df = make_genomic_dataframe(genbank_path, expression_level_csv_fid)
    genome = SeqIO.read(genbank_path, format='gb')
    cds_df.to_csv(gene_df_csv, index= False)
    #-----------------------------------

    # fasta files
    intergenic_sequences_list = extract_intergenic(cds_df.cds_start, cds_df.cds_stop, cds_df.cds_strand, prom_length=200, genome= genome, len_th=8)
    write_fasta(intergenic_fasta, intergenic_sequences_list, [str(i) for i in range(len(intergenic_sequences_list))]) #headers are just rank of placing along the genome

    prom200_dict = extract_prom(cds_df.cds_start, cds_df.cds_stop, cds_df.cds_strand, cds_df.gene_name,
                               prom_length=200, genome= genome)
    prom200_seq_list = list(prom200_dict.values())
    prom200_head_list = list(prom200_dict.keys())
    number_of_highly = round(len(prom200_head_list)/2)
    write_fasta(prom200_fasta, prom200_seq_list, prom200_head_list)
    write_fasta(prom200_highest50_fasta, prom200_seq_list[:number_of_highly], prom200_head_list[:number_of_highly])

    prom400_dict = extract_prom(cds_df.cds_start, cds_df.cds_stop, cds_df.cds_strand, cds_df.gene_name,
                                prom_length=400, genome=genome)
    prom400_seq_list = list(prom200_dict.values())
    prom400_head_list = list(prom200_dict.keys())
    write_fasta(prom400_fasta, prom400_seq_list, prom400_head_list)