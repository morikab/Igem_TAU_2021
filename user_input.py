from Bio import SeqIO
from ORF.calculating_cai import CAI
import os

### stable version!



print('##########################')
print('# USER INPUT INFORMATION #')
print('##########################')


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
    tgcn_dict = {}
    for record in SeqIO.parse(gb_path, "genbank"):
        for feature in record.features:
            if str(feature.type).lower() == "trna" and 'note' in feature.qualifiers.keys():
                note = ' '.join(feature.qualifiers['note'])

                anticodon = note[note.find("(") + 1:note.find(")")]
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
    cai_dict = {name_and_function[i]: cai_scores[i] for i in range(entry_num)}
    return prom200_dict, prom400_dict, cds_dict, intergenic_dict, cai_dict, cai_weights



def parse_single_input(val):
    '''
    create the relevant information for each organism
    :param val: the dictionary supplied for every organism:
    {'genome_path': '.gb'
                'genes_HE': '.csv'
                'optimized': False }
    :return:
    @org_name = scientific organism name
    @org_dict = {
        'tgcn': {anti codon:number of occurences}, for ORF model
        '200bp_promoters': {gene name and function: prom}, promoter model
        'third_most_HE': same as 200bp_promoters but only third most highly expressed promoters
        'gene_cds':  {gene name and function : cds}, for ORF model
        'intergenic':{position along the genome: intergenic sequence}, promoter model
        'expression_estimation_of_all_genes': {'gene name and function: estimated expression }when the expression csv is not given- the CAI is used as expression levels
        'CAI_score_of_all_genes': {'gene_name': cai} ORF and promoter
        'cai_profile': {codon:cai score},  # ORF model
        'optimized': bool- True if organism is optimized}
    '''
    gb_path = val['genome_path']
    gb_file = SeqIO.read(gb_path, format='gb')
    org_name = find_org_name(gb_file)
    print(f'\nInformation about {org_name}:')
    if  val['optimized']:
        print('Organism is optimized')
    else:
        print('Organism is deoptimized')
    tgcn_dict = find_tgcn(gb_path)
    print(f'Number of tRNA genes found: {sum(tgcn_dict.values())}, for {len(tgcn_dict)} anticodons out of 61')
    prom200_dict, prom400_dict, cds_dict, intergenic_dict, cai_dict, cai_weights = extract_gene_data(gb_path)
    print(f'Number of genes: {len(cds_dict)}, number of intregenic regions: {len(intergenic_dict)}')
    highly_expressed_gene_name_list = extract_highly_expressed_gene_names(cai_dict)

    org_dict = {
        'tgcn': tgcn_dict,  # tgcn dict {codon:number of occurences} for ORF model
        '200bp_promoters': prom200_dict,  # prom_dict {gene name and function: prom}, promoter model
        'third_most_HE': {name: seq for name, seq in prom200_dict.items()
                      if name in highly_expressed_gene_name_list},
        # '400bp_promoters': prom400_dict,  # prom_dict {gene name and function: prom}, promoter model
        'gene_cds': cds_dict,  # cds dict {gene name and function : cds}, for ORF model
        'intergenic': intergenic_dict,  # intergenic dict {position along the genome: intergenic sequence}, promoter model
        'expression_estimation_of_all_genes': cai_dict,  # when the expression csv is not given- the CAI is used as expression levels
        'CAI_score_of_all_genes': cai_dict,  # {'gene_name': expression} ORF and promoter
        'cai_profile': cai_weights,  # ORF model
        'optimized': val['optimized']}
    return org_name, org_dict



# final function
def parse_input(usr_inp):
    '''
    :param usr_inp: in the following format
        {
    'ecoli': {'genome_path': '.gb' - mandatory
                'genes_HE': '.csv'
                'optimized': True - mandatory
                }
     'bacillus': {'genome_path': '.gb'
                'genes_HE': '.csv'
                'optimized': False
                }
        }
    :return: an extended dictionary with the following items:
    @selected_prom : final used list of promoters for MAST
    @sequence : the ORF to optimize
    and for each organism- the key is the scientific organism name, and the value is parse_single_input for the organism's
    input
    '''
    full_inp_dict = {}
    for key, val in usr_inp.items() :
        if key in ['sequence', 'selected_promoters']:
            continue
        org_name, org_dict = parse_single_input(val)
        full_inp_dict[org_name] = org_dict # creating the sub dictionary for each organism- where the key is the scientific name and the value is the following dict:


    #add non org specific keys to dict
    selected_prom = {}
    if usr_inp['selected_promoters'] is not None:  #
        prom_fasta_fid = usr_inp['selected_promoters']
        print(f'promoter options ranked are given in the following file {prom_fasta_fid}')
        selected_prom = fasta_to_dict(prom_fasta_fid)
    else:
        for org, org_dict in full_inp_dict.items():
            for prom_name, prom_Seq in org_dict['third_most_HE'].items():
                selected_prom[prom_name + ' from organism: ' + org] = prom_Seq
    full_inp_dict['selected_prom'] = selected_prom
    full_inp_dict['sequence'] = usr_inp['sequence']
    return full_inp_dict


# input: from yarin
base_path = os.path.join(os.path.dirname(__file__), 'example_data')

user_inp1_raw = {
    'sequence':SeqIO.read(os.path.join(base_path, 'mCherry_original.fasta'), "fasta"),
    'selected_promoters': os.path.join(base_path, 'ORF_optimized_sequences.fasta'), #or None, # or a fasta file of promoter name and promoter
    'opt1': {'genome_path': os.path.join(base_path, 'Escherichia coli.gb'),
             'expression_csv': None, #todo: add treatment for expression data
             'optimized': True},
     'deopt1': {'genome_path': os.path.join(base_path, 'Bacillus subtilis.gb'),
                'expression_csv': None,
               'optimized': False},
     'opt2': {'genome_path': os.path.join(base_path, 'Sulfolobus acidocaldarius.gb'),
              'expression_csv': None,
             'optimized': True},
    'deopt2': {'genome_path': os.path.join(base_path, 'Pseudomonas aeruginosa.gb'),
               'expression_csv': None,
               'optimized': False}}

user_inp = parse_input(user_inp1_raw)
