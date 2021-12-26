import typing
from modules.user_IO.input_functions import *
from modules.ORF.TAI import TAI
from modules.ORF.calculating_cai import general_geomean
from modules.logger_factory import LoggerFactory
from statistics import mean, stdev

# initialize the logger object
logger = LoggerFactory.create_logger("user_input")


class UserInputModule(object):
    @staticmethod
    def get_name() -> str:
        return "User Input"

    @classmethod
    def run_module(cls, user_inp_raw: typing.Dict) -> typing.Dict:
        logger.info('##########################')
        logger.info('# USER INPUT INFORMATION #')
        logger.info('##########################')
        return cls._parse_input(user_inp_raw)

    @classmethod
    def _parse_input(cls, usr_inp):
        """
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
            @organism: contains the org names as keys, and for each organism- the key is the scientific organism name,
            and the value is _parse_single_input for the organism's input
        """
        full_inp_dict = {}

        organisms_dict = {}
        organisms_names = set()

        # TODO - should we use asyncio here?
        # creating the sub dictionary for each organism where the key is the scientific name
        for key, val in usr_inp['organisms'].items():
            try:
                organism_name, organism_dict = cls._parse_single_input(val)
            except:
                raise ValueError(f'Error in input: {key}, re-check your input')
            if organism_name in organisms_names:
                raise ValueError(f"Organism: {organism_name}'s genome is inserted twice, re-check your input.")
            organisms_dict[organism_name] = organism_dict
            organisms_names.add(organism_name)

        full_inp_dict["organisms"] = organisms_dict

        # TODO - move to another method (Use configuration value - should create plot - for running this
        # #plots for intergenic sequence analysis: change pron length to 0 and th to 0 as well in the intergenic seuqence code before applying
        #     plt.hist([len(int_seq) for int_seq in list(org_dict['intergenic'].values())],
        #              label=org_name, bins=1000, alpha=0.3, density=True)
        # print(sum([len(int_seq) for int_seq in list(org_dict['intergenic'].values())])/len([len(int_seq) for int_seq in list(org_dict['intergenic'].values())]))
        # plt.title('Intergenic sequence length histogram for promoter analysis')
        # plt.xlabel("Intergenic sequence length")
        # plt.xlim(0,1500)
        # plt.legend()
        # plt.show()

        # add non organism-specific keys to dict

        # Read ORF sequence
        orf_fasta_fid = usr_inp['sequence']
        orf_seq = None
        if orf_fasta_fid is not None:
            try:
                orf_seq = str(SeqIO.read(orf_fasta_fid, 'fasta').seq)
            except:
                raise ValueError(
                    f'Error in protein .fasta file: {orf_fasta_fid}, make sure you inserted an undamaged .fasta file containing a single recored')
        logger.info(f'\n\nSequence to be optimized given in the following file {orf_fasta_fid}')
        logger.info(f'containing this sequence: {orf_seq}')

        full_inp_dict['sequence'] = orf_seq
        full_inp_dict['tuning_param'] = usr_inp['tuning_param']
        return full_inp_dict

    @staticmethod
    def _parse_single_input(organism_input):
        """
        create the relevant information for each organism
        :param organism_input: the dictionary supplied for every organism:
            {'genome_path': '.gb', 'genes_HE': '.csv', 'optimized': False }
        :return:
        @org_name = scientific organism name
        @org_dict = {
            'tgcn': {anti codon:number of occurences}, for ORF model
            '200bp_promoters': {gene name and function: prom}, promoter model
            'third_most_HE': same as 200bp_promoters but only third most highly expressed promoters
            'gene_cds':  {gene name and function : cds}, for ORF model
            'intergenic':{position along the genome: intergenic sequence}, promoter model
            'expression_estimation_of_all_genes': {'gene name and function: estimated expression }when the expression csv is not given- the CAI is used as expression levels
            'cai_scores': {'gene_name': cai} ORF and promoter
            'cai_profile': {codon:cai score},  # ORF model
            'tai_scores': {'gene_name': tai} ORF and promoter
            'tai_profile': {codon:cai score},  # ORF model
            'optimized': bool- True if organism is optimized}
        """
        gb_path = organism_input['genome_path']
        exp_csv_fid = organism_input['expression_csv']
        try:
            gb_file = SeqIO.read(gb_path, format='gb')
        except:
            raise ValueError(
                f'Error in genome GenBank file: {gb_path}, make sure you inserted an undamaged .gb file containing '
                f'the full genome sequence and annotations'
            )

        organism_name = find_org_name(gb_file)
        logger.info(f'\nInformation about {organism_name}:')
        if organism_input['optimized']:
            logger.info('Organism is optimized')
        else:
            logger.info('Organism is deoptimized')

        cds_dict, estimated_expression = extract_gene_data(gb_path, exp_csv_fid)
        logger.info(f'Number of genes: {len(cds_dict)}')
        gene_names = list(cds_dict.keys())
        cai_weights = calculate_cai_weights_for_input(cds_dict, estimated_expression, exp_csv_fid)
        cai_scores = general_geomean(sequence_lst=cds_dict.values(), weights=cai_weights)
        cai_scores_dict = {gene_names[i]: cai_scores[i] for i in range(len(gene_names))}

        try:
            tai_weights = TAI(tai_from_tgcnDB(organism_name)).index
            tai_scores = general_geomean(sequence_lst=cds_dict.values(), weights= tai_weights)
            tai_scores_dict = {gene_names[i]: tai_scores[i] for i in range(len(gene_names))}
            tai_mean = mean(tai_scores)
            std_tai = stdev(tai_scores)
        except:
            tai_weights = None
            tai_mean = None
            std_tai = None
            tai_scores_dict = {}

        organism_dict = {
            'cai_profile': cai_weights,  # {dna_codon:cai_score}
            'tai_profile': tai_weights,  # {dna_codon:tai_score}, if not found in tgcnDB it will be an empty dict.
            'cai_scores': cai_scores_dict,  # {'gene_name': score}
            'tai_scores': tai_scores_dict,  # {'gene_name': score}, if not found in tgcnDB it will be an empty dict.
            'optimized': organism_input['optimized'],
            'cai_avg': mean(cai_scores),
            'tai_avg': tai_mean,
            'cai_std': stdev(cai_scores),
            'tai_std': std_tai
            # 'gene_cds': cds_dict,  # cds dict {gene name and function : cds}, for ORF model
            # 'expression_estimation_of_all_genes': estimated_expression, # when the expression csv is not given- the CAI is used as expression levels
        }

        return organism_name, organism_dict
