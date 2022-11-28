#all info is from ecoli model by rotem, noy and rawan
import pandas as pd
tdr_df = pd.read_excel('MTDR.xlsx')
from modules.ORF.calculating_cai import general_geomean
import matplotlib.pyplot as plt
from  scipy.stats import spearmanr
# print(dict(zip(tdr_df['codon'].to_list(), tdr_df['ecoli'].to_list())))

tai_dict = dict(AAA=0.75, AAC=0.5, AAG=0.23999999999999996, AAT=0.29500000000000004, ACA=0.125, ACC=0.25, ACG=0.29,
                ACT=0.14750000000000002, AGA=0.125, AGC=0.125, AGG=0.16499999999999998, AGT=0.07375000000000001,
                ATA=0.0, ATC=0.375, ATG=1.0, ATT=0.22125000000000003, CAA=0.25, CAC=0.125, CAG=0.32999999999999996,
                CAT=0.07375000000000001, CCA=0.125, CCC=0.125, CCG=0.16499999999999998, CCT=0.07375000000000001,
                CGA=4.999999999999449e-05, CGC=0.36, CGG=0.125, CGT=0.5, CTA=0.125, CTC=0.125, CTG=0.54,
                CTT=0.07375000000000001, GAA=0.5, GAC=0.375, GAG=0.15999999999999998, GAT=0.22125000000000003,
                GCA=0.375, GCC=0.25, GCG=0.11999999999999998, GCT=0.14750000000000002, GGA=0.125, GGC=0.5,
                GGG=0.16499999999999998, GGT=0.29500000000000004, GTA=0.625, GTC=0.25, GTG=0.19999999999999996,
                GTT=0.14750000000000002, TAA=0.0, TAC=0.375, TAG=0.0, TAT=0.22125000000000003, TCA=0.125, TCC=0.25,
                TCG=0.16499999999999998, TCT=0.14750000000000002, TGA=0.125, TGC=0.125, TGG=0.16499999999999998,
                TGT=0.07375000000000001, TTA=0.125, TTC=0.25, TTG=0.16499999999999998, TTT=0.14750000000000002)
cai_dict = dict(AAA=1.0, AAC=1.0, AAG=0.30339566929133854, AAT=0.8071286159012208, ACA=0.29234858918105244, ACC=1.0,
                ACG=0.6127854466258849, ACT=0.3751401210645998, AGA=0.0873249940675955, AGC=1.0,
                AGG=0.048171124444896435, AGT=0.5390529770276605, ATA=0.1352630800335322, ATC=0.828369248976774,
                ATG=1.0, ATT=1.0, CAA=0.5308106141331185, CAC=0.7551246860947265, CAG=1.0, CAT=1.0,
                CCA=0.35896942278692207, CCC=0.23239323643597395, CCG=1.0, CCT=0.29630699136907623,
                CGA=0.15661547849079627, CGC=1.0, CGG=0.24014373368588765, CGT=0.950642394657446,
                CTA=0.07269878905919214, CTC=0.20934597905552263, CTG=1.0, CTT=0.20600107262822143, GAA=1.0,
                GAC=0.5973479415982884, GAG=0.44915318455651126, GAT=1.0, GCA=0.5927597804921225,
                GCC=0.7563949371570189, GCG=1.0, GCT=0.4475792175606302, GGA=0.26135331516802907, GGC=1.0,
                GGG=0.3684276919971743, GGT=0.8295488949439903, GTA=0.41237054577615473, GTC=0.5818378933553964,
                GTG=1.0, GTT=0.6920201991383984, TAA=1.0, TAC=0.7622985481370108, TAG=0.10773176987907658, TAT=1.0,
                TCA=0.436615096108767, TCC=0.5356305672761369, TCG=0.5558368495077356, TCT=0.5219878105954056,
                TGA=0.4477830707218762, TGC=1.0, TGG=1.0, TGT=0.7923738665426645, TTA=0.25874333135743927,
                TTC=0.7437472575691093, TTG=0.25649928021000934, TTT=1.0)
tdr_dict = dict(AAA=0.013, AAC=0.01326, AAG=0.01307, AAT=0.01529, ACA=0.02404, ACC=0.01914, ACG=0.0194, ACT=0.01531,
                AGA=0.03089, AGC=0.01627, AGG=0.02889, AGT=0.01843, ATA=0.01793, ATC=0.01591, ATG=0.01798, ATT=0.01553,
                CAA=0.01792, CAC=0.01413, CAG=0.01693, CAT=0.01534, CCA=0.01929, CCC=0.02652, CCG=0.01551, CCT=0.01858,
                CGA=0.02292, CGC=0.01653, CGG=0.03036, CGT=0.01173, CTA=0.03195, CTC=0.0251, CTG=0.01649, CTT=0.02145,
                GAA=0.0139, GAC=0.01303, GAG=0.01436, GAT=0.01568, GCA=0.01338, GCC=0.02041, GCG=0.01348, GCT=0.01852,
                GGA=0.01938, GGC=0.01575, GGG=0.01792, GGT=0.01221, GTA=0.01778, GTC=0.017, GTG=0.01728, GTT=0.01493,
                TAC=0.01983, TAT=0.01681, TCA=0.01213, TCC=0.01415, TCG=0.01604, TCT=0.01525, TGC=0.01958, TGG=0.02147,
                TGT=0.01871, TTA=0.02154, TTC=0.01557, TTG=0.02046, TTT=0.01748)

sequence_data = pd.read_csv('cds_data_ecoli.csv')

sequences = sequence_data['cds_seq'].to_list()
exp_levels = sequence_data['mrna_level'].to_list()
tai_scores = general_geomean(sequence_lst=sequences, weights= tai_dict)
cai_scores = general_geomean(sequence_lst=sequences, weights= cai_dict)
tdr_scores = general_geomean(sequence_lst=sequences, weights= tdr_dict)

def precision_round(number, digits=3):
    power = "{:e}".format(number).split('e')[1]
    return round(number, -(int(power) - digits))

check_for = 'TDR'

graph_dict = {'CAI': cai_scores, 'TAI':tai_scores, 'TDR':tdr_scores}
color_dict = {'CAI': 'blue', 'TAI':'green', 'TDR':'purple'}

plt.scatter(exp_levels, graph_dict[check_for], edgecolors=color_dict[check_for] ,s=3)
r, pval = spearmanr(exp_levels, graph_dict[check_for])
# r = round(r, 3)
# pval = precision_round(pval, 3)
plt.suptitle(f'mRNA expression levels and {check_for} scores for E. coli')
plt.title(f'Spearman correlation: rho={r}, pval={pval}', fontsize =9)
plt.xlabel('Expression level')
plt.ylabel(f'{check_for} score')
plt.show()




