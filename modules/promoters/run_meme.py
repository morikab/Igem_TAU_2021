import glob
from modules.promoters.globals_and_shared_methods import *


"""
Calls STREME to find motifs that are enriched in one set of promoters and aren't present in the other

@param primary_path: path to a .fasta file with promoter sequences, in which you want to find enriched motifs
@param control_path: path to a .fasta file with promoter sequences, in which you don't want the motifs to occur
@param kmer: length of subsequences that will be created from the input for the first step of the algorithm (default is 3)
@param output_path: name of the folder where output files will be stored
@param minw: minimal width of motifs
@param maxw: maximal width of motifs
@param pvt: p-value threshold for when to stop searching (default is 0.05)
@param ab: alphabet of the sequences (default is dna)
"""
def streme(primary_path, control_path="", kmer=3, output_path=".", minw=8, maxw=15, pvt=0.05, ab="dna"):
    if len(control_path) > 0:
        control_path = "--n " + control_path + " "
    primary_path = "--p " + primary_path + " "
    kmer = "--order " + str(kmer - 1) + " "
    output_path = "--oc " + output_path + " "
    minw = "--minw " + str(minw) + " "
    maxw = "--maxw " + str(maxw) + " "
    pvt = "--pvt " + str(pvt)
    ab = "--" + ab + " "
    command = "streme --verbosity 1 " + output_path + ab + primary_path + control_path + kmer + minw + maxw + pvt
    os.system(command)


"""
Calls MAST to rank the presence of motifs in a sequence set

@param motif_path: path to a .txt or .html file with motif data - output of STREME in find_motifs()
@param seq_path: path to a .fasta file with promoter sequences
@param output_path: name of the folder where output files will be stored
@param ev: e-value threshold for when to stop ranking (default is 10)
"""
def mast(motif_path, seq_path, output_path=".", ev=10):
    command = "mast -oc " + output_path + " -nostatus -norc -remcorr -ev " + str(ev) + " " + motif_path + " " + seq_path
    os.system(command)


"""
Executes one run of STREME for a specific set of files

@param name1: name of the primary file
@param start1: path to the directory containing the file name1
@param name2: name of the control file
@param start2: path to the directory containing the file name2
@param out_path: (optional) directory in which to store output files from STREME
"""
def one_streme(name1, start1, name2, start2, out_path=""):
    k = 3
    minw = 6
    maxw = 20
    
    path1 = os.path.join(start1, name1 + end)
    path2 = os.path.join(start2, name2 + end)
    out_dir = '_'.join([name1, name2])
    if len(out_path) > 0:
        out_dir = os.path.join(out_path, out_dir)
    streme(primary_path=path1, control_path=path2, kmer=k, output_path=out_dir, minw=minw, maxw=maxw)


"""
Runs STREME several times to create all necessary files for subsequent motif ranking and filtering
"""
def run_streme():
    out_path = 'streme_outputs'
    create_folder(os.path.join(start, out_path))
    for opt_org in glob.glob(os.path.join(opt_path, "*", "")):

        #first run: opt organism vs. intergenic
        org_name = opt_org.split(os.sep)[-2]
        dir_name = os.path.join(opt_path, org_name)
        name1 = '_'.join([org_name, '100', '200'])
        name2 = '_'.join([org_name, 'inter'])
        one_streme(name1, dir_name, name2, dir_name, out_path)

        #second run: opt organism HE vs. all deopt organisms HE
        for deopt_org in glob.glob(os.path.join(deopt_path, "*", "")):
             de_org_name = deopt_org.split(os.sep)[-2]
             de_dir_name = os.path.join(deopt_path, de_org_name)
             name1 = '_'.join([org_name, '33', '200'])
             name2 = '_'.join([de_org_name, '33', '200'])
             one_streme(name1, dir_name, name2, de_dir_name, out_path)


"""
Runs MAST for modified motif file and copies HTML output for user

@param motif_path: path to the file containing motifs (already filtered for most selective and also most transcription inducing)
@param promoter_path200: path to the file containing promoter sequences of length 200 (endogenic by default, otherwise chosen by user)
@param promoter_path400: path to the file containing promoter sequences of length 400 (only if promoter_path200 is not chosen by user)

@return: name of MAST output directory
"""
def run_mast(motif_path, promoter_path):
    motif_name = motif_path.split(os.sep)[-1].split('.')[0] #modified motif file- not in original folder anymore (get only name of file without ext or path)
    promoter_name = promoter_path.split(os.sep)[-1].split('.')[0] #to get only file name without ext or path
    out_name = '_'.join(['motif', motif_name, 'seq', promoter_name])
    out_path = os.path.join(start, out_name)
    mast(motif_path, promoter_path, output_path=out_path)

    source_path = os.path.join(out_path, 'mast.html')
    base_path = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
    target_path = os.path.join(base_path, 'logs', 'mast_results.html')
    shutil.copyfile(source_path, target_path)
                                
    return out_path
    
