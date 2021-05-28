import pandas as pd
import assets.gen_dssp as gd
import assets.alphabet_reduction as ar
import assets.gen_Ngrams as gn
import assets.calc_AA as aa
import assets.tesselate as wt
from Bio import BiopythonWarning
from Bio import SeqIO
import warnings
from os import mkdir

def process_input(chain, start, end, alpha, filename, fasta, structure):
    # Process incoming data
    if filename != None:
        error, fasta = gd.getFastaFromPDB(filename = filename, chain = chain)
        structure = gd.getDSSPfromPDB(filename = filename, chain = chain)
    red_seq = ar.reduce_dataset_st(fasta[start:end])
    #generate algorithm data
    if red_seq:
        seq_n_grams = gn.Generate_seq_ngrams(red_seq, alpha)
    if structure:
        struct_n_grams = gn.Generate_struct_ngrams(structure)
    else:
        struct_n_grams = pd.DataFrame()
    if filename:
        tess = wt.tesselate_st(filename, fasta, chain, alpha, start, end)
    else:
        tess = pd.DataFrame()
    if fasta:
        aa_content = aa.calc_AA_st(fasta)
    result = pd.concat([seq_n_grams, struct_n_grams, tess, aa_content], axis=1)
    return(result)



