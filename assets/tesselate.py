import math
import re
import time
from collections import Counter
from os import listdir, mkdir
from os.path import isfile, join

import numpy as np
import pandas as pd
import requests
from Bio.PDB import PDBParser
from Bio.SeqUtils import IUPACData, seq1
from pyhull import qconvex, qdelaunay, qvoronoi
from pyhull.convex_hull import ConvexHull
from pyhull.delaunay import DelaunayTri
from scipy.spatial import distance


def reduce_simplex(alpha, tess):
    scheme_A_map = {"C": "B", "M": "B", "F": "B", "I": "B", "L": "B", "V": "B", "W": "B", "Y": "B",
                    "A": "J", "T": "J", "H": "J", "G": "J", "P": "J", "R": "J",
                    "D": "U", "E": "U", "S": "U", "N": "U", "Q": "U", "K": "U",
                    "X":""}

    scheme_B_map = {"C": "B", "M": "B", "F": "B", "I": "B", "L": "B", "V": "B", "W": "B", "Y": "B",
                    "G": "J", "P": "J", "A": "J", "T": "J", "S": "J",
                    "E": "U", "K": "U", "R": "U", "D": "U", "N": "U", "Q": "U", "H": "U",
                    "X":""}  

    scheme_C_map = {"A": "B", "V": "B", "F": "B", "I": "B", "L": "B", "P": "B", "M": "B", "G": "B",
                    "D": "J", "E": "J", "K": "J", "R": "J",
                    "S": "U", "T": "U", "Y": "U", "C": "U", "N": "U", "Q": "U", "H": "U", "W": "U",
                    "X":""}    

    scheme_D_map = {"M": "B", "H": "B", "V": "B", "Y": "B", "N": "B", "D": "B", "I": "B",
                    "Q": "J", "L": "J", "E": "J", "K": "J", "F": "J",
                    "W": "U", "P": "U", "R": "U", "G": "U", "S": "U", "A": "U", "T": "U", "C": "U",
                    "X":""}

    scheme_E_map = {"L": "B", "A": "B", "S": "B", "G": "B", "V": "B", "T": "B", "I": "B", "P": "B", "M": "B", "C": "B",
                    "E": "J", "K": "J", "R": "J", "D": "J", "N": "J", "Q": "J", "H": "J",
                    "F": "U", "Y": "U", "W": "U",
                    "X":""}
    
    scheme = eval(f'scheme_{alpha}_map')
    return([scheme[x] for x in tess])
    

def reduce_sequence(alpha, sequence):
    scheme_A_map = {"C": "B", "M": "B", "F": "B", "I": "B", "L": "B", "V": "B", "W": "B", "Y": "B",
                    "A": "J", "T": "J", "H": "J", "G": "J", "P": "J", "R": "J",
                    "D": "U", "E": "U", "S": "U", "N": "U", "Q": "U", "K": "U",
                    "X":""}

    scheme_B_map = {"C": "B", "M": "B", "F": "B", "I": "B", "L": "B", "V": "B", "W": "B", "Y": "B",
                    "G": "J", "P": "J", "A": "J", "T": "J", "S": "J",
                    "E": "U", "K": "U", "R": "U", "D": "U", "N": "U", "Q": "U", "H": "U",
                    "X":""} 

    scheme_C_map = {"A": "B", "V": "B", "F": "B", "I": "B", "L": "B", "P": "B", "M": "B", "G": "B",
                    "D": "J", "E": "J", "K": "J", "R": "J",
                    "S": "U", "T": "U", "Y": "U", "C": "U", "N": "U", "Q": "U", "H": "U", "W": "U",
                    "X":""}   

    scheme_D_map = {"M": "B", "H": "B", "V": "B", "Y": "B", "N": "B", "D": "B", "I": "B",
                    "Q": "J", "L": "J", "E": "J", "K": "J", "F": "J",
                    "W": "U", "P": "U", "R": "U", "G": "U", "S": "U", "A": "U", "T": "U", "C": "U",
                    "X":""}

    scheme_E_map = {"L": "B", "A": "B", "S": "B", "G": "B", "V": "B", "T": "B", "I": "B", "P": "B", "M": "B", "C": "B",
                    "E": "J", "K": "J", "R": "J", "D": "J", "N": "J", "Q": "J", "H": "J",
                    "F": "U", "Y": "U", "W": "U",
                    "X":""}
    scheme = eval(f'scheme_{alpha}_map')

    reduced_sequence = ""
    for residue in sequence:
        reduced_sequence += scheme[residue]
    return(reduced_sequence)


def make_simp_list():
    letters = ['B', 'J', 'U']
    all_simp_list = []
    for i in range(len(letters)):
        for j in range(len(letters)):
            for k in range(len(letters)):
                for l in range(len(letters)):
                    all_simp_list.append(letters[i] + letters[j] + letters[k] + letters[l])
    return all_simp_list




def get_pdb_file(filename, chain_needed):
    '''
    get_pdb_file

    Input: 
    PDBID - string for the peptide id on PDB
    chain_needed - char for the chain for an ID with more than one chain

    output:
    points: List of xyz points for each alpha carbon in the desired chain

    ''' 
    p = PDBParser()
    s = p.get_structure("", filename)

    d = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
        'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
        'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
        'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

    
    points = []
    seq = []
    for chains in s:
        if chains.id == 0:
            for chain in chains:
                if chain.id == chain_needed:
                    for residue in chain:
                        for atom in residue:
                            if atom.get_name() == 'CA':
                                vector = atom.get_vector().get_array().tolist()
                                points.append(vector)
    return points
    

def calc_likelihoods(reduced_sequence, simplex_counts):
    #Generates frequencies for the entire reduced sequence
    prob_B = reduced_sequence.count('B')/len(reduced_sequence)
    prob_J = reduced_sequence.count('J')/len(reduced_sequence)
    prob_U = reduced_sequence.count('U')/len(reduced_sequence)

    number_simplexes = sum(simplex_counts.values())
    simplex_LOR = {}

    for simplex, count in simplex_counts.items():
        simplex_probability = (eval(f'prob_{simplex[0]}')*\
                               eval(f'prob_{simplex[1]}')*\
                               eval(f'prob_{simplex[2]}')*\
                               eval(f'prob_{simplex[3]}'))
        simplex_frequency = count / number_simplexes                      
        simplex_liklihood = simplex_frequency / simplex_probability
        simplex_log_odds_ratio = math.log(simplex_liklihood, 2)
        simplex_LOR[simplex] = simplex_log_odds_ratio
    return(simplex_LOR)
        

def tesselate_st(filename, fasta, chain, alpha, start, end):
    reduced_sequence = reduce_sequence(alpha, fasta)
    points = get_pdb_file(filename, chain)
    points = points[start:end]

    if points == []:
        return(pd.DataFrame([[np.nan]* len(df_columns)], columns=df_columns))
    if len(points) < 4:
        return(pd.DataFrame([[np.nan]* len(df_columns)], columns=df_columns))
    
    residue_dict = dict(zip(np.arange(len(fasta)), fasta))
    tri = DelaunayTri(points)
    reduced_simp = []
    for j, tess in enumerate(tri.vertices):
        simp = [residue_dict[x] for x in tess]
        reduced_simp.append(''.join(reduce_simplex(alpha, simp)))
    simplex_counts = Counter(reduced_simp)
    simplex_liklihoods = calc_likelihoods(reduced_sequence, simplex_counts)

    # Put the data into a dataframe with all of the possible simplexes
    all_simp_list = make_simp_list()
    simplex_dataframe = pd.DataFrame(data = simplex_liklihoods, columns = all_simp_list, index=[0])
    simplex_dataframe.fillna(0, inplace=True)
    #rename columns to differentiate alphabet
    column_alpha_names = dict(zip(all_simp_list, [f'reducedSeq_{alpha}{simp}' for simp in all_simp_list]))
    simplex_dataframe.rename(columns=column_alpha_names, inplace = True)  
    return(simplex_dataframe)
