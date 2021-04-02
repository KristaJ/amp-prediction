import pandas as pd
import math
from os import listdir
from os.path import isfile, join
import numpy as np
from collections import Counter

##############################################################
###   This code will generate Ngram liklihoods for each    ###
###   reduced alphabet                                     ###
##############################################################
def gen_all_seq_ngrams():
    letters = ['B', 'U', 'J']
    n_gram_names = []
    for letter1 in letters:
        for letter2 in letters:
            for letter3 in letters:
                n_gram_names.append(letter1+letter2+letter3)
    return n_gram_names


def gen_all_struct_ngrams():
    letters = ['B', 'H', 'C']
    n_gram_names = []
    for letter1 in letters:
        for letter2 in letters:
            for letter3 in letters:
                n_gram_names.append(letter1+letter2+letter3)
    return n_gram_names


def Generate_seq_ngrams(sequence, alpha):
    #Generates frequencies for the entire reduced sequence
    prob_B = sequence.count('B')/len(sequence)
    prob_J = sequence.count('J')/len(sequence)
    prob_U = sequence.count('U')/len(sequence)

    n_grams = [sequence[x : x+3] for x in range(len(sequence)-3)]
    number_n_grams = len(n_grams)
    n_gram_counts = Counter(n_grams)
    
    n_gram_LOR = {}
    for n_gram, count in n_gram_counts.items():
        n_gram_probability = (eval(f'prob_{n_gram[0]}')*\
                            eval(f'prob_{n_gram[1]}')*\
                            eval(f'prob_{n_gram[2]}'))
        n_gram_frequency = count / number_n_grams                      
        n_gram_liklihood = n_gram_frequency / n_gram_probability
        n_gram_log_odds_ratio = math.log(n_gram_liklihood, 2)
        n_gram_LOR[n_gram] = n_gram_log_odds_ratio
    
    all_seq_n_grams = gen_all_seq_ngrams()
    seq_n_gram_dataframe = pd.DataFrame(data = n_gram_LOR, columns = all_seq_n_grams, index=[0])
    seq_n_gram_dataframe.fillna(0, inplace=True)
    #rename columns to differentiate alphabet
    column_alpha_names = dict(zip(all_seq_n_grams, [f'{alpha}_{n_gram}' for n_gram in all_seq_n_grams]))
    seq_n_gram_dataframe.rename(columns=column_alpha_names, inplace = True)  
    return(seq_n_gram_dataframe)


def Generate_struct_ngrams(structure):
    #Generates frequencies for the entire reduced structure
    prob_B = structure.count('B')/len(structure)
    prob_H = structure.count('H')/len(structure)
    prob_C = structure.count('C')/len(structure)

    n_grams = [structure[x : x+3] for x in range(len(structure)-3)]
    number_n_grams = len(n_grams)
    n_gram_counts = Counter(n_grams)
    
    n_gram_LOR = {}
    for n_gram, count in n_gram_counts.items():
        n_gram_probability = (eval(f'prob_{n_gram[0]}')*\
                            eval(f'prob_{n_gram[1]}')*\
                            eval(f'prob_{n_gram[2]}'))
        n_gram_frequency = count / number_n_grams                      
        n_gram_liklihood = n_gram_frequency / n_gram_probability
        n_gram_log_odds_ratio = math.log(n_gram_liklihood, 2)
        n_gram_LOR[n_gram] = n_gram_log_odds_ratio
    
    all_struct_n_grams = gen_all_struct_ngrams()
    struct_n_gram_dataframe = pd.DataFrame(data = n_gram_LOR, columns = all_struct_n_grams, index=[0])
    struct_n_gram_dataframe.fillna(0, inplace=True)
    #rename columns to differentiate as structure
    column_alpha_names = dict(zip(all_struct_n_grams, [f'st_{n_gram}' for n_gram in all_struct_n_grams]))
    struct_n_gram_dataframe.rename(columns=column_alpha_names, inplace = True)  
    return(struct_n_gram_dataframe)



