import dash_core_components as dcc
import dash_html_components as html
import dash_bootstrap_components as dbc

def intro_text():
    helpText='''
        nGAMPP makes use of sequence and structure peptide information to predict the probability of 
        antimicrobial activity
        ###### ALPHABET REDUCTION
        The encoding process first utilizes the selected alphabet reduction scheme to reduce the amino acid alphabet from 
        20 residues to 3 representative residues. This reduces sparsity in the resulting dataset.  
        ###### NGRAMS
        Sequences with reduced alphabets are then processed into N-grams.  Trigrams are used for sequence and 
        secondary structure information while four-grams are used for tirtiary structures.  N-grams allow for data encoding 
        while maintaining the integrity of the sequence.
        ###### LOG-ODDS-RATIO
        N-gram sequences are transformed into log-odds-ratios using the following equation:
        '''
    helpText2 = '''
        ###### PREDICTION
        Once the peptide has been encoded one of 3 trained prediction algorithms can be run using either complete feature
        sets or a selection of optimized feature sets
        '''
    return helpText, helpText2
