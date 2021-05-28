import plotly.graph_objects as go
import pandas as pd
import Bio
from Bio.Seq import Seq
from Bio.SeqUtils.ProtParam import ProteinAnalysis

def window_calc_AA(peptides, alphas):    
    # peptides.reset_index(drop = True, inplace = True)
    for index, row in peptides.iterrows():
        for alpha in alphas: 
            sequence = row['Sequence']
            ana_seq = ProteinAnalysis(sequence)
            percent = ana_seq.get_amino_acids_percent()
            for key, value in percent.items():

                peptides.at[index, key] = value
    # Just to aleviate confusion we are going to change the column names
    new_keys = ['pct_'+key for key in percent.keys()]
    col_update = dict(zip(percent.keys(), new_keys))
    peptides.rename(columns = col_update, inplace = True)
    return peptides

def calc_AA_st(seq):    
    ana_seq = ProteinAnalysis(seq)
    percent = ana_seq.get_amino_acids_percent()
    aa_pcts = pd.DataFrame([percent])
    new_keys = ['pct_'+key for key in percent.keys()]
    col_update = dict(zip(percent.keys(), new_keys))
    aa_pcts.rename(columns = col_update, inplace = True)
    return aa_pcts

def plot_aa_content(aa_pcts):
    plot_data = aa_pcts.T
    plot_data.reset_index(inplace = True)
    plot_data = plot_data.rename(columns = {'index':'aa', plot_data.columns[-1]:'pct'})
    fig = go.Figure()
    data = fig.add_trace(go.Bar(
            x=plot_data['aa'],
            y=plot_data['pct'],

        ))
    fig.update_layout(
        height = 500,
        width = 700
    )
    return fig