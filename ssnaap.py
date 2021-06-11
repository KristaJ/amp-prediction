import dash
import dash_core_components as dcc
import dash_html_components as html
import dash_bootstrap_components as dbc
import dash_daq as daq
from dash.dependencies import Input, Output, State
import pandas as pd
import numpy as np
import assets.process_input as pi
import assets.calc_AA as aa
import assets.classify as classify
import assets.gen_dssp as gd
import assets.helpText as ht
import assets.layout as lo
import flask
import base64
import os
import plotly.graph_objects as go
from dash_extensions import Download
from dash_extensions.snippets import send_data_frame
from scipy.optimize import least_squares
import fontawesome as fa
import json


UPLOAD_DIRECTORY = "./PDB_FILES"
FONT_AWESOME = "https://use.fontawesome.com/releases/v5.7.2/css/all.css"
default_css = "https://codepen.io/chriddyp/pen/bWLwgP.css"

external_stylesheets = [dbc.themes.SPACELAB, FONT_AWESOME, default_css]
server = flask.Flask(__name__)
app = dash.Dash(__name__, external_stylesheets=external_stylesheets, server=server, prevent_initial_callbacks=True, suppress_callback_exceptions=True)

default_fig = go.Figure(
    data=[
        go.Bar(
            x=[], 
            y=[]
        )
    ],
    layout=go.Layout(
        xaxis={
            'showticklabels': False,
            'ticks': '',
            'showgrid': False,
            'zeroline': False
        },
        yaxis={
            'showticklabels': False,
            'ticks': '',
            'showgrid': False,
            'zeroline': False
        },
        height=10,
        width=10
    ))


app.layout = lo.layout()

def format_model_output(algo, alpha, num_features, feature_string):
    feature_map = {'aa': 'amino acid composition percentages', 
                   'seq': 'primary sequence tri-grams',
                   'simp': '3D simplexes',
                   'struct': 'secondary structure tri-grams'}
    model_string = f"The {algo} algorithm using the {alpha} alphabet reduction scheme with "
    if num_features != "":
        model_string += f'an optimized features set consisting of {num_features} features predicts'
    else:
        feature_list = []
        if 'aa' in feature_string:
            feature_list.append("amino acid composition")
        if 'seq' in feature_string:
            feature_list.append('primary sequence tri-grams')
        if 'str' in feature_string:
            feature_list.append('secondary structure tri-grams')
        if 'simp'in feature_string:
            feature_list.append('3D simplexes')
        if len(feature_list) == 1:
            feature_output = feature_list[0]
        else:
            feature_output = ""
            for feature in feature_list:
                if feature == feature_list[-1]:
                    feature_output += f' and {feature} predicts'
                else:
                    feature_output += f' {feature},'
        model_string = model_string + feature_output

    return model_string
        

def format_df(df):
    if df.shape[0]<1:
        return ""
    return [f"{col} \t : \t {df[col].values[0]:0.1%}\n" for col in df.columns]

#save the input file name
def save_file(contents, filename):
    content_type, content_string = contents.split(',')
    decoded = base64.b64decode(content_string)
    with open(os.path.join(UPLOAD_DIRECTORY, filename), "wb") as fp:
        fp.write(decoded)

#get model file
def getModel(algo, alpha, num_features, features_list):
    
    if num_features:
        features = num_features
    else:
        features = ""
        if 'aa' in features_list:
            features += 'aa'
        if 'seq' in features_list:
            features += 'seq'
        if 'str' in features_list:
            features += 'struct'
        if 'simp'in features_list:
            features += 'simp'                        
    model_path = "models/{}/alpha_{}/{}.sav".format(algo, alpha, features)
    return(model_path)

# Update peptide file field
@app.callback(Output(component_id='upload-data', component_property='children'),
                    [Input("upload-data", "filename")])
def update_filename(filename):
    return (filename)


#process peptide input
@app.callback(
    [Output(component_id='Fasta', component_property='children'),
    Output(component_id='Structure', component_property='children'),
    Output(component_id='final_Start', component_property='children'),
    Output(component_id='final_End', component_property='children'),
    Output(component_id='final_filename', component_property='children'),
    Output('final_chain', 'children'),
    Output('graphContainer', 'style'),
    Output('aa_content', 'figure')],
    [Input(component_id='PDB_submit_button', component_property='n_clicks'),
    Input(component_id='local_submit_button', component_property='n_clicks'),
    Input(component_id='fasta_submit_button', component_property='n_clicks'),
    ],
    [State(component_id="PDBID" , component_property="value"),
    State(component_id="CHAIN" , component_property="value"),
    State(component_id="Start" , component_property="value"),
    State(component_id="End" , component_property="value"),
    State("upload-data", "filename"), 
    State("upload-data", "contents"),
    State("local_CHAIN", "value"),
    State("local_Start", "value"),
    State("local_End", "value"),
    State("fasta_sequence", "value"),
    State('direct_structure', "value")] 
)
def getPeptideParameters(PDB_n_clicks, local_n_clicks, fasta_n_clicks,
                           ID, chain, start, end,  
                           filename, contents, local_chain, local_start, local_end,
                           fasta_sequence, direct_structure):
    ctx = dash.callback_context
    if ctx.triggered:
        button_id = ctx.triggered[0]['prop_id'].split('.')[0]
    # No input
    else:
        return ("", "", "", "", "", "", "none", default_fig)
    if button_id == "PDB_submit_button":
        # PDBID entry is empty
        if ID == None:
            return ('missing PDBID', "", "", "", "", "", {"display":"none"}, default_fig)
        if chain == None:
            chain = "A"
        chain = chain.upper()
        error, fasta = gd.getFastaFromPDB(ID = ID, chain = chain)
        fasta = fasta.upper()
        structure = gd.getDSSPfromPDB(ID = ID, chain = chain).upper()
        if start == None:
            start = 0
        if end == None:
            end = len(fasta)
        final_filename = f"./PDB_FILES/pdb{ID.lower()}.ent"
    if button_id == "local_submit_button":
        if filename == None:
            return ('missing file', "", "", "", "", "", {"display":"none"}, default_fig)
        if local_chain == None:
            local_chain = "A"
        local_chain = local_chain.upper()
        chain = local_chain
        save_file(contents, filename)
        local_filename = f"./PDB_FILES/{filename}"
        error, fasta = gd.getFastaFromPDB(filename=local_filename, chain = local_chain)
        fasta = fasta.upper()
        structure = gd.getDSSPfromPDB(filename=local_filename, chain = local_chain).upper()
        if local_start == None:
            start = 0
        if local_end == None:
            end = len(fasta)
        final_filename = local_filename
    #if just just a fasta sequence is provided
    if button_id == "fasta_submit_button":
        if fasta_sequence == None:
            return ('missing Sequence', "", "", "", "", "", {"display":"none"}, default_fig)
        fasta = fasta_sequence.upper()
        if direct_structure:
            structure = direct_structure.upper()
        else:
            structure = None
        start = 0
        end = len(fasta)
        final_filename = None
        chain = "A" 
    # Error checking
    if start > end:
        error_text = f"The start index can't be larger than the end index"
        structure = ""
        start = None
        end = None
        final_filename = None
        chain = None
        aa_content = ""
        return(error_text, structure, start, end, final_filename, chain, {"display":"none"}, default_fig)        
    # if error:
    #     structure = ""
    #     start = None
    #     end = None
    #     final_filename = None
    #     chain = None
    #     aa_content = ""
    #     return(fasta, structure[start:end], start, end, final_filename, chain, {"display":"none"}, default_fig)
    if fasta[start:end].find('X') > -1:
        error_text = f"This Sequence contains non-standard amino acids.\n\
        Please trim the sequence using the start and end parameters\n{fasta[start:end]}"
        structure = ""
        start = None
        end = None
        final_filename = None
        chain = None
        aa_content = ""
        return(error_text, structure, start, end, final_filename, chain, {"display":"none"}, default_fig)
    if structure and len(fasta[start:end]) != len(structure[start:end]):
        error_text = 'Sequence and structure lengths are not the same\n\
            Please double check your inputs'
        structure = ""
        start = None
        end = None
        final_filename = None
        chain = None
        aa_content = ""
        return(error_text, structure, start, end, final_filename, chain, {"display":"none"}, default_fig)
    fasta_string = fasta[start:end]
    aa_content = aa.calc_AA_st(fasta_string)
    aa_fig = aa.plot_aa_content(aa_content)
    
    if structure:
        structure_string = structure[start:end]
    else:
        structure_string = None
    return(fasta_string, structure_string, start, end, final_filename, chain, {"display":"block"}, aa_fig)
        # return(fasta[start:end], structure[start:end], start, end, final_filename, chain, {"display":"block"}, aa_fig)
  
#toggle collapse buttons
@app.callback(
    Output("alpha_collapse", "is_open"),
    [Input("alpha_button", "n_clicks")],
    [State("alpha_collapse", "is_open")],
)
def toggle_alpha_collapse(n, is_open):
    if n:
        return not is_open
    return is_open
@app.callback(
    Output("fasta_collapse", "is_open"),
    [Input("fasta_button", "n_clicks")],
    [State("fasta_collapse", "is_open")],
)
def toggle_fasta_collapse(n, is_open):
    if n:
        return not is_open
    return is_open
@app.callback(
    Output("upload_collapse", "is_open"),
    [Input("upload_button", "n_clicks")],
    [State("upload_collapse", "is_open")],
)
def toggle_upload_collapse(n, is_open):
    if n:
        return not is_open
    return is_open
@app.callback(
    Output("PDBID_collapse", "is_open"),
    [Input("PDBID_button", "n_clicks")],
    [State("PDBID_collapse", "is_open")],
)
def toggle_PDBID_collapse(n, is_open):
    if n:
        return not is_open
    return is_open
@app.callback(
    Output("feature_collapse", "is_open"),
    [Input("feature_button", "n_clicks")],
    [State("feature_collapse", "is_open")],
)
def toggle_feature_collapse(n, is_open):
    if n:
        return not is_open
    return is_open
@app.callback(
    Output("algo_collapse", "is_open"),
    [Input("algo_button", "n_clicks")],
    [State("algo_collapse", "is_open")],
)
def toggle_algo_collapse(n, is_open):
    if n:
        return not is_open
    return is_open

@app.callback(
    Output("window_collapse", "is_open"),
    [Input("window_toggle_switch", "value")],
    [State("window_collapse", "is_open")],
)
def toggle_window_slider(window_toggle, is_open):
    return window_toggle

# Format feature group selection
@app.callback(
    [Output('feature_groups', 'options'),
    Output('feature_groups', 'value')],
    [Input(component_id='PDB_submit_button', component_property='n_clicks'),
    Input(component_id='local_submit_button', component_property='n_clicks'),
    Input(component_id='fasta_submit_button', component_property='n_clicks'),
    Input('Structure', 'children'),
    Input('final_filename', 'children')]
)
def toggle_inputs(n_click1, n_click2, n_click3, structure, filename):
    if (structure == None):
        return([
            {'label': 'Amino Acid Content', 'value': 'aa'},
            {'label': 'Sequence n-grams', 'value': 'seq'},
            {'label': 'Secondary structure n-grams', 'value': 'struct', 'disabled': 'true'},
            {'label': 'Tertiary structure simplexes', 'value': 'simp', 'disabled': 'true'}
        ],
        "")
    if (filename == None):
        return([
            {'label': 'Amino Acid Content', 'value': 'aa'},
            {'label': 'Sequence n-grams', 'value': 'seq'},
            {'label': 'Secondary structure n-grams', 'value': 'struct'},
            {'label': 'Tertiary structure simplexes', 'value': 'simp', 'disabled': 'true'}
        ],
        "")
    else:
        return([
            {'label': 'Amino Acid Content', 'value': 'aa'},
            {'label': 'Sequence n-grams', 'value': 'seq'},
            {'label': 'Secondary structure n-grams', 'value': 'struct'},
            {'label': 'Tertiary structure simplexes', 'value': 'simp'}
        ],
        "")

# Format optimized feature selection
@app.callback(
    [Output('optimized_feature', 'options'),
    Output('optimized_feature', 'value')],
    [Input('feature_groups', 'value')],
    [State('Structure', 'children'),
    State('final_filename', 'children')]
)
def toggle_inputs(feature_groups, structure, filename):
    if (structure == "") or (filename == None):
        return ([{'label': '15 features', 'value': '15', 'disabled': 'true'},
                    {'label': '25 features', 'value': '25', 'disabled': 'true'},
                    {'label': '35 features', 'value': '35', 'disabled': 'true'},
                    {'label': '50 features', 'value': '50', 'disabled': 'true'}],
                     "")
    if len(feature_groups)>0:
        return ([{'label': '15 features', 'value': '15', 'disabled': 'true'},
                    {'label': '25 features', 'value': '25', 'disabled': 'true'},
                    {'label': '35 features', 'value': '35', 'disabled': 'true'},
                    {'label': '50 features', 'value': '50', 'disabled': 'true'}],
                     "")
    if len(feature_groups)==0:
        return ([{'label': '15 features', 'value': '15'},
                    {'label': '25 features', 'value': '25'},
                    {'label': '35 features', 'value': '35'},
                    {'label': '50 features', 'value': '50'}],
                    "35")


# Classify
@app.callback(
    [Output("string_output", "value"),
    Output("explaingraphContainer", 'style'),
    Output("explain_output_graph", "figure"),
    Output("saved_data", "data")],
    [Input('predict_button', 'n_clicks')],
    [State('Fasta', 'children'),
     State('Structure', 'children'),
     State("final_filename" , "children"),
     State("final_chain" , "children"),
     State("final_Start" , "children"),
     State("final_End" , "children"),
     State('window_toggle_switch', "value"),
     State('window_size', "value"),
     State('alphabet', 'value'),
     State('algorithm', 'value'),
     State('feature_groups', 'value'),
     State('optimized_feature', 'value'),
]
)       
def getModelParams(n_click, Fasta, Structure, filename, chain, start, end, window_toggle, window_size, alpha = None, algo=None, 
feature_list=None, num_features=None):
    if (algo== ""):
        return("no algorithm selected",  {"display":"none"}, default_fig,"")
    if (alpha == ""):
        return("no alphabet selected", {"display":"none"}, default_fig, "")
    if (len(feature_list) == 0 and num_features == ""):
        return('no features selected', {"display":"none"}, default_fig, "")


    feature_names = "".join(feature_list)
    model = getModel(algo, alpha, num_features, feature_names)
    features = classify.get_features(alpha, feature_list, num_features)
    if not window_toggle:
        processed_data = pi.process_input(chain, start, end, alpha, filename, Fasta, Structure)
        pred, pred_prob = classify.classify_input(processed_data, features, model)
        explain_graph = classify.explain_prediction(features, model)
        output_string = format_model_output(algo, alpha, num_features, feature_list)

        if pred == 0:
            output_string += " a low probability antimicrobial activity."
        if pred == 1:
            output_string += " a high probability of antimicrobial activity."
        json_data = processed_data.to_json(orient='index')

        return output_string, {"display":"block"}, explain_graph, json_data
    if window_toggle:
        window_start = start

        window_end = window_start + window_size
        if window_end > len(Fasta):
            return('Window larger than Peptide', {"display":"none"}, default_fig, "")
        window_preds = {}
        processed_data_all = pd.DataFrame()
        while window_end <= end:
            processed_data = pi.process_input(chain, window_start, window_end, alpha, filename, Fasta, Structure)
            pred, pred_prob = classify.classify_input(processed_data, features, model)
            window_preds[window_start] = round(pred_prob[0][1], 6)
            processed_data['start_residue'] = window_start
            window_start += 1
            window_end = window_start + window_size
            processed_data_all = processed_data_all.append(processed_data, ignore_index = True)
        fit_line = {}
        lowest_res = 100
        for i in range(0,5):
            fit_line = np.polyfit(list(window_preds.keys()), list(window_preds.values()), i, full=True)
            if fit_line[1] < lowest_res:
                lowest_res = fit_line[1]
                lowest_poly = np.poly1d(fit_line[0])
        
        fig_window = go.Figure(data=go.Scatter(x=list(window_preds.keys()), y=list(window_preds.values()), mode='markers'))
        fig_window.add_trace(go.Scatter(x=list(window_preds.keys()), y=lowest_poly(list(window_preds.keys())),
                    mode='lines',
                    name='lines'))
        fig_window.update_layout(title_text="AMP prediction Probabilities",
                          yaxis=dict(title = "AMP Probability",
                          range = [0, 1.0]),
                          xaxis=dict(title= "Start Residue"),
                          title_font_size=30,
                          showlegend=False,
                          height = 700,
                          width = 700
                          )  
        json_data_window = processed_data_all.to_json(orient='index')
        return (f"Prediction using a window size of { window_size}", {"display":"block"}, fig_window, json_data_window)

@app.callback(Output("download", "data"),
             [Input("btn", "n_clicks")],
            [State("saved_data", "data")]
)
def generate_csv(n_nlicks, processed_data):

    df = pd.read_json(processed_data, orient='index')
    return send_data_frame(df.to_csv, filename="output.csv")

if __name__ == '__main__':
    layout2 = lo.layout()
    app.layout = layout2
    app.run_server(debug=True,
    port = 8000
    )