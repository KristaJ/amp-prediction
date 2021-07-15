import dash_core_components as dcc
import dash_html_components as html
import dash_bootstrap_components as dbc
import dash_daq as daq
import plotly.graph_objects as go
from dash_extensions import Download
from . import helpText

introtext1, introtext2 = helpText.intro_text()

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
# def make_tip_button(button_id):
#         button = html.I(
#                 className="fas fa-question-circle",
#                 id = button_id
#             )
#   return button


def make_tip_button(button_id):
        button = html.Button(
                "\U00002754",
                id = button_id,
                className = "tip_button"
            )
        return button
def make_tip(button_id, tooltip_text):
        tip = dbc.Tooltip(
            tooltip_text,
            target=button_id,
            placement="right",
        )
        return tip
def make_alpha_tip(button_id, tooltip_text):
        tip = dbc.Tooltip(
            tooltip_text,
            target=button_id,
            placement="top-start",
        )
        return tip       
def peptide_input_card():
    PDBInputCard = dbc.Card([
        dbc.CardBody([
            dbc.Input(id="PDBID", placeholder="PDBID", type="text"),
            dbc.Input(id="CHAIN", placeholder="chain", type="text"),
            dbc.Input(id="Start", placeholder="start", type="number"),
            dbc.Input(id="End", placeholder="end", type="number"),
            dbc.Button("submit", id = 'PDB_submit_button', outline=True, color="primary", size = 'sm'),
            ])
    ])
    localInputCard = dbc.Card([
        dbc.CardBody([
            dcc.Upload(
                id='upload-data',
                children= html.A('Select Files'),
                style={
                    'width': '100%',
                    'height': '30px',
                    'lineHeight': '60px',
                    'borderWidth': '1px',
                    'borderStyle': 'dashed',
                    'borderRadius': '5px',
                    'textAlign': 'center',
                    'margin': '2px'
                },
            ),
            dbc.Input(id="local_CHAIN", placeholder="chain", type="text"),
            dbc.Input(id="local_Start", placeholder="start", type="number"),
            dbc.Input(id="local_End", placeholder="end", type="number"),
            dbc.Button("submit", id = 'local_submit_button', outline=True, color="primary", size = 'sm'),
            ])
    ])
    FastaInputCard = dbc.Card([
        dbc.CardHeader("3D simplex features will not be available"),
        dbc.CardHeader("If secondary structure is not included structure features will not be available"),
        dbc.CardBody([
            dbc.Input(id="fasta_sequence", placeholder="Primary Structure", type="text"),
            dbc.Input(id="direct_structure", placeholder="secondary structure (optional)", type="text"),
            dbc.Button("submit", id = 'fasta_submit_button', outline=True, color="primary", size = 'sm'),
            ])
    ])
    input_card = dbc.Card([
        dbc.CardBody([
            dbc.Button(
                "Input Sequence and Structure using a PDB ID",
                id="PDBID_button",
                color="primary"
            ),
            make_tip_button("PDBID_tip_button"),
            make_tip("PDBID_tip_button", "Uses an ID from the PDB archive"),
            dbc.Collapse(
                id='PDBID_collapse',
                is_open = False,
                children = PDBInputCard
            )
        ]),
        html.Br(),
        dbc.CardBody([
            dbc.Button(
                "Upload PBD file",
                id="upload_button",
                color="primary"
            ),
            make_tip_button("upload_tip_button"),
            make_tip("upload_tip_button", "Upload a custom PDB formatted file"),
            dbc.Collapse(
                id='upload_collapse',
                is_open = False,
                children = localInputCard
            )
        ]),
        html.Br(),
        dbc.CardBody([
            dbc.Button(
                "input primary and secondary structures",
                id="fasta_button",
                color="primary"
            ),
            make_tip_button("fasta_tip_button"),
            make_tip("fasta_tip_button", 
                    [html.P("Paste a custom fasta sequence with an optional structure sequence"),
                    html.P("Simplex features will not be available"),
                    html.P("Structure features will not be available if no structure is provided")
                    ]),
            dbc.Collapse(
                id='fasta_collapse',
                is_open = False,
                children = FastaInputCard
            ),
        ]),
    ])
    return input_card
def algo_input_card():

    algo_select_card = dbc.Card(
        dbc.CardBody(
            [
                dbc.Alert('Classification Algorithm', color="info"),
                dbc.RadioItems(
                    id = 'algorithm',
                    options=[
                        {'label': 'ExtraTrees', 'value': 'ExtraTrees'},
                        {'label': 'Gradient Boosted Clasifier', 'value': 'GBC'},
                        {'label': 'AdaBoost', 'value': 'AdaBoost'},
                        
                    ],
                    value='ExtraTrees'
                )  
            ]
        )
    )
    Alpha_selection = dbc.RadioItems(
        id = "alphabet",
        options=[
            {'label': "Minimized Mismatch \U00002754", 'value': 'A', 'label_id': "min_mismatch"},
            {'label': 'Monte Carlo Reduction \U00002754', 'value': 'B', 'label_id': "MC_reduction"},
            {'label': 'Chemical Property Similarity \U00002754', 'value': 'C', 'label_id': "chem"},
            {'label': 'Molecular Recognition Theory \U00002754', 'value': 'D', 'label_id': "MRT"},
            {'label': 'Deterministic Reduction \U00002754', 'value': 'E', 'label_id': "Determ"},
        ],
        value='A'
    )
    alphabet_card = dbc.Card(
        dbc.CardBody(
            [
                dbc.Alert('Alphabet Reduction Scheme', color="info"),
                Alpha_selection,
                make_alpha_tip("min_mismatch", [html.Em("Minimized Mismatch"),
                                        html.Br(),
                                        "B maps to CMFILVWY",
                                        html.Br(),
                                        "J maps to ATHGPR",
                                        html.Br(),
                                        "U maps to DESNQK"]),
                make_alpha_tip("MC_reduction", [html.Em("Monte Carlo Reduction"),
                                        html.Br(),
                                        "B maps to CMFILVWY",
                                        html.Br(),
                                        "J maps to GPATS",
                                        html.Br(),
                                        "U maps to EKRDNQH"]),
                make_alpha_tip("chem", [html.Em("Chemical Property Similarity"),
                                        html.Br(),
                                        "B maps to AVFILPMG",
                                        html.Br(),
                                        "J maps to DEKR",
                                        html.Br(),
                                        "U maps to STYCNQHW"]),
                make_alpha_tip("MRT", [html.Em("Molecular Recognition Theory"),
                                        html.Br(),
                                        "B maps to MHVYNDI",
                                        html.Br(),
                                        "J maps to QLEKF",
                                        html.Br(),
                                        "U maps to WPRGSATC"]),
                make_alpha_tip("Determ", [html.Em("Deterministic Reduction"),
                                        html.Br(),
                                        "B maps to LASGVTIPMC",
                                        html.Br(),
                                        "J maps to EKRDNQH",
                                        html.Br(),
                                        "U maps to FYW"]),
            ]
        )
    )
    feature_group_card = dbc.Card(
        dbc.CardBody(
            [
                dbc.Alert(["Select Feature Groups  ", make_tip_button("group_tip_button")], color="info"),
                make_tip("group_tip_button", 
                        ["Select one or more features groups",
                        html.Br(),
                        "Feature groups can't be combined with optimized feature sets"
                        ]
                ),
                dbc.Checklist(
                    id = 'feature_groups',
                    options=[
                        {'label': 'Amino Acid Content', 'value': 'aa'},
                        {'label': 'Sequence n-grams', 'value': 'seq'},
                        {'label': 'Secondary structure n-grams', 'value': 'struct'},
                        {'label': 'Tertiary structure simplexes', 'value': 'simp'}
                    ],
                    value=[]
                )  
            ]
        )
    )
    optimized_feature_card = dbc.Card(
        dbc.CardBody(
            [
                dbc.Alert(["Select optimized Feature Set  ", make_tip_button("optimized_tip_button")], color="info"),
                make_tip("optimized_tip_button", 
                        ["Select from a predefined set of optimized features",
                        html.Br(),
                        "Features sets are optimized by alphabet reduction scheme and algorithm"
                        ]
                ),
                dbc.RadioItems(
                    id = 'optimized_feature',
                    options=[
                        {'label': '15 features', 'value': '15'},
                        {'label': '25 features', 'value': '25'},
                        {'label': '35 features', 'value': '35'},
                        {'label': '50 features', 'value': '50'}
                    ],
                    value='25',
                )  
            ]
        )
    )
    window_card = dbc.Card(
        dbc.CardBody(
            [ 
                daq.ToggleSwitch(
                    id='window_toggle_switch',
                    value=False,
                    label='Use Prediction Window',
                    labelPosition='top',
                    color = 'blue',
                    size = 30
                )
            ]
        ),
        outline=False
    )
    window_slider = dcc.Slider(
                id='window_size',
                min=5,
                max=25,
                step=1,
                value=20, 
                tooltip = {'always_visible': False,
                            'placement': "top"}
            )
    
    algo_card = dbc.Card([
        dbc.CardBody([
            dbc.Button(
                "Select Algorithm",
                id="algo_button",
                color="primary"
            ),
            make_tip_button("algo_tip_button"),
            make_tip("algo_tip_button", "Select the Algorithm to use for prediction"),
            dbc.Collapse(
                id='algo_collapse',
                is_open = False,
                children = algo_select_card
            )
        ]),
        html.Br(),
        dbc.CardBody([
            dbc.Button(
                "Select Alphabet Reduction Scheme",
                id="alpha_button",
                color="primary"
            ),
            dbc.Collapse(
                id='alpha_collapse',
                is_open = False,
                children = alphabet_card
            )
        ]),

        html.Br(),
        dbc.CardBody([
            dbc.Button(
                "Select Features",
                id="feature_button",
                color="primary"
            ),
            dbc.Collapse(
                id='feature_collapse',
                is_open = False,
                children = [optimized_feature_card,
                feature_group_card]
            )
        ]),
        html.Br(),
        dbc.CardBody([
            window_card,
            dbc.Collapse(
                id='window_collapse',
                is_open = False,
                children = [
                    "Select Window Size",
                    window_slider]
            )
        ]),
        dbc.Button(
            "PREDICT",
            id="predict_button",
            className="mb-3",
            color="primary"
        )
    ])
    return algo_card
def output_card():
    peptideOutputCard = dbc.Card([
        dbc.CardBody([
            html.H6("Peptide Characteristics"),
            html.Br(),
            html.U("Primary Structure"),
            html.Div(id='Fasta', style={'whiteSpace': 'pre-line'}),
            html.U("Secondary Structure"),
            html.Div(id='Structure', style={'whiteSpace': 'pre-line'}),
            html.U("Amino Acid Content"),
            html.Div(
                id="graphContainer",
                children = [
                    dcc.Graph(id="aa_content", figure=default_fig),
                ],
                style={"display":"none"}
            ),
            html.H5(id = 'final_Start', children = "", style={"display":"none"}),
            html.H5(id = 'final_End', style={"display":"none"}),
            html.H5(id = 'final_filename', style={"display":"none"}),
            html.H5(id = 'final_chain', style={"display":"none"}),
        ]),
    ])
    prediction_card= dbc.Container([
        html.Br(),
        dbc.Spinner(
            dcc.Textarea(id = 'string_output',
                        value = "",
                        style = {'width': '90%', 'height': 100}
            )
        ),   
        dbc.Spinner([
            html.Div(
                id="explaingraphContainer",
                children = [
                    dcc.Graph(id="explain_output_graph", figure=default_fig),
                ],
                style={"display":"none"}
            )
        ]), 
        dcc.Store(id='saved_data'),
        html.Div([html.Button("Download csv", id="btn"), Download(id="download")])   
    ])
    output_card = html.Div([
        peptideOutputCard,
        prediction_card
    ])
    return output_card


def layout():
    layout = html.Div(
        id = "app-body",
        className = "app-body",
        children = [
            html.Div(
                id="center",
                children=[
                    html.Img(src="assets/ssnaap_logo.png",
                             height = 400)
                ],
            ),
            html.Div(
                id = "input-tabs",
                className = "input-tabs",
                children = [
                    dbc.Tabs([
                        dbc.Tab(
                            label='Introduction', 
                            children=[
                                dcc.Markdown(introtext1),
                                html.Img(src="assets/equation.png",
                                        width="80%"),
                                dcc.Markdown(introtext2),
                                html.A(href="assets/supplemental.html",
                                children="Supplemental Informationßß"),
                            ]
                        ),
                        dbc.Tab(
                            label='Peptide Information', 
                            children=[peptide_input_card()]),
                        dbc.Tab(
                            label='Algorithm Information', 
                            children=[algo_input_card()]),
                    ])
                ]
            ),
            html.Div(
                id='output-card',
                className="output-card",
                children=[
                    output_card()
                ]
            )
        ]
    )

    return layout        

