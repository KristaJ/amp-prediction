import dash
import dash_core_components as dcc
import dash_html_components as html
import dash_bootstrap_components as dbc
import dash_daq as daq
from dash.dependencies import Input, Output, State
import pandas as pd
import numpy as np
# import assets.process_input as pi
# import assets.calc_AA as aa
# import assets.classify as classify
# import assets.gen_dssp as gd
# import assets.helpText as ht
# import assets.layout as lo
import assets.layout as lo
import flask
import base64
import os
import plotly.graph_objects as go
from dash_extensions import Download
from dash_extensions.snippets import send_data_frame
from scipy.optimize import least_squares
import fontawesome as fa


UPLOAD_DIRECTORY = "./PDB_FILES"
FONT_AWESOME = "https://use.fontawesome.com/releases/v5.7.2/css/all.css"
default_css = "https://codepen.io/chriddyp/pen/bWLwgP.css"

external_stylesheets = [dbc.themes.SPACELAB, FONT_AWESOME, default_css]
server = flask.Flask(__name__)
app = dash.Dash(__name__, 
        server=server, 
        external_stylesheets = external_stylesheets, 
        prevent_initial_callbacks=True)

if __name__ == '__main__':
    layout2 = lo.layout()
    app.layout = layout2
    app.run_server(debug=True,
    port = 8000
    )