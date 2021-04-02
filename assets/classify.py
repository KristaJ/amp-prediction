import pickle
import pandas as pd
import eli5
from sklearn.linear_model import LogisticRegression
import warnings
warnings.filterwarnings('ignore')
# from shapash.explainer.smart_explainer import SmartExplainer

import plotly.graph_objects as go


def classify_input(processed_data, features, model):
    model_file = f"./assets/{model}"
    clf = pickle.load(open(model_file, 'rb'))
    X_sample = processed_data[features]
    y_pred = clf.predict(X_sample)
    y_prob = clf.predict_proba(X_sample)
    return y_pred[0], y_prob

def get_features(alpha, feature_list= None, num_features=None):
    simp_noalpha = ['reducedSeq_{}BBBB',
                    'reducedSeq_{}BBBJ',
                    'reducedSeq_{}BBBU',
                    'reducedSeq_{}BBJB',
                    'reducedSeq_{}BBJJ',
                    'reducedSeq_{}BBJU',
                    'reducedSeq_{}BBUB',
                    'reducedSeq_{}BBUJ',
                    'reducedSeq_{}BBUU',
                    'reducedSeq_{}BJBB',
                    'reducedSeq_{}BJBJ',
                    'reducedSeq_{}BJBU',
                    'reducedSeq_{}BJJB',
                    'reducedSeq_{}BJJJ',
                    'reducedSeq_{}BJJU',
                    'reducedSeq_{}BJUB',
                    'reducedSeq_{}BJUJ',
                    'reducedSeq_{}BJUU',
                    'reducedSeq_{}BUBB',
                    'reducedSeq_{}BUBJ',
                    'reducedSeq_{}BUBU',
                    'reducedSeq_{}BUJB',
                    'reducedSeq_{}BUJJ',
                    'reducedSeq_{}BUJU',
                    'reducedSeq_{}BUUB',
                    'reducedSeq_{}BUUJ',
                    'reducedSeq_{}BUUU',
                    'reducedSeq_{}JBBB',
                    'reducedSeq_{}JBBJ',
                    'reducedSeq_{}JBBU',
                    'reducedSeq_{}JBJB',
                    'reducedSeq_{}JBJJ',
                    'reducedSeq_{}JBJU',
                    'reducedSeq_{}JBUB',
                    'reducedSeq_{}JBUJ',
                    'reducedSeq_{}JBUU',
                    'reducedSeq_{}JJBB',
                    'reducedSeq_{}JJBJ',
                    'reducedSeq_{}JJBU',
                    'reducedSeq_{}JJJB',
                    'reducedSeq_{}JJJJ',
                    'reducedSeq_{}JJJU',
                    'reducedSeq_{}JJUB',
                    'reducedSeq_{}JJUJ',
                    'reducedSeq_{}JJUU',
                    'reducedSeq_{}JUBB',
                    'reducedSeq_{}JUBJ',
                    'reducedSeq_{}JUBU',
                    'reducedSeq_{}JUJB',
                    'reducedSeq_{}JUJJ',
                    'reducedSeq_{}JUJU',
                    'reducedSeq_{}JUUB',
                    'reducedSeq_{}JUUJ',
                    'reducedSeq_{}JUUU',
                    'reducedSeq_{}UBBB',
                    'reducedSeq_{}UBBJ',
                    'reducedSeq_{}UBBU',
                    'reducedSeq_{}UBJB',
                    'reducedSeq_{}UBJJ',
                    'reducedSeq_{}UBJU',
                    'reducedSeq_{}UBUB',
                    'reducedSeq_{}UBUJ',
                    'reducedSeq_{}UBUU',
                    'reducedSeq_{}UJBB',
                    'reducedSeq_{}UJBJ',
                    'reducedSeq_{}UJBU',
                    'reducedSeq_{}UJJB',
                    'reducedSeq_{}UJJJ',
                    'reducedSeq_{}UJJU',
                    'reducedSeq_{}UJUB',
                    'reducedSeq_{}UJUJ',
                    'reducedSeq_{}UJUU',
                    'reducedSeq_{}UUBB',
                    'reducedSeq_{}UUBJ',
                    'reducedSeq_{}UUBU',
                    'reducedSeq_{}UUJB',
                    'reducedSeq_{}UUJJ',
                    'reducedSeq_{}UUJU',
                    'reducedSeq_{}UUUB',
                    'reducedSeq_{}UUUJ',
                    'reducedSeq_{}UUUU']
    seq_noalpha = ['{}_BBB',
                  '{}_BBJ',
                  '{}_BBU',
                  '{}_BJB',
                  '{}_BJJ',
                  '{}_BJU',
                  '{}_BUB',
                  '{}_BUJ',
                  '{}_BUU',
                  '{}_JBB',
                  '{}_JBJ',
                  '{}_JBU',
                  '{}_JJB',
                  '{}_JJJ',
                  '{}_JJU',
                  '{}_JUB',
                  '{}_JUJ',
                  '{}_JUU',
                  '{}_UBB',
                  '{}_UBJ',
                  '{}_UBU',
                  '{}_UJB',
                  '{}_UJJ',
                  '{}_UJU',
                  '{}_UUB',
                  '{}_UUJ',
                  '{}_UUU']
    simp = [x.format(alpha) for x in simp_noalpha]
    seq = [x.format(alpha) for x in seq_noalpha]
    struct = ['st_CCC',
                      'st_CCH',
                      'st_CCB',
                      'st_CHC',
                      'st_CHH',
                      'st_CHB',
                      'st_CBC',
                      'st_CBH',
                      'st_CBB',
                      'st_HCC',
                      'st_HCH',
                      'st_HCB',
                      'st_HHC',
                      'st_HHH',
                      'st_HHB',
                      'st_HBC',
                      'st_HBH',
                      'st_HBB',
                      'st_BCC',
                      'st_BCH',
                      'st_BCB',
                      'st_BHC',
                      'st_BHH',
                      'st_BHB',
                      'st_BBC',
                      'st_BBH',
                      'st_BBB']
    aa = ['pct_A',
                          'pct_C',
                          'pct_D',
                          'pct_E',
                          'pct_F',
                          'pct_G',
                          'pct_H',
                          'pct_I',
                          'pct_K',
                          'pct_L',
                          'pct_M',
                          'pct_N',
                          'pct_P',
                          'pct_Q',
                          'pct_R',
                          'pct_S',
                          'pct_T',
                          'pct_V',
                          'pct_W',
                          'pct_Y']
    
    if num_features:
        opt_features = pd.read_csv('./assets/models/features.csv', index_col = None)
        features = opt_features[(opt_features.alpha == alpha) & (opt_features.num_features == int(num_features))]['features_list'].\
            str.\
            slice(2, -2).\
            str.split("', '").\
            tolist()[0]
    else:
        features = []
        for feature in feature_list:
            features += eval(feature)
    return features


def explain_prediction(features, model):
    model_file = f"./assets/{model}"
    clf = pickle.load(open(model_file, 'rb'))

    temp = eli5.format_as_dict(eli5.explain_weights(clf, feature_names=features))
    features = [x['feature'] for x in temp['feature_importances']['importances']][::-1]
    weights = [x['weight'] for x in temp['feature_importances']['importances']][::-1]
    
    explain_fig = go.Figure(go.Bar(
            x=weights,
            y=features,
            orientation='h'))

    explain_fig.update_layout(
        title="Features with highest prediction weights",
        xaxis_title="Prediction Weight",
        yaxis = dict(
            title = "Feature",
            dtick =  1),
        height = 300,
        width = 500
    )
    return explain_fig
    
    
    