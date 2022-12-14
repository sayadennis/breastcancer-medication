import re
import numpy as np
import pandas as pd
from scipy.stats import chisquare
from sklearn.utils import resample
from sklearn.linear_model import LogisticRegression

## Load data 
din = '/share/fsmresfiles/breast_cancer_medication/data/02_statistics_input'
dout = '/share/fsmresfiles/breast_cancer_medication/data/03_statistics_output'

data = pd.read_csv(f'{din}/labeled_cohort_data.csv')

cts = {}
for med in ['metformin', 'hydrostatin', 'lipostatin']:
    cts[med] = {}
    for outcome in ['biocat', 'recur']:
        cts[med][outcome] = pd.read_csv(f'{din}/{med}_{outcome}.csv', index_col=0)

##########################
#### Define functions ####
##########################

def get_oddsratio_ci(X, y, alpha=0.95, rep=5000):
    oddsratio = {}
    for colnum in range(X.shape[1]):
        oddsratio[colnum] = []
    # 
    for i in range(rep):
        X_bs, y_bs = resample(X, y, random_state=i) # create bootstrap (bs) sample
        if ~np.all(y_bs==0):
            lrm = LogisticRegression(penalty='l2', solver='lbfgs')
            lrm.fit(X_bs, y_bs)
            for colnum in range(X.shape[1]):
                oddsratio[colnum].append(np.exp(lrm.coef_[0][colnum]))
        else:
            continue
    mean_or = [np.mean(oddsratio[colnum]) for colnum in range(X.shape[1])]
    ci = []
    for colnum in range(X.shape[1]):
        # first get ci1
        p = ((1.0-alpha)/2.0) * 100
        lower = max(0.0, np.percentile(oddsratio[colnum], p))
        p = (alpha+((1.0-alpha)/2.0)) * 100
        upper = np.percentile(oddsratio[colnum], p)
        ci.append((lower, upper))
    return mean_or, ci


##########################################################
#### Stratified analysis 1 - three biomarker subtypes ####
##########################################################

#### Odds-ratio and confidence intervals calculated via logistic regression ####
meds = ['metformin', 'lipostatin', 'hydrostatin']
bm_cats = ['ER/PR+ HER2-', 'Triple Negative', 'HER2+']

record_ors = pd.DataFrame(None, index=meds, columns=bm_cats)

for bm_cat in bm_cats:
    sub_data = data.iloc[data.biomarker_subtypes.values==bm_cat,:]
    for med in meds:
        X = sub_data[f'using_{med}'].to_numpy(dtype=float).reshape(-1,1)
        y = sub_data['recurrence'].map(
            {'None' : 0, 'Both' : 1, 'Local' : 1, 'Distant' : 1}
        ).fillna(0).to_numpy()
        # perform LR analysis 
        or_mean, ci = get_oddsratio_ci(X, y)
        record_ors.loc[med, bm_cat] = f'{or_mean[0]:.2f} ({ci[0][0]:.2f}-{ci[0][1]:.2f})'

record_ors.to_csv(f'{dout}/logistic_oddsratio_biomarkercategory.csv', index=True)

#########################################################
#### Stratified analysis 2 - ER/PR positive vs. rest ####
#########################################################

#### Odds-ratio and confidence intervals calculated via logistic regression ####
erpr_cats = ['ER/PR+', 'ER/PR-']
record_ors = pd.DataFrame(None, index=meds, columns=erpr_cats)

for erpr_cat in erpr_cats:
    sub_data = data.iloc[data.er_pr.values==erpr_cat,:]
    for med in meds:
        X = sub_data[f'using_{med}'].to_numpy(dtype=float).reshape(-1,1)
        y = sub_data['recurrence'].map(
            {'None' : 0, 'Both' : 1, 'Local' : 1, 'Distant' : 1}
        ).fillna(0).to_numpy()
        # perform LR analysis 
        or_mean, ci = get_oddsratio_ci(X, y)
        record_ors.loc[med, erpr_cat] = f'{or_mean[0]:.2f} ({ci[0][0]:.2f}-{ci[0][1]:.2f})'

record_ors.to_csv(f'{dout}/logistic_oddsratio_erpr_category.csv', index=True)

## Do the same, but further stratify for menopausal status 
for meno_status in ['pre', 'post']:
    erpr_cats = ['ER/PR+', 'ER/PR-']
    record_ors = pd.DataFrame(None, index=meds, columns=erpr_cats)
    for erpr_cat in erpr_cats:
        sub_data = data.iloc[((data.er_pr.values==erpr_cat) & (data.menopause_status.values==meno_status)),:]
        for med in meds:
            X = sub_data[f'using_{med}'].to_numpy(dtype=float).reshape(-1,1)
            y = sub_data['recurrence'].map(
                {'None' : 0, 'Both' : 1, 'Local' : 1, 'Distant' : 1}
            ).fillna(0).to_numpy()
            # perform LR analysis 
            or_mean, ci = get_oddsratio_ci(X, y)
            record_ors.loc[med, erpr_cat] = f'{or_mean[0]:.2f} ({ci[0][0]:.2f}-{ci[0][1]:.2f})'
    record_ors.to_csv(f'{dout}/logistic_oddsratio_erpr_category_{meno_status}menopausal.csv', index=True)
