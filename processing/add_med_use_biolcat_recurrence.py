import re
import numpy as np
import pandas as pd

## Load data 
din = '/share/fsmresfiles/breast_cancer_medication/data/01_ssms'
dout = '/share/fsmresfiles/breast_cancer_medication/data/02_statistics_input'

# just the cohort's labeled data 
data = pd.read_csv(f'{din}/6800_cohort.csv', header=None)

with open(f'{din}/6800_cohort_colnames.txt', 'r') as f:
    lines = f.readlines()

colnames = [line.strip() for line in lines]
colnames.remove('diagnosis')
for i in range(len(colnames)):
    if colnames[i]=='total_mamography_before':
        colnames[i] = 'total_mamography_before_diagnosis'

data.columns = colnames

# drug use data
metformin = pd.read_csv(f'{din}/test_6800_metformin_use.csv', header=None)
statins = pd.read_csv(f'{din}/test_6800_statins_use.csv', header=None)

with open(f'{din}/test_6800_metformin_use_colnames.txt', 'r') as f:
    lines = f.readlines()

# colnames = [line.strip() for line in lines]
colnames = lines[0].split()
colnames.remove('diagnosis')
for i in range(len(colnames)):
    if colnames[i]=='total_mamography_before':
        colnames[i] = 'total_mamography_before_diagnosis'

metformin.columns = colnames
statins.columns = colnames

metformin.order_start_date_key = pd.to_datetime(metformin.order_start_date_key)
metformin.order_end_date_key = pd.to_datetime(metformin.order_end_date_key)
metformin.cancer_directed_surgery_date = pd.to_datetime(metformin.cancer_directed_surgery_date)

statins.order_start_date_key = pd.to_datetime(statins.order_start_date_key)
statins.order_end_date_key = pd.to_datetime(statins.order_end_date_key)
statins.cancer_directed_surgery_date = pd.to_datetime(statins.cancer_directed_surgery_date)

def drug_during_surg(df_row):
    surg_date = df_row['cancer_directed_surgery_date']
    start = df_row['order_start_date_key']
    end = df_row['order_end_date_key']
    if (pd.isnull(start) | pd.isnull(end)):
        return False
    elif ((start <= surg_date) & pd.isnull(end)):
        return True
    elif ((start <= surg_date) & (end > surg_date)):
        return True
    else:
        return False

## Add metformin use 
metformin['using'] = metformin.apply(drug_during_surg, axis=1)
metformin_ir_id = metformin.iloc[metformin.using.values,:].patient_ir_id.unique()

data['using_metformin'] = False
metformin_ix = data.iloc[[ir_id in metformin_ir_id for ir_id in data.patient_ir_id.values],:].index
data.loc[metformin_ix,'using_metformin'] = True

## Add statin use - separated by lipophillic and hydrophillic
statin_types = {
    'atorvastatin' : 'lipophilic',
    'simvastatin' : 'lipophilic',
    'lovastatin' : 'lipophilic',
    'fluvastatin' : 'lipophilic',
    'cerivastatin' : 'lipophilic',
    'pitavastatin' : 'lipophilic',
    'rosuvastatin' : 'hydrophilic',
    'pravastatin' : 'hydrophilic',
}

def get_statin_type(df_row):
    generic_name = df_row['generic_name']
    if pd.isnull(generic_name):
        return None
    else:
        patterns = [re.findall(x, generic_name.lower(), re.IGNORECASE) for x in statin_types.keys()]
        patterns = [item for sublist in patterns for item in sublist]
        if len(patterns)==0:
            return 'unknown'
        elif len(patterns)>1:
            return 'found multiple'
        else: # only one match 
            return statin_types[patterns[0]]

statins['type'] = statins.apply(get_statin_type, axis=1)

statins['using'] = statins.apply(drug_during_surg, axis=1)
lipo_statins_ir_id = statins.iloc[statins.using.values & (statins['type'].values=='lipophilic'),:].patient_ir_id.unique()
hydro_statins_ir_id = statins.iloc[statins.using.values & (statins['type'].values=='hydrophilic'),:].patient_ir_id.unique()

data['using_lipostatin'] = False; data['using_hydrostatin'] = False;
lipo_statin_ix = data.iloc[[ir_id in lipo_statins_ir_id for ir_id in data.patient_ir_id.values],:].index
hydro_statin_ix = data.iloc[[ir_id in hydro_statins_ir_id for ir_id in data.patient_ir_id.values],:].index
data.loc[lipo_statin_ix,'using_lipostatin'] = True
data.loc[hydro_statin_ix,'using_hydrostatin'] = True

recur = []
for i in data.index:
    if ((data.loc[i,'Local_Recurrence']=='YES') & (data.loc[i,'Distant_Recurrence']=='YES')):
        recur.append('Both')
    elif data.loc[i,'Local_Recurrence']=='YES':
        recur.append('Local')
    elif data.loc[i,'Distant_Recurrence']=='YES':
        recur.append('Distant')
    else:
        recur.append('None')

bm_cat = []
for i in data.index:
    if ((data.loc[i,'ER']=='POSITIVE') & (data.loc[i,'PR']=='POSITIVE') & (data.loc[i,'HER2']=='NEGATIVE')):
        bm_cat.append('ER/PR+ HER2-')
    elif data.loc[i,'HER2']=='POSITIVE':
        bm_cat.append('HER2+')
    else:
        bm_cat.append('Triple Negative')

erpr_cat = []
for i in data.index:
    if ((data.loc[i,'ER']=='POSITIVE') | (data.loc[i,'PR']=='POSITIVE')):
        erpr_cat.append('ER/PR+')
    else:
        erpr_cat.append('ER/PR-')


## Summary count tables
data['recurrence'] = recur
data['biomarker_subtypes'] = bm_cat
data['er_pr'] = erpr_cat

data.to_csv(f'{dout}/labeled_cohort_data.csv')

## counts for chi-squared tests
crosstab = pd.crosstab(data.using_metformin, data.biomarker_subtypes)
crosstab.columns.name = None
crosstab.to_csv(f'{dout}/metformin_biocat.csv')

crosstab = pd.crosstab(data.using_metformin, data.recurrence)
crosstab.columns.name = None
crosstab.to_csv(f'{dout}/metformin_recur.csv')

crosstab = pd.crosstab(data.using_lipostatin, data.biomarker_subtypes)
crosstab.columns.name = None
crosstab.to_csv(f'{dout}/lipostatin_biocat.csv')

crosstab = pd.crosstab(data.using_lipostatin, data.recurrence)
crosstab.columns.name = None
crosstab.to_csv(f'{dout}/lipostatin_recur.csv')

crosstab = pd.crosstab(data.using_hydrostatin, data.biomarker_subtypes)
crosstab.columns.name = None
crosstab.to_csv(f'{dout}/hydrostatin_biocat.csv')

crosstab = pd.crosstab(data.using_hydrostatin, data.recurrence)
crosstab.columns.name = None
crosstab.to_csv(f'{dout}/hydrostatin_recur.csv')
