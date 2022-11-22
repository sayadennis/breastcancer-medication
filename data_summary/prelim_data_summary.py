import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

## Load data 
dn = '/share/fsmresfiles/breastcancer_medication/'

data = pd.read_csv(f'{dn}/6800_cohort.csv', header=None)

with open(f'{dn}/6800_cohort_colnames.txt', 'r') as f:
    lines = f.readlines()

colnames = [line.strip() for line in lines]
colnames.remove('diagnosis')
for i in range(len(colnames)):
    if colnames[i]=='total_mamography_before':
        colnames[i] = 'total_mamography_before_diagnosis'

data.columns = colnames

drug_use = pd.read_csv(f'{dn}/test_6800_metformin_use.csv', header=None)

with open(f'{dn}/test_6800_metformin_use_colnames.txt', 'r') as f:
    lines = f.readlines()

colnames = [line.strip() for line in lines]
colnames.remove('diagnosis')
for i in range(len(colnames)):
    if colnames[i]=='total_mamography_before':
        colnames[i] = 'total_mamography_before_diagnosis'

drug_use.columns = colnames

## Get recurrence & drug use categorical data 
ptid_using_metformin = drug_use.iloc[~pd.isnull(drug_use.generic_name).values,:].patient_ir_id.unique()
using_metfo = [x in ptid_using_metformin for x in data.patient_ir_id.values]
local_recur = (data.Local_Recurrence=='YES').values
dista_recur = (data.Distant_Recurrence=='YES').values

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

biol_cat = []
for i in data.index:
    if ((data.loc[i,'ER']=='POSITIVE') & (data.loc[i,'PR']=='POSITIVE') & (data.loc[i,'HER2']=='NEGATIVE')):
        biol_cat.append('ER/PR+,HER2-')
    elif data.loc[i,'HER2']=='POSITIVE':
        biol_cat.append('HER2+')
    else:
        biol_cat.append('TN')

## Summary count tables
data['metformin_use'] = using_metfo
data['local_recurrence'] = local_recur
data['distant_recurrence'] = dista_recur
data['recurrence'] = recur
data['biologic_category'] = biol_cat

pd.crosstab(data.metformin_use, data.biologic_category)
pd.crosstab(data.metformin_use, data.recurrence)

## Plot distributions of tumor characteristics summary 
