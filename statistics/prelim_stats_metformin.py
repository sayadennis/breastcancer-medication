import numpy as np
import pandas as pd
from scipy.stats import chisquare

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

## Get necessary categorical data 
ptid_using_metformin = drug_use.iloc[~pd.isnull(drug_use.generic_name).values,:].patient_ir_id.unique()
using_metfo = [x in ptid_using_metformin for x in drug_use.patient_ir_id.values]
local_recur = (drug_use.Local_Recurrence=='YES').values
dista_recur = (drug_use.Distant_Recurrence=='YES').values

## Perform statistical test
