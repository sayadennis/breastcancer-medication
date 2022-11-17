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

## Plot distributions of tumor characteristics summary 

