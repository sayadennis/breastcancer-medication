import numpy as np
import pandas as pd
from scipy.stats import chisquare

## Load data 
dn = '/share/fsmresfiles/breastcancer_medication/'

data = pd.read_csv(f'{dn}/6800_cohort.csv', header=None)
drug_use = pd.read_csv(f'{dn}/test_6800_metformin_use.csv', header=None)

## Get necessary categorical data 

## Perform statistical test
