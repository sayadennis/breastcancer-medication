import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib

## Load data 
din = '/share/fsmresfiles/breastcancer_medication/data/01_ssms'
dout = '/share/fsmresfiles/breastcancer_medication/plots'

data = pd.read_csv(f'{din}/6800_cohort.csv', header=None)

with open(f'{din}/6800_cohort_colnames.txt', 'r') as f:
    lines = f.readlines()

colnames = [line.strip() for line in lines]
colnames.remove('diagnosis')
for i in range(len(colnames)):
    if colnames[i]=='total_mamography_before':
        colnames[i] = 'total_mamography_before_diagnosis'

data.columns = colnames

drug_use = pd.read_csv(f'{din}/test_6800_metformin_use.csv', header=None)

with open(f'{din}/test_6800_metformin_use_colnames.txt', 'r') as f:
    lines = f.readlines()

# colnames = [line.strip() for line in lines]
colnames = lines[0].split()
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

## Plot distributions of features of interest, separated by recurrence status (outcome) 

# Age at diagnosis 
fig, ax = plt.subplots(5,1, figsize=(4,8))
cmap = matplotlib.cm.get_cmap('Pastel1')
ax[0].hist(data.age_at_diagnosis, bins=20, color='gray', label='All')
ax[0].set_xlim((np.min(data.age_at_diagnosis), np.max(data.age_at_diagnosis)))
# ax[0].set_title(f'All patients combined')
ax[0].legend()

for i, recur_status in enumerate(['None', 'Local', 'Distant', 'Both']):
    ax[i+1].hist(
        data.iloc[data.recurrence.values==recur_status,:].age_at_diagnosis, 
        bins=20, color=cmap(i), label=recur_status # , alpha=0.25
    )
    ax[i+1].set_xlim((np.min(data.age_at_diagnosis), np.max(data.age_at_diagnosis)))
    ax[i+1].legend()
    # ax[i+1].set_title(f'Recurrence status: {recur_status}')

ax[2].set_ylabel('Number of patients')
ax[4].set_xlabel('Age')
fig.suptitle('Age at Dx: by Recurrence Status')
plt.tight_layout()
fig.savefig(f'{dout}/histogram_age_at_dx.png')
plt.close()

# biologic subtype 
pd.crosstab(data.biologic_category, data.recurrence)

# histology 
def hist_cat(row):
    if row['Histology']=='DUCT':
        if row['Invasive']=='YES':
            return 'IDC'
        else:
            return 'DCIS'
    elif row['Histology']=='LOBULAR':
        if row['Invasive']=='YES':
            return 'ILC'
        else:
            return 'LCIS'
    else:
        if row['Invasive']=='YES':
            return 'Mixed Invasive'
        else:
            return 'Mixed Noninvasive'

data['histology_category'] = data.apply(hist_cat, axis=1)
pd.crosstab(data.histology_category, data.recurrence)

# drug (metformin, statin) use status (as determined by current method) 
