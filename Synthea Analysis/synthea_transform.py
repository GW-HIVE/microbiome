#Python Script to Transform Raw Synthea COVID-19 Dataset:
#Author(s): James Ziegler, Alex Coleman
#Date: 09/20/2021
#Version: 1.0
#Libraries Needed:
import pandas as pd
import numpy as np
import datetime
from numpy.core.numeric import NaN
import os
import time
from functools import reduce
from datetime import datetime, date
#Importing Datasets with Pandas:
conditions = pd.read_csv("conditions.csv")
patients = pd.read_csv("patients.csv")
observations = pd.read_csv("observations.csv")
care_plans = pd.read_csv("careplans.csv")
encounters = pd.read_csv("encounters.csv")
devices = pd.read_csv("devices.csv")
supplies = pd.read_csv('supplies.csv')
procedures = pd.read_csv("procedures.csv")
medications = pd.read_csv("medications.csv")
#Creating arrays of unique patient Ids for patients with COVID, Diabetes, and Prediabetes to reference later:
covid_patient_ids = conditions[conditions.CODE == 840539006].PATIENT.unique()
diabetes_patient_ids = conditions[conditions.CODE == 44054006].PATIENT.unique()
prediabetes_patient_ids = conditions[conditions.CODE == 15777000].PATIENT.unique()
negative_covid_patient_ids = observations[(observations.CODE == '94531-1') & 
                                          (observations.VALUE == 'Not detected (qualifier value)')
                                          ].PATIENT.unique()
deceased_patients = patients[patients.DEATHDATE.notna()].Id
completed_isolation_patients = care_plans[(care_plans.CODE == 736376001) & (care_plans.STOP.notna()) & (care_plans.REASONCODE == 840539006)].PATIENT
survivor_ids = np.union1d(completed_isolation_patients, negative_covid_patient_ids)
inpatient_ids = encounters[(encounters.REASONCODE == 840539006) & (encounters.CODE == 1505002)].PATIENT
isolation_ids = care_plans[(care_plans.CODE == 736376001) & (care_plans.REASONCODE == 840539006)].PATIENT
icu_ids = encounters[encounters.CODE == 305351004].PATIENT
vent_ids = procedures[procedures.CODE == 26763009].PATIENT
#Grabbing lab observation values for six biomarkers:
#(Note: to add more biomarkers, expand the code here)
lab_obs = observations[(observations.CODE == '48065-7') | (observations.CODE == '26881-3') | 
                       (observations.CODE == '89579-7') | (observations.CODE == '731-0') |
                       (observations.CODE == '1988-5') | (observations.CODE == '751-8') | 
                       (observations.CODE == '14804-9') |(observations.CODE == '2276-4')]
#Indexing the conditions dataset to just COVID patients:
covid_conditions = conditions[(conditions.CODE == 840539006)]
#Merging all of the COVID conditions with the total patients dataframe:
covid_patients = covid_conditions.merge(patients, how='left', left_on='PATIENT', right_on='Id')
#Create new columns within the COVID Patients dataframe:
covid_patients['Diabetes'] = covid_patients.PATIENT.isin(diabetes_patient_ids)
covid_patients['PreDiabetes'] = covid_patients.PATIENT.isin(prediabetes_patient_ids)
covid_patients['Survivor'] = covid_patients.PATIENT.isin(survivor_ids)
covid_patients['HomeIsolation'] = covid_patients.PATIENT.isin(isolation_ids)
covid_patients['HospitalAdmit'] = covid_patients.PATIENT.isin(inpatient_ids)
covid_patients['ICUAdmit'] = covid_patients.PATIENT.isin(icu_ids)
covid_patients['Ventilated'] = covid_patients.PATIENT.isin(vent_ids)
#Merging the lab observations for the six biomarkers on the COVID patient information dataframe:
covid_patients_obs = covid_patients.merge(lab_obs, on='PATIENT')
covid_patients_obs['START'] = pd.to_datetime(covid_patients_obs.START)
covid_patients_obs['DATE'] = pd.to_datetime(covid_patients_obs.DATE)
covid_patients_obs['lab_days'] = covid_patients_obs.DATE - covid_patients_obs.START
covid_patients_obs['days'] = covid_patients_obs.lab_days / np.timedelta64(1, 'D')
covid_patients_obs = covid_patients_obs.drop('lab_days', axis=1)
covid_patients_obs['VALUE'] = pd.to_numeric(covid_patients_obs['VALUE'], errors='coerce')
#Setting up the raw machine learning table:
ml_table = covid_patients_obs[['PATIENT',
                               'Diabetes', 
                               'PreDiabetes',
                               'HomeIsolation', 
                               'HospitalAdmit',
                               'ICUAdmit',
                               'Ventilated',
                               'Survivor',
                               'DESCRIPTION_y',
                               'VALUE',
                               'days']]
ml_table = ml_table.rename(columns={"DESCRIPTION_y":"Biomarker"})
ml_table.to_csv('ml_table.csv')
#Transforming the raw machine learning table to be patient centric:
df = pd.read_csv('ml_table.csv', index_col=False)
grouped = df.groupby('Biomarker')
for biomarker in grouped.groups:

    # Get the group
    temp_table = grouped.get_group(biomarker)
    
    # Pivot the group such that days becomes the column, with the value being VALUE, and the index being all other columns
    patient_grouped = pd.pivot_table(temp_table, index = list(temp_table.columns)[1:-2], columns = ["days"], values = ["VALUE"])
    # https://stackoverflow.com/questions/33290374/pandas-pivot-table-column-names
    
    # rename day columns to have "day_" appended in front 
    patient_grouped.columns = ["day_" + str(s2) for (s1,s2) in patient_grouped.columns.tolist()]

    # for each day in range 0, 27 that doesn't exist, create column and fill with NaN values
    for i in range (0, 27):
        if "day_" + str(i) + ".0" not in patient_grouped.columns:
            patient_grouped.insert(i, "day_" + str(i) + ".0", NaN)

    # Reset index to put other columns inplace with day columns
    patient_grouped.reset_index(inplace=True)

    # Save dataframe in CSV in biomarker directory, with (partial) name of biomarker as filename
    filename = os.path.join(biomarker.split(" [")[0].replace(" ", "-") + ".csv")
    patient_grouped.to_csv(filename,index=False,na_rep="NA")
#Allowing the script to pause in order to load in the new dataframes:
time.sleep(7)
#Loading in patient-centric daily biomarker data and other related patient data:
d_d = pd.read_csv('Fibrin-D-dimer-FEU.csv')
crp = pd.read_csv('C-reactive-protein.csv')
il6 = pd.read_csv('Interleukin-6.csv')
lymp = pd.read_csv('Lymphocytes.csv')
neut = pd.read_csv('Neutrophils.csv')
tnn = pd.read_csv('Troponin-I.cardiac.csv')
lac_deh = pd.read_csv('Lactate-dehydrogenase.csv')
ferr = pd.read_csv('Ferritin.csv')
#Renaming the columns:
d_d = d_d.replace('Fibrin D-dimer FEU [Mass/volume] in Platelet poor plasma', 'D-dimer_BID:A0020')
crp = crp.replace('C reactive protein [Mass/volume] in Serum or Plasma', 'C_reactive_protien_BID:A0009')
il6 = il6.replace('Interleukin 6 [Mass/volume] in Serum or Plasma','Interleukin_6_BID:A0001')
lymp = lymp.replace('Lymphocytes [#/volume] in Blood by Automated count', 'Lymphocytes_BID:A0004')
neut = neut.replace('Neutrophils [#/volume] in Blood by Automated count', 'Neutrophils_BID:A0055')
tnn = tnn.replace('Troponin I.cardiac [Mass/volume] in Serum or Plasma by High sensitivity method','Troponin_cardiac_BID:A0029')
lac_deh = lac_deh.replace('Lactate dehydrogenase [Enzymatic activity/volume] in Serum or Plasma by Lactate to pyruvate reaction', 'Lactate_dehydrogenase_BID:A0021')
ferr = ferr.replace('Ferritin [Mass/volume] in Serum or Plasma','Ferritin_BID:A0006')
#Creating a seperate dataframe that exludes biomarker data to use to merge later on:
patient_stats = d_d[['PATIENT',
                     'Diabetes',
                     'PreDiabetes',
                     'HomeIsolation',
                     'HospitalAdmit',
                     'ICUAdmit',
                     'Ventilated',
                     'Survivor']].copy()
patient_stats = patient_stats.rename(columns={'Diabetes':'Diabetes_DOID:9352',
                                              'PreDiabetes':'PreDiabetes_DOID:11716',
                                              })
patient_dob = patients.filter(['Id','BIRTHDATE'], axis=1)
patient_dob = patient_dob.rename(columns={'Id':'PATIENT'})
patient_stats = patient_stats.merge(patient_dob, on='PATIENT')
#Adding patient age:
def age(born):
    born = datetime.strptime(born, "%Y-%m-%d").date()
    today = date.today()
    return today.year - born.year - ((today.month, 
                                      today.day) < (born.month, 
                                                    born.day))
patient_stats['Age'] = patient_stats['BIRTHDATE'].apply(age)
del patient_stats['BIRTHDATE']
#Adding biomarker names to column headers, and removing information that was put into "patient_stats" dataframe:
dataframes = [d_d, crp, il6, lymp, neut, tnn, lac_deh, ferr]
for df in dataframes:
    cols = df.columns[~df.columns.isin(['PATIENT',
                                        'Diabetes',
                                        'PreDiabetes',
                                        'HomeIsolation',
                                        'HospitalAdmit',
                                        'ICUAdmit',
                                        'Ventilated',
                                        'Survivor',
                                        'Biomarker'])]
    df.rename(columns = dict(zip(cols, cols + '_' + df.Biomarker[1])), inplace=True)
    del df['Diabetes']
    del df['PreDiabetes']
    del df['HomeIsolation']
    del df['HospitalAdmit']
    del df['ICUAdmit']
    del df['Ventilated']
    del df['Survivor']
    del df['Biomarker']
#Merging all biomarker levels and patient stats:
dataframes = [patient_stats, d_d, crp, lymp, tnn, lac_deh, ferr]
df_merged = reduce(lambda  left,right: pd.merge(left,right,on=['PATIENT',],
                                            how='outer'), dataframes).fillna('NA')
#Extracting Day Zero (admission day) data and joining it with the patient stats:
day_0 = df_merged.filter(like='day_0.0_', axis=1)
day_0 = patient_stats.join(day_0, how='outer')
#Exporting ML-ready datasets:
day_0.to_csv('day_0_levels.csv', index=False)
df_merged.to_csv('day_all_levels.csv', index=False)
