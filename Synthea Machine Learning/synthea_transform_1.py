#Python Script to Transform Raw Synthea COVID-19 Dataset (1/2):
#Author(s): James Ziegler, Alex Coleman
#Date: 09/20/2021
#Version: 1.0
#Libraries Needed:
import pandas as pd
import numpy as np
import datetime
from numpy.core.numeric import NaN
import os
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
biomarkerDir = os.path.join(os.path.abspath(os.getcwd()), "biomarkerFiles")
if not os.path.isdir('biomarkerFiles'):
    os.mkdir(biomarkerDir)
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

    # Remove column containing biomarker
    #del patient_grouped["Biomarker"]

    # Save dataframe in CSV in biomarker directory, with (partial) name of biomarker as filename
    filename = os.path.join(biomarkerDir, biomarker.split(" [")[0].replace(" ", "-") + ".csv")
    patient_grouped.to_csv(filename,index=False,na_rep="NA")
                       
