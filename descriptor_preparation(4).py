import pandas as pd

df3 = pd.read_csv('data/bioactivity_data_pIC50.csv')
#print(df3)

selection = ['canonical_smiles','molecule_chembl_id']
df3_selection = df3[selection]
df3_selection.to_csv('molecule.smi', sep='\t', index=False, header=False)
#Get-Content molecule.smi | Select-Object -First 5
#(Get-Content molecule.smi).Count


#Calculate fingerprint descriptors
#Calculate PaDEL descriptors
#bash padel.sh


#Preparing the X and Y Data Matrices
#X data matrix

df3_X = pd.read_csv('data/descriptors_output.csv')
#print(df3_X)
df3_X = df3_X.drop(columns=['Name'])
#print(df3_X)

#Y variable
#Convert IC50 to pIC50
df3_Y = df3['pIC50']
#print(df3_Y)


#Combining X and Y variable
dataset3 = pd.concat([df3_X,df3_Y], axis=1)
#print(dataset3)

dataset4 = dataset3.dropna()


dataset4.to_csv('data/bioactivity_data_pIC50_pubchem_fp.csv', index=False)


