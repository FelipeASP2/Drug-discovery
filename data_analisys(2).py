import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski

df = pd.read_csv('data/bioactivity_curated_data.csv')
#print(df)


#Smiles to be the last column
df_no_smiles = df.drop(columns='canonical_smiles')
smiles = []
for i in df.canonical_smiles.tolist():
  cpd = str(i).split('.')
  cpd_longest = max(cpd, key = len)
  smiles.append(cpd_longest)

smiles = pd.Series(smiles, name = 'canonical_smiles')

df_clean_smiles = pd.concat([df_no_smiles,smiles], axis=1)
#print(df_clean_smiles)


#Calculate Lipinski descriptors
#Christopher Lipinski, a scientist at Pfizer, came up with a set of rule-of-thumb for evaluating the druglikeness of compounds.
#Such druglikeness is based on the Absorption, Distribution, Metabolism and Excretion (ADME) that is also known as the pharmacokinetic profile.
# Lipinski analyzed all orally active FDA-approved drugs in the formulation of what is to be known as the Rule-of-Five or Lipinski's Rule.

#The Lipinski's Rule stated the following:

# -Molecular weight < 500 Dalton
# -Octanol-water partition coefficient (LogP) < 5
# -Hydrogen bond donors < 5
# -Hydrogen bond acceptors < 10
def lipinski(smiles, verbose=False):

    moldata= []
    for elem in smiles:
        mol=Chem.MolFromSmiles(elem)
        moldata.append(mol)

    baseData= np.arange(1,1)
    i=0
    for mol in moldata:

        desc_MolWt = Descriptors.MolWt(mol)
        desc_MolLogP = Descriptors.MolLogP(mol)
        desc_NumHDonors = Lipinski.NumHDonors(mol)
        desc_NumHAcceptors = Lipinski.NumHAcceptors(mol)

        row = np.array([desc_MolWt,
                        desc_MolLogP,
                        desc_NumHDonors,
                        desc_NumHAcceptors])

        if(i==0):
            baseData=row
        else:
            baseData=np.vstack([baseData, row])
        i=i+1

    columnNames=["MW","LogP","NumHDonors","NumHAcceptors"]
    descriptors = pd.DataFrame(data=baseData,columns=columnNames)

    return descriptors

df_lipinski = lipinski(df_clean_smiles.canonical_smiles)
#print(df_lipinski)


#Combine the two dataframes
df_combined = pd.concat([df,df_lipinski], axis=1)
#print(df_combined)


#Convert IC50 to pIC50
#To allow IC50 data to be more uniformly distributed, we will convert IC50 to the negative logarithmic scale which is essentially -log10(IC50).

#This custom function pIC50() will accept a DataFrame as input and will:

#Take the IC50 values from the standard_value column and converts it from nM to M by multiplying the value by 10 −9
#Take the molar value and apply -log10
#Delete the standard_value column and create a new pIC50 column
def pIC50(input):
    pIC50 = []

    for i in input['standard_value_norm']:
        molar = i*(10**-9) # Converts nM to M
        pIC50.append(-np.log10(molar))

    input['pIC50'] = pIC50
    x = input.drop('standard_value_norm', axis=1)

    return x

#print(df_combined.standard_value.describe())
def norm_value(input):
    norm = []

    for i in input['standard_value']:
        if i > 100000000:
          i = 100000000
        norm.append(i)

    input['standard_value_norm'] = norm
    x = input.drop('standard_value', axis=1)

    return x

df_norm = norm_value(df_combined)
#print(df_norm)

df_final = pIC50(df_norm)
#print(df_final)
#print(df_final.pIC50.describe())

#Removing the 'intermediate' bioactivity class
df_2class = df_final[df_final['bioactivity class'] != 'intermediate']
#print(df_2class)


df_2class.to_csv('data/bioactivity_data_pIC50.csv')