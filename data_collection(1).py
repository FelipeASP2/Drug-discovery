# Import necessary libraries
import pandas as pd
from chembl_webresource_client.new_client import new_client


# Target search for coronavirus
target = new_client.target
target_query = target.search('coronavirus')
targets = pd.DataFrame.from_dict(target_query)
for index, row in targets.iterrows():
    print(f"{index} {row['pref_name']} - {row['target_chembl_id']}\n")

#Select and retrieve bioactivity data for Replicase polyprotein 1ab
selected_target = targets.target_chembl_id[9]
print(f"{selected_target}\n")

#Here, we will retrieve only bioactivity data for Replicase polyprotein 1ab (CHEMBL4523582) that are reported as IC 50  values in nM (nanomolar) unit
activity = new_client.activity
res = activity.filter(target_chembl_id=selected_target).filter(standard_type="IC50")
df = pd.DataFrame.from_dict(res)
#Showing only the first 3
print(df.head(3))

#Finally we will save the resulting bioactivity data to a CSV file bioactivity_data.csv.
df.to_csv('data/bioactivity_data.csv', index=False)


#If any compounds has missing value for the standard_value and canonical_smiles column then drop it.
df2 = df[df.standard_value.notna()]
df2 = df2[df.canonical_smiles.notna()]
print(df2)
#print(len(df2.canonical_smiles.unique()))
df2_nr = df2.drop_duplicates(['canonical_smiles'])


#The bioactivity data is in the IC50 unit.
# Compounds having values of less than 1000 nM will be considered to be active while those greater than 10,000 nM will be considered to be inactive.
# As for those values in between 1,000 and 10,000 nM will be referred to as intermediate

selection = ['molecule_chembl_id', 'canonical_smiles', 'standard_value']
df3 = df2_nr[selection]
df3.to_csv('data/bioactivity_preprocessed_data.csv', index=False)
df4 = pd.read_csv('data/bioactivity_preprocessed_data.csv')

bioactivity_class = []
for i in df4.standard_value:
  if float(i) >= 10000:
    bioactivity_class.append("inactive")
  elif float(i) <= 1000:
    bioactivity_class.append("active")
  else:
    bioactivity_class.append("intermediate")


bioactivity_series = pd.Series(bioactivity_class, name="bioactivity class")
df5 = pd.concat([df4, bioactivity_series], axis=1)
print(df5)

df5.to_csv('data/bioactivity_curated_data.csv', index=False)
