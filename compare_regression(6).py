import pandas as pd
import seaborn as sns
from sklearn.model_selection import train_test_split
from lazypredict.Supervised import LazyRegressor
import matplotlib.pyplot as plt
import seaborn as sns


df = pd.read_csv('data/bioactivity_data_pIC50_pubchem_fp.csv')
X = df.drop('pIC50', axis=1)
Y = df.pIC50

# Examine X dimension
#print(X.shape)

# Remove low variance features
from sklearn.feature_selection import VarianceThreshold
selection = VarianceThreshold(threshold=(.8 * (1 - .8)))
X = selection.fit_transform(X)
#print(X.shape)

# Perform data splitting using 80/20 ratio
X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=0.2, random_state=42)

# Defines and builds the lazyclassifier
clf = LazyRegressor(verbose=0,ignore_warnings=True, custom_metric=None)
models_train,predictions_train = clf.fit(X_train, X_train, Y_train, Y_train)
models_test,predictions_test = clf.fit(X_train, X_test, Y_train, Y_test)
print(predictions_train)
# Performance table of the test set (20% subset)
#print(predictions_test)

#Data visualization of model performance
#train["R-Squared"] = [0 if i < 0 else i for i in train.iloc[:,0] ]
palette = sns.color_palette("husl", len(predictions_train.index))
# Bar plot of R-Squared values
plt.figure(figsize=(8, 10))
sns.set_theme(style="whitegrid")
ax = sns.barplot(y=predictions_train.index, x="R-Squared", data=predictions_train, palette=palette)
ax.set(xlim=(0, 1))
plt.yticks(fontsize=12)  # Rotate and adjust font size
plt.tight_layout()  # Adjust layout
plt.savefig('data/r_squared_plot.pdf')
plt.clf()

# Bar plot of RMSE values
plt.figure(figsize=(8, 10))
sns.set_theme(style="whitegrid")
ax = sns.barplot(y=predictions_train.index, x="RMSE", data=predictions_train, palette=palette)
ax.set(xlim=(0, 10))
plt.yticks(fontsize=12)  # Rotate and adjust font size
plt.tight_layout()  # Adjust layout
plt.savefig('data/rmse_plot.pdf')
plt.clf()

# Bar plot of calculation time
plt.figure(figsize=(8, 10))
sns.set_theme(style="whitegrid")
ax = sns.barplot(y=predictions_train.index, x="Time Taken", data=predictions_train, palette=palette)
ax.set(xlim=(0, 10))
plt.yticks(fontsize=12)  # Rotate and adjust font size
plt.tight_layout()  # Adjust layout
plt.savefig('data/calculation_time_plot.pdf')
plt.clf()
