import pandas as pd
import seaborn as sns
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.model_selection import cross_val_score
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import VarianceThreshold


df = pd.read_csv('data/bioactivity_data_pIC50_pubchem_fp.csv')


# 2. Remove low-variance features
selection = VarianceThreshold(threshold=(.8 * (1 - .8)))
X_clean = selection.fit_transform(df.drop('pIC50', axis=1))
Y_clean = df['pIC50']

# 3. Scale the features
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X_clean)

# 4. Split the data
X_train, X_test, Y_train, Y_test = train_test_split(X_scaled, Y_clean, test_size=0.2, random_state=100)

# 5. Train the model using cross-validation

model = RandomForestRegressor(n_estimators=100, random_state=100)
cv_scores = cross_val_score(model, X_train, Y_train, cv=5, scoring='r2')
print(f'Mean R²: {cv_scores.mean()}')
print(f'Standard Deviation: {cv_scores.std()}')


# 6. Train on the entire training set and evaluate on the test set
model.fit(X_train, Y_train)
r2_test = model.score(X_test, Y_test)
print("Test R²:", r2_test)
Y_pred = model.predict(X_test)

#Scatter Plot of Experimental vs Predicted pIC50 Values
sns.set(color_codes=True)
sns.set_style("white")

ax = sns.regplot(x=Y_test, y=Y_pred, scatter_kws={'alpha': 0.4})
ax.set_xlabel('Experimental pIC50', fontsize='large', fontweight='bold')
ax.set_ylabel('Predicted pIC50', fontsize='large', fontweight='bold')
ax.set_xlim(0, 12)
ax.set_ylim(0, 12)
ax.figure.set_size_inches(5, 5)
plt.savefig('data/plot_regression_model.pdf')
