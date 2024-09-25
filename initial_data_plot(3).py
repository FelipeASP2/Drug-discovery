import seaborn as sns
sns.set(style='ticks')
import matplotlib.pyplot as plt
import pandas as pd

df_2class = pd.read_csv('data/bioactivity_data_pIC50.csv')

#Frequency plot of the 2 bioactivity classes
plt.figure(figsize=(7, 5.5))

sns.countplot(x='bioactivity class', data=df_2class, hue="bioactivity class", edgecolor='black')

plt.xlabel('Bioactivity class', fontsize=14, fontweight='bold')
plt.ylabel('Frequency', fontsize=14, fontweight='bold')

plt.savefig('data/plot_bioactivity_class.pdf')

#Scatter plot of MW versus LogP
plt.figure(figsize=(7, 5.5))

sns.scatterplot(x='MW', y='LogP', data=df_2class, hue='bioactivity class', size='pIC50', edgecolor='black', alpha=0.7)

plt.xlabel('MW', fontsize=14, fontweight='bold')
plt.ylabel('LogP', fontsize=14, fontweight='bold')
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0)
plt.subplots_adjust(right=0.7)  # Adjust as needed
plt.savefig('data/plot_MW_vs_LogP.pdf')


#BOX PLOTS

#Statistical analysis | Mann-Whitney U Test
def mannwhitney(descriptor, verbose=False):
  # https://machinelearningmastery.com/nonparametric-statistical-significance-tests-in-python/
  from numpy.random import seed
  from numpy.random import randn
  from scipy.stats import mannwhitneyu

# seed the random number generator
  seed(1)

# actives and inactives
  selection = [descriptor, 'bioactivity class']
  df = df_2class[selection]
  active = df[df['bioactivity class'] == 'active']
  active = active[descriptor]

  selection = [descriptor, 'bioactivity class']
  df = df_2class[selection]
  inactive = df[df['bioactivity class'] == 'inactive']
  inactive = inactive[descriptor]

# compare samples
  stat, p = mannwhitneyu(active, inactive)
  #print('Statistics=%.3f, p=%.3f' % (stat, p))

# interpret
  alpha = 0.05
  if p > alpha:
    interpretation = 'Same distribution (fail to reject H0)'
  else:
    interpretation = 'Different distribution (reject H0)'

  results = pd.DataFrame({'Descriptor':descriptor,
                          'Statistics':stat,
                          'p':p,
                          'alpha':alpha,
                          'Interpretation':interpretation}, index=[0])
  filename = 'data/mannwhitneyu_' + descriptor + '.csv'
  results.to_csv(filename)

  return results

#pIC50 value
plt.figure(figsize=(7, 5.5))

sns.boxplot(x = 'bioactivity class', y = 'pIC50', hue='bioactivity class', data = df_2class)

plt.xlabel('Bioactivity class', fontsize=14, fontweight='bold')
plt.ylabel('pIC50 value', fontsize=14, fontweight='bold')

plt.savefig('data/plot_ic50.pdf')

mannwhitney('pIC50')

#MW
plt.figure(figsize=(7, 5.5))

sns.boxplot(x = 'bioactivity class', y = 'MW', hue='bioactivity class', data = df_2class)

plt.xlabel('Bioactivity class', fontsize=14, fontweight='bold')
plt.ylabel('MW', fontsize=14, fontweight='bold')

plt.savefig('data/plot_MW.pdf')

mannwhitney('MW')

#LogP
plt.figure(figsize=(7, 5.5))

sns.boxplot(x = 'bioactivity class', y = 'LogP', hue='bioactivity class', data = df_2class)

plt.xlabel('Bioactivity class', fontsize=14, fontweight='bold')
plt.ylabel('LogP', fontsize=14, fontweight='bold')

plt.savefig('data/plot_LogP.pdf')

mannwhitney('LogP')

#NumHDonors
plt.figure(figsize=(7, 5.5))

sns.boxplot(x = 'bioactivity class', y = 'NumHDonors', hue="bioactivity class", data = df_2class)

plt.xlabel('Bioactivity class', fontsize=14, fontweight='bold')
plt.ylabel('NumHDonors', fontsize=14, fontweight='bold')

plt.savefig('data/plot_NumHDonors.pdf')

mannwhitney('NumHDonors')

#NumHAcceptors
plt.figure(figsize=(7, 5.5))

sns.boxplot(x = 'bioactivity class', y = 'NumHAcceptors',  hue="bioactivity class", data = df_2class)

plt.xlabel('Bioactivity class', fontsize=14, fontweight='bold')
plt.ylabel('NumHAcceptors', fontsize=14, fontweight='bold')

plt.savefig('data/plot_NumHAcceptors.pdf')

mannwhitney('NumHAcceptors')