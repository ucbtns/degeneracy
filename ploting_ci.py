import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

%matplotlib qt
sns.set(font_scale=2)


df = pd.read_csv('red_deg.csv')
df['Group Type'] = ''
n = 50

df.iloc[0:n, 1] = 'Control'
df.iloc[n:n*2, 1] = 'Control'
df.iloc[n*2:n*3, 1] = 'B Lesion'
df.iloc[n*3:n*4, 1] = 'B Lesion'
df.iloc[n*4:n*5, 1] = 'A lesion'
df.iloc[n*5:n*6, 1] = 'A lesion'
df.iloc[n*6:n*7, 1] = 'A & B Lesion'
df.iloc[n*7:n*8, 1] = 'A & B Lesion'


df['Redundant Specification'] = ''
df.iloc[0:n, 2] = 'Y'
df.iloc[n:n*2, 2] = 'N'
df.iloc[n*2:n*3, 2] = 'Y'
df.iloc[n*3:n*4, 2] = 'N'
df.iloc[n*4:n*5, 2] = 'Y'
df.iloc[n*5:n*6, 2] = 'N'
df.iloc[n*6:n*7, 2] = 'Y'
df.iloc[n*7:n*8, 2] = 'N'

df.iloc[:,0] = df.iloc[:,0]*100
df.columns = ['Behavioural Accuracy (%)', 'Group Type \n (n=50)', 'Duplicate Specification']


df2 = pd.read_csv('accur.csv')
df2['Duplicate Specification'] = ['Y', 'N']*4
df2['Group Type \n (n=50)'] = ['Control']*2 + ['B Lesion']*2 + ['A Lesion']*2 + ['A & B Lesion']*2

sns.set_style("white")

meanlineprops = dict(linestyle='-', linewidth=1.5, color='black')
meanmarker=dict(marker='d',markerfacecolor="dimgrey", markeredgecolor="dimgrey", markersize=10)

fig, ax = plt.subplots()
sns.boxplot(x=df['Group Type \n (n=50)'], y=df["Behavioural Accuracy (%)"],   color="lightsteelblue", 
                 hue = df['Duplicate Specification'], 
             meanline=False, showmeans=False, linewidth=1.5 )
sns.despine(offset=10, trim=True)
ax.set_ylim(0,110)

ax2 = ax.twinx()  # instantiate a second axes that shares the same x-axis

color = 'white'
ax2.set_ylabel('Statistical Accuracy (nats)')  # we already handled the x-label with ax1
ax2.plot(df2['Group Type \n (n=50)'], df2['Var1'],color=color)

plt.show()


df3 = pd.read_excel('tabl2.xlsx')
df3 = df3[df3['Duplicate'] == 'N']

fig, ax = plt.subplots()
sns.barplot(x= df3['Quantity'], y=df3["Value"],   color="lightsteelblue", 
                 hue =  df3['Group'])

df4 = pd.read_csv('factor_breakdown.csv')
df4 = df4[df4['Quantity'] == 'Redundancy']

fig, ax = plt.subplots()
sns.barplot(x=df4['Factor'], y=df4["Natural Units (nats)"],   color="lightsteelblue", 
                 hue = df4['Specification'])
sns.despine(offset=10, trim=True)
plt.show()

df4 = pd.read_csv('factor_breakdown.csv')
df4.sort_values(by = 'Quantity')
df4 = df4[df4['Quantity'] == 'Free Energy']

fig, ax = plt.subplots()
sns.barplot(x=df4['Hidden State'], y=df4["Natural Units (nats)"],  color="lightsteelblue", 
                 hue = df4['Duplicate Specification'])
sns.despine(offset=10, trim=True)
ax.set_ylim(0,3.5)
plt.show()
