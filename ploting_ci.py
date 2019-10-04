# -*- coding: utf-8 -*-
"""
Created on Fri Aug 23 14:19:01 2019

@author: vilmarith
"""

import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

%matplotlib qt
sns.set(font_scale=2)

os.getcwd()
os.chdir('D:\\PhD\Draft Papers\Current\degeneracy')

df = pd.read_csv('red_deg_output2.csv')
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
df.columns = ['Number of Correct Reponses (%)', 'Group Type \n (n=50)', 'Redundant Specification']


sns.set_style("white")

meanlineprops = dict(linestyle='-', linewidth=1.5, color='black')
meanmarker=dict(marker='d',markerfacecolor="dimgrey", markeredgecolor="dimgrey", markersize=10)

fig, ax = plt.subplots()
sns.boxplot(x=df['Group Type \n (n=50)'], y=df["Number of Correct Reponses (%)"],   color="lightsteelblue", 
                 hue = df['Redundant Specification'], 
             meanline=False, showmeans=False, linewidth=1.5 )
sns.despine(offset=10, trim=True)

ax.set_ylim(0,110)
plt.show()
