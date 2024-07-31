# -*- coding: utf-8 -*-
"""
Created on Thu Jan 26 10:04:54 2023

@author: Angel.BAUDON
"""


import numpy as np, pandas as pd, matplotlib.pyplot as plt, scipy, seaborn as sns, scipy


ctrl = r"C:\Angel.BAUDON\Exp\Comportement\Social Interaction\Ctrl.xlsx"
ko = r"C:\Angel.BAUDON\Exp\Comportement\Social Interaction\GFAP OTR KO.xlsx"
labels = ("KO", "Ctrl")
Param = ('Distance', 'SI : time', 'Loin : time', 'Balek : time')
Out = pd.DataFrame()


for data, name in zip((ko, ctrl), labels):
    data = pd.read_excel(data)
    
    for col in data.columns:
        try: data.rename(columns={col:col.split('; ')[-1]}, inplace=True)
        except: pass
    
    fem, mal = [data[data['Test'] == x] for x in ('FEM', 'MAL')]
    
    for s, sex in zip((' FEM', ' MAL'), (fem, mal)):
        out_sex = pd.DataFrame(columns=('Grp', 'Sex', 'ID', *Param))
        
        for a, animal in sex.iterrows():
            out_animal = {'Distance': sum(animal.loc['Distance']),'Sex': s,
                          'Grp': name, 'ID': animal.loc['Animal']}
            
            for param in Param[1:]:
                p = animal.loc[param]
                splt = (sum(p[:len(p)//2]), sum(p[len(p)//2:]))
                delta = splt[1] - splt[0]
                out_animal[param] = delta
            out_sex.loc[a] = out_animal
            
        Out = pd.concat((Out, out_sex))
        
                
plt.figure()
sns.barplot(data=Out, x='Grp', y='Distance', hue='Sex')
