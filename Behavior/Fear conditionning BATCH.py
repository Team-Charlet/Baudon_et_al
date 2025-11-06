# -*- coding: utf-8 -*-
"""
Created on Tue Jul 19 11:57:37 2022

@author: Angel.BAUDON
"""

import pandas as pd, matplotlib.pyplot as plt, numpy as np, glob, scipy, pandas as pd
from scipy.signal import find_peaks

folder = r"C:\Users\33666\Desktop\GFAP OTR KO FC Males.pdf"
files = glob.glob(rf'{folder}\*.xlsx')
Yggdrasil = {}
show_fig = False

for file in files:
    drug = file.split('\\')[-1].split('.')[0]
    print(drug)
    
    data = pd.read_excel(file)
    for col in data.columns:
        try: data.rename(columns={col:col.split('; ')[-1]}, inplace=True)
        except: pass
    
    n_stage, sound_len = len(set(data['Stage'])), 20
    Parameters = ('Freezing', 'Speed', 'PSTH Freezing', 'PSTH speed',
                  'Distance', 'Center', 'Periphery')
    Output = {x:[] for x in Parameters}
    
    for animal in set(data['Animal']):
        datanimal = data[data['Animal'] == animal]
        
        Output_mouse = {x:[] for x in Parameters}
        plt.figure(), plt.title(f'Mouse {animal}')
        for s, stage in enumerate(datanimal['Stage']):
            anistage = datanimal[datanimal['Stage'] == stage]
            
            X = {}
            for P in [x for x in set(anistage.columns)
                      if x not in ('Animal', 'Stage', 'Test')]:
                X[P] = np.asarray(anistage[P])[0]
            cs, gzz,  = X['Cage speaker : time active'], X['Floor shock : number']
            speed, freezing = [X[x] for x in ('Mean speed', 'Time freezing')]
            periph, center = [X[x] for x in ('Periphery : time', 'Center : time')]
            distance = X['Distance']
            
            CS, _ = find_peaks(np.diff(cs), height=.05, distance=30)
            Gzz, _ = find_peaks(gzz, height=.5, distance=30)
    
    
            for a, b in zip(('Freezing', 'Speed'), (freezing, speed)):
                Output_mouse[a].append([np.nanmean(b[int(cs+sound_len) : int(cs+sound_len*1.5+1)]) for cs in CS])
            
            Output_mouse['Distance'].append(np.nansum(distance))
            for a, b in zip(('Center', 'Periphery'), (center, periph)):
                Output_mouse[a].append(np.nansum(b)/len([x for x in b if str(x) != 'nan']))
    
            red_cs = CS[2:] if 'FC' in stage else CS[:3]
            for a, b in zip(('PSTH Freezing', 'PSTH speed'), (freezing, speed)):
                Output_mouse[a].append([b[cs+1-sound_len : cs+(sound_len*2)] for cs in red_cs])
            
    
            rows, loc = n_stage*2+1, (n_stage*2)**2
            for i, (j, label) in enumerate(zip((freezing, speed), ('Freezing', 'Speed'))):
                plt.subplot(rows,1,s+(n_stage*i)+1)
                plt.plot(j, zorder=1), plt.ylabel(label)
                for idx in CS: plt.axvspan(idx, idx+sound_len, color='orange', alpha=.5, zorder=0)
                if len(Gzz): [plt.axvline(gz, color='r') for gz in Gzz]
    
        for k, l in enumerate((Output_mouse['Freezing'], Output_mouse['Speed'])):
            for i, j in enumerate(l):
                plt.subplot(rows, n_stage*2, loc+1+i+(n_stage*k))
                plt.ylim(0,1), plt.plot(j, c='r')
        if not show_fig: plt.close()
        
        for k in Output.keys(): Output[k].append(Output_mouse[k])
    
    
    col, n = len(set(data['Stage'])), len(set(data['Animal']))
    for k in (('Freezing', 'PSTH Freezing'), ('Speed', 'PSTH speed')):
        val = [[x[y] for x in Output[k[0]]] for y in range(col)]
        psth = [[np.nanmean(np.asarray(x[y]), axis=0) for x in Output[k[1]]]
                for y in range(col)]
        plt.figure()
        for i, (va, pst) in enumerate(zip(val, psth)):
            ar = np.zeros((n, max([len(x) for x in va])))
            ar[:], pst = np.nan, np.asarray(pst)
            
            for y, v in enumerate(va):
                for z, d in enumerate(v): ar[y,z] = d
    
            for a, (b, label) in enumerate(zip((ar, pst), ('Expositions', 'Seconds'))):
                mean, sem = np.nanmean(b, axis=0), scipy.stats.sem(b, axis=0)
                
                plt.subplot(2, col, col*a+i+1)
                plt.plot(mean, c='b')
                plt.fill_between(np.arange(len(mean)), mean-sem, mean+sem,
                                 color='lightblue', alpha=0.25, zorder=1)
                
                plt.xlabel(label), plt.ylabel(k[a])
                if a: plt.axvspan(sound_len, sound_len*2, color='orange', alpha=.5)
                plt.ylim(0,1) if k[0] == 'Freezing' else plt.ylim(0,.2)
        # if not show_fig: plt.close()
    
    for k in ('Distance', 'Center', 'Periphery'):
        plt.figure(), plt.title(k)
        import seaborn as sns
        df = pd.DataFrame(data = Output[k], columns = list(datanimal['Stage']))
        
        sns.barplot(data=df, color='lightgreen', capsize=.2)
        sns.swarmplot(data=df, size=5, color='w', edgecolor='k', linewidth=1)
        if not show_fig: plt.close()
        
    Yggdrasil[drug] = Output


df = pd.DataFrame(columns=('Treatment', 'Stage', 'Animal', 'Exposition', 'Score'))
plt.figure()
for i, drug in enumerate(Yggdrasil):
    Ygg, stages = Yggdrasil[drug], set(data['Stage'])
    col, n = len(stages), len(Ygg['Center'])
    val = [[x[y] for x in Ygg['Freezing']] for y in range(col)]
    
    for i, (va, stage) in enumerate(zip(val, stages)):
        ar = np.zeros((n, max([len(x) for x in va])))
        ar[:] = np.nan
        
        for animal_ID, animal_data in enumerate(va):
            for expo_ID, expo_data in enumerate(animal_data):
                ar[animal_ID,expo_ID] = expo_data
                df.loc[len(df.index)] = (drug, stage, animal_ID, expo_ID, expo_data)
    
        mean, sem = np.nanmean(ar, axis=0), scipy.stats.sem(ar, axis=0)
        plt.subplot(1, col, i+1), plt.plot(mean, label=drug), plt.legend()
        plt.fill_between(np.arange(len(mean)), mean-sem, mean+sem,
                         color='lightblue', alpha=0.25, zorder=1)
        plt.xlabel('Expositions'), plt.ylabel('Freezing'), plt.ylim(0,1)
        
 plt.figure()
 sns.barplot(data=df, x='Treatment', y='Score', hue='Stage', color='lightgreen', capsize=.2)
 sns.swarmplot(x='Treatment', y='Score', hue='Stage', dodge=True, data=df, size=5, color='w',
               edgecolor='k', linewidth=1)
