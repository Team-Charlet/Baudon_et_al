# -*- coding: utf-8 -*-
"""
Created on Fri Dec 22 10:33:48 2023

@author: Angel.BAUDON
"""

import pandas as pd, matplotlib.pyplot as plt, numpy as np, glob, scipy
import seaborn as sns, scipy.stats as stat
from scipy.signal import find_peaks

folder = r"C:\Angel.BAUDON\Exp\Comportement\Freezing\Reanalyse stress\Femelles"
folder_list = [x for x in glob.glob(rf'{folder}\*.xlsx') if not 'analysis' in x]
cond = [x.split('\\')[-1][:-5] for x in folder_list]
show_fig = True

Output = {}
for f, file in enumerate(folder_list):
    
    data = pd.read_excel(file)
    for col in data.columns:
        try: data.rename(columns={col:col.split('; ')[-1]}, inplace=True)
        except: pass
    
    n_stage, sound_len = len(set(data['Stage'])), 20
    n = len(set(data['Animal']))
    output = []
    
    plt.figure(figsize=(20,15))
    for animal in set(data['Animal']):
        datanimal = data[data['Animal'] == animal]
        
        output_mouse = []
        for s, stage in enumerate(datanimal['Stage']):
            anistage = datanimal[datanimal['Stage'] == stage]
            
            X = {}
            for P in [x for x in set(anistage.columns)
                      if x not in ('Animal', 'Stage', 'Test')]:
                X[P] = np.asarray(anistage[P])[0]
            
            cs = X['Cage speaker : time active']
            gzz = X['Floor shock : number']
            freezing = X['Time freezing']
    
            
            CS, _ = find_peaks(np.diff(cs), height=.05, distance=30)
            Gzz, _ = find_peaks(gzz, height=.5, distance=30)
    
    
            output_mouse.append([np.nanmean(freezing[int(cs):int(cs+sound_len)]) for cs in CS])
             
            
            plt.subplot(n, 2, animal*2+s-1)
            plt.plot(freezing, zorder=1), plt.ylabel('Freezing')
            for idx in CS: plt.axvspan(idx, idx+sound_len, color='orange', alpha=.5, zorder=0)
            if len(Gzz): [plt.axvline(gz, color='r') for gz in Gzz]
            plt.savefig(rf'{file[:-5]} all mice.pdf')
        
        output.append(output_mouse)
    plt.close()
    
    
    val = np.asarray([[x[y] for x in output] for y in range(2)])
    
    plt.figure(figsize=(20, 10))
    for i, va in enumerate(val):
        mean, sem = np.nanmean(va, axis=0), scipy.stats.sem(va, axis=0)
            
        plt.subplot(1, 2, i+1)
        plt.plot(mean), plt.xlabel('Expositions'), plt.ylabel('Freezing')
        plt.fill_between(np.arange(len(mean)), mean-sem, mean+sem,
                          color='lightblue', alpha=0.25, zorder=1)
        plt.ylim(0,1)

    Output[cond[f]] = val
    
    plt.savefig(rf'{file[:-5]}.pdf')
    if not show_fig: plt.close()
    

for i, label in enumerate(('Hab', 'Fear')):
    
    vals = [Output[k][i] for k in Output.keys()]
    final_writer = pd.ExcelWriter(f'{folder}\\analysis {label}.xlsx')

    
    plt.figure(), plt.title(label)
    for val, l in zip(vals, cond):
        means, sems = np.nanmean(val, axis=0), stat.sem(val, axis=0)
        plt.errorbar(np.arange(5), means, yerr=sems, label=l)
        plt.legend()
    
        pd.DataFrame(val).to_excel(final_writer, sheet_name=l)
        
    final_writer.save()