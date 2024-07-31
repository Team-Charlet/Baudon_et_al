# -*- coding: utf-8 -*-
"""
Created on Fri Dec 22 13:41:37 2023

@author: Angel.BAUDON
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Mar 14 17:37:52 2023

@author: Angel.BAUDON
"""
import pyabf, matplotlib.pyplot as plt, numpy as np, glob, pandas as pd, os
import scipy.stats as stat, seaborn as sns
from scipy.signal import find_peaks
from pyabf.filter import gaussian
from ToolKit.Tollbox_TG_shorten import  toto_filter

folder = r"C:\Angel.BAUDON\Exp\Patch\Data\2023_Puff OTR Ago on OTR GFP"
sub_folders = [x for x in glob.glob(rf'{folder}\*')
               if x.split('\\')[-1] not in ('analysis', 'ISteps', 'Raw')]
if not os.path.exists(rf'{folder}\analysis'): os.makedirs(rf'{folder}\analysis')
n_bin, Raw_bl = 60, []
Output_spike, Output_vm = {}, {}
Files = []

for sub_folder in sub_folders:
    drug = sub_folder.split('\\')[-1]
    print('\n'*2, '='*30, '\n', drug, '\n', '='*30, '\n'*2)
    files = glob.glob(rf'{sub_folder}\*.abf')
    Files.append([x.split('\\')[-1] for x in files])
    n = len(files)
    
    sub_folder_analysis = rf'{sub_folder}\analysis'
    if not os.path.exists(sub_folder_analysis): os.makedirs(sub_folder_analysis)
    
    matrix_spike, matrix_vm, raw_bl = np.zeros((n, n_bin)), np.zeros((n, n_bin)), []
    for f, file in enumerate(files):
        file_name = (file.split('\\')[-1]).split('.')[0]
        print(file_name, '\n'*2)
        
        abf = pyabf.ABF(file)
        raw = list(abf.sweepY)[:600 * abf.sampleRate]  
        y_ax = toto_filter(raw, sample_rate=abf.sampleRate, freq_low=1, freq_high=500)    
        
        gaussian(abf, 100)
        abf.setSweep(0)
        envlp = abf.sweepY[:600 * abf.sampleRate]  

        #Find frequencies
        indexes, _ = find_peaks(y_ax, height = 10)
        binary = np.zeros(len(y_ax))
        for ind in indexes: binary[ind] = 1
        splt = np.sum(np.split(binary, n_bin), axis=1)
        matrix_spike[f,:] = splt/max(splt)
        raw_bl.append(np.mean(splt[:6])/10)
        
        #Find Membrane potential
        matrix_vm[f,:] = np.nanmean(np.split(envlp, n_bin), axis=1)
        
        
        plt.figure(), plt.title(file_name)
        plt.subplot(511), plt.plot(raw[::2]), plt.title('Raw trace')
        plt.subplot(512), plt.plot(envlp[::2]), plt.title('Envlp')
        plt.subplot(513), plt.plot(y_ax[::2]), plt.title('Spike trace')
        plt.subplot(514), plt.plot(binary), plt.title('Binary')
        plt.subplot(515), plt.plot(matrix_spike[f,:]), plt.title('Normalized')
        plt.savefig(rf'{sub_folder_analysis}\{file_name}.pdf'), plt.close()
        
    Output_spike[drug] = matrix_spike
    Output_vm[drug] = matrix_vm
    Raw_bl.append(raw_bl)

    plt.figure(), plt.suptitle(drug)
    for i, (matrix, label) in enumerate(zip((matrix_spike, matrix_vm), ('AP', 'Vm'))):
        plt.subplot(4,1,(i*2)+1), plt.ylabel(f'Time course of {label}')
        mean, sem = np.nanmean(matrix, axis=0), stat.sem(matrix, axis=0)
        plt.fill_between(np.arange(len(mean)), mean-sem, mean+sem,
                          color='lightblue', alpha=0.25, zorder=1)
        plt.plot(mean, c='b', zorder=2), plt.xlabel('Time (10s bins)')
        plt.axvline(6, c='gold', lw = 2)
        plt.subplot(4, 1, (i*2)+2), sns.heatmap(matrix, cmap="coolwarm")
    plt.savefig(rf'{sub_folder_analysis}\Heat course and Time map.pdf'), plt.close()

    

plt.figure()
writer = pd.ExcelWriter(f'{folder}/analysis/all data spike.xlsx')

for d, (drug, files) in enumerate(zip(Output_spike.keys(), Files)):
    data = [(cell[:6], cell[8:14], cell[-6:]) for cell in Output_spike[drug]]
    data = np.asarray([[np.nanmean(x) for x in cell] for cell in data])
    
    Cell_id = [(i,)*3 for i in range(data.shape[0])]
    df = pd.DataFrame({'Cell_ID': [x for y in Cell_id for x in y],
                      'Time': ('Bl', drug, 'Wsh')*data.shape[0],
                      'Score':[x for y in data for x in y]})
    
    mean, sem = np.mean(data, axis=0), stat.sem(data, axis=0)
    x_ax = np.arange(d*3, d*3+3)
    
    plt.bar(x_ax, mean, yerr=sem, width=.9, capsize=3, label=drug)
    [plt.plot(x_ax, d, lw=.5) for d in data]
    plt.legend()
    
    df2 = pd.DataFrame(data)
    df2.index = files
    df2.to_excel(writer, sheet_name=drug) 

plt.savefig(rf'{folder}\analysis\Histo Hz.pdf'), plt.close()
writer.save()

plt.figure()
writer = pd.ExcelWriter(f'{folder}/analysis/all data vm.xlsx')

for d, (drug, files) in enumerate(zip(Output_vm.keys(), Files)):
    data = [(cell[:6], cell[8:14], cell[-6:]) for cell in Output_vm[drug]]
    data = np.asarray([[np.nanmean(x) for x in cell] for cell in data])
    
    Cell_id = [(i,)*3 for i in range(data.shape[0])]
    df = pd.DataFrame({'Cell_ID': [x for y in Cell_id for x in y],
                      'Time': ('Bl', drug, 'Wsh')*data.shape[0],
                      'Score':[x for y in data for x in y]})
    
    mean, sem = np.mean(data, axis=0), stat.sem(data, axis=0)
    x_ax = np.arange(d*3, d*3+3)
    
    plt.bar(x_ax, mean, yerr=sem, width=.9, capsize=3, label=drug)
    [plt.plot(x_ax, d, lw=.5) for d in data]
    plt.legend()
    
    df2 = pd.DataFrame(data)
    df2.index = files
    df2.to_excel(writer, sheet_name=drug) 

plt.savefig(rf'{folder}\analysis\Membrane voltage.pdf'), plt.close()
writer.save()




