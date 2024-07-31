# -*- coding: utf-8 -*-
"""
Created on Mon Apr 22 14:43:26 2024

@author: Angel.BAUDON
"""

import matplotlib.pyplot as plt, numpy as np, glob, pandas as pd, os, scipy.stats as stat
from scipy.signal import find_peaks, savgol_filter
from ToolKit.Tollbox_TG_shorten import Rawtrace


folder = r"C:\Angel.BAUDON\Biblio\Publi\MoprhAstro in prep\v2\Data\Figure 6\Rec Neu CeL OptoStim BLA\Ctrl DNQX AP5"
neurons = [x for x in glob.glob(rf'{folder}\*') if x.split('\\')[-1]!='analysis']
if not os.path.exists(rf'{folder}\analysis'): os.makedirs(rf'{folder}\analysis')
show_fig, save_fig = False, True

Params = ('CPSE freqeuncy (Hz)', 'CPSE amplitude (pA)')
Output = {x:{} for x in Params}

for neuron in neurons:
    files = glob.glob(rf'{neuron}\*.wcp')
    print('\n'*2, '='*30, '\n', neuron, '\n', '='*30, '\n'*2)
    
    CPSE, CPSEa = np.zeros((len(files),3)), np.zeros((len(files),3))
    
    for f, file in enumerate(files):
        file_name = (file.split('\\')[-1]).split('.')[0]
        print(file_name, '\n'*2)
        
        raw = Rawtrace(file)
        y_ax, sampling = np.concatenate(raw.matrix), raw.sampling_rate
        rec_len = len(y_ax)/sampling 
        trace = savgol_filter(y_ax, 101, 1)[::10]
        sampling = int(sampling/10)
        envlp = savgol_filter(trace, 5001, 1)
        trace = trace - envlp
        
        indx = (20, 30, 70, 99.8, 109.8, 149.8, 179.6, 189.6, 230)
        indx = np.array([int(x*sampling) for x in indx])

                
        sd = np.std(trace[:indx[0]])
        ppse, _ = find_peaks(trace-envlp, height=5*sd, prominence=5*sd)        
        bin_ppse, bin_amp = np.zeros(len(trace)), np.zeros(len(trace))
        bin_amp[:] = np.nan
        
        for i in ppse:
            bin_ppse[i] += 1
            bin_amp[i] = trace[i]
         
        splt_ppse, splt_amp = np.split(bin_ppse, indx)[:-1], np.split(bin_amp, indx)[:-1]
        
        #Select the 5 firsts seconds of the 10 last seconds of each run        
        med_splt_ppse = [run[-(10*sampling):-(5*sampling)] for run in splt_ppse]
        med_splt_amp = [run[-(10*sampling):-(5*sampling)] for run in splt_amp]
        
        CPSE[f,:] = [np.nanmean(med_splt_ppse[i::3])*(5*sampling) for i in (0,1,2)]
        CPSEa[f,:] = [np.nanmean(med_splt_amp[i::3]) for i in (0,1,2)]

        x_ax = np.linspace(0, rec_len, len(trace))
        plt.figure(), plt.title(f'{file_name} \n {np.std(trace[:indx[0]])}')
        plt.plot(x_ax, trace, c='b', lw=.5, label='Raw')
        [plt.axvline(i/sampling, c='r') for i in indx]
        plt.plot(ppse/sampling, [trace[i] for i in ppse], 'xr')
        plt.legend(loc='upper right')
        if save_fig: plt.savefig(rf'{folder}\analysis\{file_name}.pdf')
        if not show_fig: plt.close()


    for k, data in zip(Output.keys(), (CPSE, CPSEa)):
        Output[k][neuron[-8:]] = data




'''           Output section         '''
for k in Output.keys():
    data = Output[k]
    plt.figure(), plt.title(k), plt.ylabel(k)
    
    before = [data[k][0][1] for k in data.keys()]
    after = [data[k][1][1] for k in data.keys()]
    
    means = (np.nanmean(before), np.nanmean(after))
    sems = (stat.sem(before), stat.sem(after))
    plt.bar((0,1), means, yerr=sems, width=.9, capsize=3)
    
    for b, a in zip(before, after): plt.plot((0, 1), (b, a), lw=.5, c='k')
        
        
    writer = pd.ExcelWriter(f'{folder}/analysis/{k}.xlsx')
    df = pd.DataFrame({'Before': before, 'After': after})
    df.to_excel(writer)
    writer.save()
    

