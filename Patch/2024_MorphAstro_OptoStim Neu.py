# -*- coding: utf-8 -*-
"""
Created on Sat Apr 22 11:32:26 2023

@author: Angel.BAUDON
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Mar  3 10:59:14 2023

@author: Angel.BAUDON
"""

import matplotlib.pyplot as plt, numpy as np, glob, pandas as pd, os
import scipy.stats as stat, seaborn as sns
from scipy.signal import find_peaks
from scipy.optimize import curve_fit
from ToolKit.Tollbox_TG_shorten import Rawtrace, toto_filter


folder = r"C:\Angel.BAUDON\Biblio\Publi\MoprhAstro in prep\v2\Data\Figure 6\Rec Neu CeL OptoStim BLA\VClamp"
sub_folders = [x for x in glob.glob(rf'{folder}\*') if x.split('\\')[-1]!='analysis']
if not os.path.exists(rf'{folder}\analysis'): os.makedirs(rf'{folder}\analysis')
show_fig, save_fig = False, True

Params = ('CPSE freqeuncy (Hz)', 'CPSE amplitude (pA)', 'Amplitude(pA) TC', 'Failure rate TC', 'Rises TC',
          'Decays TC', 'Delta CPSE freqeuncy (Hz)', 'Delta CPSE amplitude (Hz)',
          'Failure rate', 'Rises', 'Decays', 'PPR')
Output = {x:{} for x in Params}

for sub_folder in sub_folders:
    files, drug = glob.glob(rf'{sub_folder}\*.wcp'), sub_folder.split('\\')[-1]
    print('\n'*2, '='*30, '\n', drug, '\n', '='*30, '\n'*2)
    
    sub_folder_analysis = rf'{sub_folder}\analysis'
    if not os.path.exists(sub_folder_analysis): os.makedirs(sub_folder_analysis)
    
    CPSE, dCPSE, Rises, Decays, Fails, PPR = np.zeros((len(files),3)), [], [], [], [], []
    CPSEa, dCPSEa = np.zeros((len(files),3)), []
    Amps, Fails_tc, Rises_tc, Decays_tc = [], [], [], []
    
    for f, file in enumerate(files):
        file_name = (file.split('\\')[-1]).split('.')[0]
        print(file_name, '\n'*2)
        
        raw = Rawtrace(file)
        y_ax, sampling = np.concatenate(raw.matrix), raw.sampling_rate
        rec_len = len(y_ax)/sampling 
        
        print('Wollen sie mit mir schlaffen? \n\n')
        trace = toto_filter(y_ax, sample_rate=sampling, freq_low=1, freq_high=500)[::10]
        sampling = int(sampling/10)

        indx = (20, 30, 70, 99.8, 109.8, 149.8, 179.6, 189.6, 230)
        indx = np.array([int(x*sampling) for x in indx])
        sub_indx = [np.linspace(start, stop, 200, dtype=int)
                    for (start, stop) in np.split(np.delete(indx, (2,5,8)), 3)]
                
        sd = np.std(trace[:indx[0]])
        cpse, _ = find_peaks(-trace, height=3*sd, prominence=sd)
        
        x_ax = np.linspace(0, rec_len, len(trace))
        # plt.figure(), plt.title(sd)
        # plt.plot(x_ax, trace)
        # plt.plot(cpse/sampling, [trace[i] for i in cpse], 'xr')
        # [plt.axvline(i/sampling, color='gold') for i in indx]

        
        bin_cpse, bin_amp = np.zeros(len(trace)), np.zeros(len(trace))
        bin_amp[:] = np.nan
        
        for i in cpse:
            bin_cpse[i] += 1
            bin_amp[i] = trace[i]

        splt_cpse, splt_amp = np.split(bin_cpse, indx)[:-1], np.split(bin_amp, indx)[:-1]
        
        #Select the 5 firsts seconds of the 10 last seconds of each run        
        med_splt_cpse = [run[-(10*sampling):-(5*sampling)] for run in splt_cpse]
        dcpse = np.nanmean([np.nanmean(x)*(5*sampling) for x in med_splt_cpse[1::3]])

        #Select the 5 firsts seconds of the 10 last seconds of each run
        med_splt_amp = [run[-(10*sampling):-(5*sampling)] for run in splt_amp]
        damp = np.nanmean([np.nanmean(x) for x in med_splt_amp[1::3]])
        
        print('Tatinos il sent le caca \n\n')
            
        CPSE[f,:] = [np.nanmean(med_splt_cpse[i::3])*(5*sampling) for i in (0,1,2)]
        CPSEa[f,:] = [np.nanmean(med_splt_amp[i::3]) for i in (0,1,2)]

        def expo(x, a, b, c): return a*np.exp(-b*x)+c #a=size, b=angle, c=intercept
        def reglin(x, a, b): return a*x + b
        
        
        print('Bouchouuuuuuuuuuuuuuuuuuuuuuuuuuuurhein')
        
        cell_fails, cell_amps = np.zeros((3, 199)), np.zeros((3, 199))
        cell_rise, cell_decay = np.zeros((3, 199)), np.zeros((3, 199))
        cell_rise[:], cell_decay[:] = np.nan, np.nan
        
        for si, sub_i in enumerate(sub_indx):
            stims = np.split(trace, sub_i)[1:-1]
            for s, stim in enumerate(stims):
                start, stop = sub_i[s], sub_i[s+1]
                stim_cpse = cpse[np.where(np.logical_and(cpse>start, cpse<stop))]
                
                if len(stim_cpse):
                    indx_stim_cpse = stim_cpse[0] - start -10
                    stim = stim[10:]-stim[0]
                    rise, decay = stim[:indx_stim_cpse], stim[indx_stim_cpse:]
                    amp, fail = stim[indx_stim_cpse], 0
                    
                    # plt.figure(), plt.title(indx_stim_cpse)
                    # plt.subplot(2,1,1), plt.plot(stim), plt.title(indx_stim_cpse)
                    
                    for k, kin in enumerate((rise, decay)):
                        x = np.arange(len(kin))

                        try:
                            # fig = plt.subplot(2, 2, 3+k)
                            # fig.plot(x, kin)
                            func, p0 = (reglin, (-.1, 5)) if k else (expo, (-5, 0.01, 5))
                            popt, _ = curve_fit(func, x, kin, p0=p0)
                            fit = func(x, *popt)-min(func(x, *popt))
                            low, high = max(fit)*.1, max(fit)*.9
                            tau = len(np.where(np.logical_and(fit>=low, fit<=high))[0])/sampling
                            # fig.plot(x, func(x, *popt)), fig.set_title(tau)
                            if not k: cell_rise[si,s] = tau
                            else: cell_decay[si,s] = tau
                        except: continue

                else: amp, fail = np.nan, 1
                
                cell_amps[si,s], cell_fails[si,s] = amp, fail

        for all_, cell_ in zip((Amps, Fails_tc, Rises_tc, Decays_tc),
                                (cell_amps, cell_fails, cell_rise, cell_decay)):
            all_.append(np.nanmean(cell_, axis=0))

        plt.figure(), plt.title(f'{file_name} \n {np.std(trace)}')
        plt.plot(x_ax, trace, c='b', lw=.5, label='Raw')
        [plt.axvline(i/sampling, c='r') for i in indx]
        [plt.axvline(i/sampling, c='lightblue', zorder=0)
          for i in np.concatenate(sub_indx)]
        plt.plot(cpse/sampling, [trace[i] for i in cpse], 'xr')
        print("J'ai envie de manger du moumous mais il est 9h du matin...")
        plt.legend(loc='upper right')
        if save_fig: plt.savefig(rf'{sub_folder_analysis}\{file_name}.pdf')
        if not show_fig: plt.close()
        

        
        PPR.append(np.nanmean([x[1]/x[0]for x in cell_amps]))
        Rises.append(np.nanmean(cell_rise)), Decays.append(np.nanmean(cell_decay))
        Fails.append(np.nanmean(cell_fails)), dCPSE.append(dcpse), dCPSEa.append(damp)
        
    for k, data in zip(Output.keys(), (CPSE, CPSEa, Amps, Fails_tc, Rises_tc, Decays_tc,
                                        dCPSE, dCPSEa, Fails, Rises, Decays, PPR)):
        Output[k][drug] = data




'''           Output section         '''
for k in Output.keys():
    data = Output[k]
    plt.figure(), plt.title(k), plt.ylabel(k)
    TRT = [x.split('\\')[-1] for x in sub_folders]
    Trt = [a for b in [(T,)*len(Output[k][O]) for T, O in zip(TRT, Output[k])] for a in b]

    if k in Params[:2]:
        mean = [np.nanmean(data[x], axis=0) for x in data]
        sem = [stat.sem(data[x], axis=0, nan_policy='omit') for x in data]
        x_ax = np.arange(len(sub_folders)*3)
        plt.bar(x_ax, np.hstack(mean), yerr=np.hstack(sem), width=.9, capsize=3)
        for k_data, x in zip(data, np.split(x_ax, len(data))):
            for cell in data[k_data]: plt.plot(x, cell, lw=.5, c='k')
        
        
        writer = pd.ExcelWriter(f'{folder}/analysis/{k}.xlsx')
        for l in data.keys():
            sub_data = data[l]
            df = pd.DataFrame(sub_data)
            df.to_excel(writer, sheet_name = l)
        writer.save()
        

    elif k in Params[2:6]:
        means = [np.nanmean(data[x], axis=0) for x in data]
        sems = [stat.sem(data[x], axis=0, nan_policy='omit') for x in data]
        x_ax = np.arange(199)
        
        from scipy.signal import savgol_filter
        means = [savgol_filter(mean, 15, 1) for mean in means]
        sems = [savgol_filter(mean, 15, 1) for mean in sems]
        
        
        for mean, sem, label in zip(means, sems, TRT):
            plt.plot(x_ax, mean, label=label)
            plt.fill_between(x_ax, mean-sem, mean+sem, alpha=0.5, zorder=1)
            plt.legend()
            
        writer = pd.ExcelWriter(f'{folder}/analysis/{k}.xlsx')
        for l in data.keys():
            sub_data = data[l]
            df = pd.DataFrame(sub_data)
            df.to_excel(writer, sheet_name = l)
        writer.save()

    else:
        df = pd.DataFrame({'Treatment': Trt,
                            k: [a for b in [Output[k][O] for O in Output[k]] for a in b]})
        sns.barplot(data=df, x='Treatment', y=k, color='lightgreen', capsize=.2)
        sns.swarmplot(x='Treatment', y=k, data=df, size=5, color='w',
                      edgecolor='k', linewidth=1)
        # if not show_fig: plt.close()

        #df_stat = IntraGrpStat(df)
        
        writer = pd.ExcelWriter(f'{folder}/analysis/{k}.xlsx')
        df.to_excel(writer, sheet_name = k)
        #df_stat.to_excel(writer, sheet_name = 'Stats')
        writer.save()
    
    # if save_fig: plt.savefig(rf'{folder}\analysis\{k}.pdf')
    # if not show_fig: plt.close()

writer = pd.ExcelWriter(f'{folder}/analysis/Amp TC.xlsx')
for k in Output['Amplitude(pA) TC'].keys():
    bjr = Output['Amplitude(pA) TC'][k]
    pd.DataFrame(bjr).to_excel(writer, sheet_name = k)
writer.save()


