# -*- coding: utf-8 -*-
"""
Created on Mon Mar 20 10:22:20 2023

@author: Angel.BAUDON
"""

import matplotlib.pyplot as plt, numpy as np, glob, pandas as pd, os, scipy.stats as stat
from scipy.optimize import curve_fit
from scipy.signal import savgol_filter
from ToolKit.Tollbox_TG_shorten import Rawtrace


folder = r"C:\Angel.BAUDON\Exp\Neuro\Patch\Data\2023_03_PatchAstro\Optostim Astro tout bien rangé\Stim V clamp"
sub_folders = [x for x in glob.glob( rf'{folder}\*') if x.split('\\')[-1] != 'analysis']
if not os.path.exists(rf'{folder}\analysis'): os.makedirs(rf'{folder}\analysis')
show_fig, save_fig = False, True
Params = ('Raw', 'Current (pA)', 'Delta current', 'Rise (s)', 'Decay (s)', 'Substracted current')
Output = {x: {} for x in Params}
Time_course = {}

for sub_folder in sub_folders:
    drug = sub_folder.split('\\')[-1]
    print('\n'*2, '='*30, '\n', drug, '\n', '='*30, '\n')
    cells = [x for x in glob.glob(rf'{sub_folder}\*') if x != 'analysis']
    Output_drug = {x: {} for x in Params}
    time_course_drug = []
    
    for cell in cells:
        files, cell_name = glob.glob(rf'{cell}\*.wcp'), cell.split('\\') [-1]
        print('\n', ':'*30, '\n', cell_name, '\n', ':'*30, '\n'*2)
    
        cell_analysis = rf'{cell}\analysis'
        if not os.path.exists(cell_analysis): os.makedirs(cell_analysis)
    
        Raw, Vm, dVm, Rises, Decays = [[] for i in range(5)]
        for f, file in enumerate(files):
            file_name = (file.split('\\')[-1]).split('.')[0]
            print(file_name, '\n'*2)
    
            raw = Rawtrace(file)
            y_ax = raw.matrix[0][::100]
            sampling = raw.sampling_rate/100
            rec_len = y_ax.shape[0]/sampling
            trace = savgol_filter(y_ax, 101, 1)
            indx = np.array([int(x*sampling) for x in (20, 30)])
            
            splt_vm = np.split(np.asarray(trace), indx)
            med_splt_vm = [min(x) for x in splt_vm]

            def expo(x, a, b, c): return a*np.exp(-b*x)+c  # a=size, b=angle, c=intercept
    
            lambdas = []
            for s, splt in enumerate(splt_vm[1:]):
                try:
                    x = np.linspace(0, len(splt), len(splt))
                    if folder.split('\\')[-1] == 'VClamp':
                        (size, loc, mark) = (-50, -1, .632) if s else (50, 0, .368)
                    else: (size, loc, mark) = (1, 0, .632) if s else (-1, -1, .368)
                    popt, _ = curve_fit(expo, x, splt, p0=[size, 2e-04, 0])
                    fit = expo(x, *popt)-min(expo(x, *popt))
                    tau = len(fit[np.where(np.logical_and(fit>.1*max(fit), fit<.9*max(fit)))])
                    lambdas.append(tau/sampling)
                    # plt.figure(), plt.plot(x, splt), plt.plot(x, expo(x, *popt)), plt.title(tau)
                except RuntimeError:
                    lambdas.append(np.nan)
                    continue
            dvm = ((med_splt_vm[1])-(med_splt_vm[0]))
            Raw.append(trace), Rises.append(lambdas[0]), Decays.append(lambdas[1])
            Vm.append(med_splt_vm), dVm.append(dvm)
            
            plt.figure(), plt.title(file_name)
            x_ax = np.linspace(0, rec_len, len(y_ax))
            plt.plot(x_ax, y_ax, label='Raw'), plt.plot(x_ax, trace, label='Fltr')
            plt.axvspan(*indx/sampling, color='lightblue'), plt.legend(loc='upper right')
            if save_fig: plt.savefig(rf'{cell_analysis}\{file_name}.pdf')
            if not show_fig: plt.close()

        for k, data in zip(Params, (Raw, Vm, dVm, Rises, Decays, dVm[0]-dVm[1])):
            Output_drug[k][cell_name] = data
            
            
        plt.figure(), plt.title(f'{cell_name}, {dVm}, {dVm[1]-dVm[0]}')
        x_ax = np.linspace(0, rec_len, len(y_ax))
        raw_fig = [x-np.mean(x[:indx[0]]) for x in Raw]
        for raw, label in zip(raw_fig, ('Ctrl', 'Ba2+')): plt.plot(x_ax, raw, label=label)
        plt.plot(x_ax, raw_fig[0]-raw_fig[1], label='Substraction')
        plt.axvspan(*indx/sampling, color='lightblue'), plt.legend(loc='upper right')
        if save_fig: plt.savefig(rf'{cell_analysis}\{file_name} substraction.pdf')
        if not show_fig: plt.close()
        Output_drug[k][cell_name] = data
        time_course_drug.append(raw_fig[0]-raw_fig[1])
            
    for k in Params: Output[k][drug] = Output_drug[k]
    Time_course[drug] = time_course_drug
    
    

writer = pd.ExcelWriter(f'{folder}/analysis/all_data.xlsx')
for k in Params[1:-3]:
    data, out = Output[k], np.zeros((100,100))
    out[:] = np.nan
    
    plt.figure()
    for d, drug in enumerate(data):
        matrix = np.asarray([data[drug][cell] for cell in data[drug]])        
        mean, sem = np.nanmean(matrix, axis=0), stat.sem(matrix, axis=0, nan_policy='omit')
        fig = plt.subplot(1, len(sub_folders), d+1)
        fig.set_title(drug), fig.set_ylabel(k)

        if k == 'Current (pA)':
            fig.bar(np.arange(6), np.hstack(mean), yerr=np.hstack(sem), width=.9, capsize=3)
            for a, x_bar in enumerate(((0,1,2), (3,4,5))):
                for m in matrix:
                    fig.plot(x_bar, m[a], lw=.5, c='k', marker='o', mfc='w', mec='k')
        else:
            out[:matrix.shape[0],3*d : matrix.shape[1]+3*d] = matrix
            pd.DataFrame(out).to_excel(writer, sheet_name = k)

            fig.bar(np.arange(2), mean, yerr=sem, width=.9, capsize=3)
            for m in matrix:
                fig.plot(np.arange(2), m, lw=.5, c='k', marker='o', mfc='w', mec='k')
                
    if save_fig: plt.savefig(rf'{folder}\analysis\{k}.pdf')
    # if not show_fig: plt.close()
    
plt.figure(), plt.title(Params[-1]), plt.ylabel(Params[-1])
data = Output[Params[-1]]
Ik = [[data[y][x] for x in data[y]]for y in data]
mean, sem = [np.nanmean(x) for x in Ik], [stat.sem(x) for x in Ik]
plt.bar(np.arange(4), mean, yerr=sem, width=.9, capsize=3)
for i, ik in enumerate(Ik):
    [plt.plot(i, x, lw=.5, c='k', marker='o', mfc='w', mec='k') for x in ik]
if save_fig: plt.savefig(rf'{folder}\analysis\Substracted current.pdf')
# if not show_fig: plt.close()


out_ik = np.zeros((20,4))
out[:] = np.nan

for i, ik in enumerate(Ik): out_ik[:len(ik),i] = ik

pd.DataFrame(out_ik).to_excel(writer, sheet_name = 'Ik')
writer.save()



plt.figure()
for k in list(Time_course.keys())[2:]: Time_course.pop(k)
means = [np.nanmean(np.asarray(Time_course[k]), axis=0) for k in Time_course.keys()]
sems = [stat.sem(np.asarray(Time_course[k]), axis=0) for k in Time_course.keys()]
time = np.linspace(0, 80, 11999)
for mean, sem, label in zip(means, sems, Time_course.keys()):
    plt.plot(time, mean, label = label)
    plt.fill_between(time, mean-sem, mean+sem, alpha=0.5, zorder=1)
    plt.legend()





