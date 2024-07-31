# -*- coding: utf-8 -*-
"""
Created on Tue Feb  9 10:10:10 2021

@author: Angel.BAUDON
"""

import os
from xlrd import XLRDError
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter



class Manager():
    
    def __init__(self, folder_exp, select):
        self.folder_exp = folder_exp
        self.S = select


    def FileMaker(self, drug, names):
        F1 = f'{self.folder_exp}\{names[0]}'
        F2, F3 = f'{F1}\{names[1]}', f'{F1}\{names[2]}'
        [os.makedirs(path) for path in [F1, F2, F3] if not os.path.exists(path)]
        
        try:
            w = f'{F1}\{[x for x in os.listdir(F1) if x.split(".")[-1] == "xlsx"][0]}'
            wr = w if os.path.exists(w) else f'{F1}/{drug}.xlsx'
        except IndexError:
            wr = f'{F1}/{drug}.xlsx'
        writer = pd.ExcelWriter(wr) if self.S else wr
        return F1, F2, F3, writer


    
    def DataExtractor(self, file, Select_cell, writer_slice, screen_position, K=False):
        print('non, toi vas niquer tes morts')
        if not Select_cell:
            print(writer_slice)
            data = (pd.read_excel(writer_slice, sheet_name='Raw',
                                  header=None)).to_numpy()
            data, ks = data[1:,1:].T, 'y'
            writer_slice = pd.ExcelWriter(writer_slice)
            
        if Select_cell:
            OGB, SR = [(pd.read_excel(file, sheet_name=i, header=None)).to_numpy() for i in [1,0]]
            if not SR.size: SR = np.zeros((OGB.shape[0], OGB.shape[1]))
            
            r, c, data = (np.ma.size(OGB, 1)+1)//2, 2, []
            if r == 1: r, c = 2, 1
            
            if not K:
                mosa = plt.figure()
                mosa.subplots_adjust(hspace = 1, wspace = 0.5)
                for i, (ogb, sr) in enumerate(zip(OGB.T, SR.T)):
                    fig_plot = plt.subplot(r, c, (i + 1))
                    [plt.plot(x, c = y) for x, y in zip([ogb, sr], ['g', 'r'])]
                    fig_plot.xaxis.set_visible(False)
                
                move = 2000 if screen_position == 'right' else -1000
                mosa.canvas.manager.window.move(move, 300)
                mng = plt.get_current_fig_manager()
                mng.window.showMaximized(), plt.pause(.5)
                ks = input('Keep this slice ? (y/n)')
                while ks not in ('y', 'n', 'yy'):
                    ks = input('Answer not understood ! (y/n)')
                plt.close()
                
            if K: ks = 'yy'
            
            if ks == 'yy':
                for i, (O, S) in enumerate(zip(OGB.T, SR.T)):
                    data.append(O-savgol_filter(S, 51, 1))

            elif ks == 'y':
                for i, (O, S) in enumerate(zip(OGB.T, SR.T)):
                    trend = savgol_filter(S, 51, 1)
                    sub = O-trend
                    
                    fig = plt.figure()
                    [plt.plot(x, c=y) for x, y in zip([O,S,trend,sub], ['g', 'r', 'k', 'purple'])]
                    fig.canvas.manager.window.move(move, 300)
                    mng = plt.get_current_fig_manager()
                    mng.window.showMaximized(), plt.pause(.5)
                    
                    k = input('Keep this cell ? (y/n)')
                    while k not in ['y', 'n']: k = input('Answer not understood ! (y/n)')
                    # plt.savefig(rf'C:\Users\Angel.BAUDON\Desktop\Nouveau dossier/Cell {i}.pdf')
                    plt.close()
                    if k == 'y': data.append(sub)
        
        #Define the duration of the recording and the number of cells
        try:
            data = np.asarray(data).T
            n_cell = np.ma.size(data, 1)
        except IndexError: pass
        
        if ks != 'n' and len(data): return data, n_cell, writer_slice
        else: return
    
    
def line_plot(x, y, pval):
    plt.plot((x[0]+.1, x[1]-.1), (y,)*2, c='k', lw=1)
    plt.text(np.mean(x), y, Q(pval), size=10, weight='bold')


def Q(pval):
    if pval<=0.001: star = '***'
    elif 0.001<pval<=0.01: star = '**'
    elif 0.01<pval<=0.05: star = '*'
    else: star = 'ns'
    return star
