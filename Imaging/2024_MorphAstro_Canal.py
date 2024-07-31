# -*- coding: utf-8 -*-
"""
Created on Mon Apr 19 08:48:04 2021

@author: Angel.BAUDON

To do on this script:
    - Add smth that exclude bursts begining before drug application
    - More precise peak foots
"""
import numpy as np, pandas as pd, os, winsound, matplotlib.pyplot as plt
from ToolKit import Toolz1, Analysor, PlotFactory

#Locate the experiment folder
folder_exp = r"C:\Users\Angel.BAUDON\Desktop\ImCa Clem\PFC\dOVT + TGOT"

#Define parameters
Select_cell, remove_inactive = False, False
drug = folder_exp.split('\\')[-1]

#Create folder to stock the final analysis
M = Toolz1.Manager(folder_exp, Select_cell)
Analysis, General, Response, Writer = M.FileMaker(drug, ('Analysis', 'General', 'Response'))

#Define parameters of the analysis & create empty lists
drug, drug_time, sampling_Hz, exp_duration = drug, 270, 2, 900
drug_time, rec_duration = drug_time*sampling_Hz, exp_duration*sampling_Hz

#Create the Output data frame
Parameters = ('Slice', 'Cell n°', 'AUC', 'mHz', 'Height', 'HMFW', 'Rise', 'Decay',
              'Area', 'Active', 'Resp', 'Binary', 'AUC time course',
              'Resp border', 'Resp duration', 'Resp delay')
Output = pd.DataFrame(columns = Parameters)
Raw = []

#Combine the root and the file to obtain the file location
folder_exp_list = [x for x in os.listdir(folder_exp) if x != 'Analysis']
for file_day in folder_exp_list:
    folder_day = f'{folder_exp}\{file_day}'
    
    #Create folder to stock the slice analysis
    M = Toolz1.Manager(folder_day, Select_cell)
    analysis, analysis_cell, analysis_raw, writer = M.FileMaker(drug, ('Analysis', 'Final traces',
                                                                        'Raw traces'))

    #Find Excel file and exctract raw data
    file_name = [x for x in os.listdir(folder_day) if x.split('.')[-1]=='xlsx'][0]
    file = (f'{folder_day}\{file_name}')
    print('='*30, '\n'*2, file_name, '\n'*2, '='*30)
    try: data, n_cell, writer = M.DataExtractor(file, Select_cell, writer, 'right')
    except TypeError: continue
    
    #Generate a DataFrame to stock measured parameters
    Paramz = (*Parameters[1:], 'Raw', 'RegLin', 'dFF0', 'Indexes', 'Bursts', 'Peak borders')
    OutP = pd.DataFrame(columns = Paramz)
    
    '''                      Begining of the analysis                   '''
    '''                          Select one cell                        '''

    for c in range(n_cell):
        A = Analysor.Analysor(data, c, rec_duration, sampling_Hz, drug_time)
        
        #Isolate a cell
        raw, rl, cell, coef = A.CellReader(RegLin=True, Envlp=(1,5))
        Raw.append(cell/max(cell))
        
        #Find events
        indexes, heights = A.PeakDetector(analysis_raw, show_fig=False)
        
        #Create a binary list of events
        cell_binary = A.Binary()
        
        #Find bursts
        bursts_cell = A.FindBurst()
        
        #Peak Properties Analysis
        Mimir, Hz, AUC, AUC_cell, borders = A.PeakAnalysor()
        
        #Label the cell
        active = True if 1 in cell_binary else False
        resp, resp_param = A.Response_finder()

        #Stock section
        #General Stock
        Out = (file_name, c, AUC, Hz, *Mimir, active, resp, cell_binary, AUC_cell, *resp_param)
        Output = Output.append({x:y for x, y in zip(Parameters, Out)}, ignore_index = True)
        
        #Stock by slice
        O = (*Out[1:], raw, rl, cell, indexes, bursts_cell, borders)
        OutP = OutP.append({x:y for x, y in zip(Paramz, O)}, ignore_index = True)

    #Fill and save excel file
    for key in OutP.columns[1:]:
        if key in ('Raw', 'RegLin', 'dFF0'):
            df = pd.DataFrame({x:y for x,y in zip(OutP.index, OutP[key])})
            df.to_excel(writer, sheet_name = key)
        else: OutP[key].to_excel(writer, sheet_name = key)
    writer.save()
    
    #Plot section
    PlotFactory.PlotMaker(OutP, drug, drug_time, file_name, analysis, analysis_cell, analysis_raw)
   

'''                     Final figures                                       '''

Output_plot = Output[Output['Active'] == True] if remove_inactive else Output
FP = PlotFactory.FinalPlots(Output_plot, drug, General, rec_duration, drug_time)

#Time courses
FP.TimeCoursesMaker(exp_duration, bin_len = 10)

#Histos
Output2 = FP.HistoMaker(Response, General, violin = False)

#Mosaic
Output2 = FP.MosaMaker(General, violin = False)

#Radar
normalizer = [1]*22
FP.RadarMaker(normalizer)

#Raster
FP.RasterMaker()

#Time course
import scipy.stats as stat
mean, sem = np.nanmean(Raw, axis=0), stat.sem(Raw, axis=0)
x_ax = np.arange(len(mean))
plt.figure(), plt.plot(x_ax, mean)
plt.fill_between(x_ax, mean-sem, mean+sem, alpha=0.5, zorder=1)


'''                      Output to excel                                    '''

#Create an excel file & stock the final Values
writer_exp = pd.ExcelWriter(f'{folder_exp}/Analysis/{drug}.xlsx')
Output.to_excel(writer_exp, sheet_name = 'All data')
Output2.to_excel(writer_exp, sheet_name = 'Stats')

writer = pd.ExcelWriter(f'{folder_exp}/Analysis/{drug} all data.xlsx')
for key in list(Output2.keys())[:-3]:
    data = list(Output2[key])
    pd.DataFrame({'Baseline':data[0], key:data[1]}).to_excel(writer, sheet_name=key)

writer_ratio = pd.ExcelWriter(f'{folder_exp}/Analysis/{drug} all ratios.xlsx')

for key in list(Output2.keys())[:-3]:
    ratio = list(Output2[key][2])
    pd.DataFrame({'Ratio':ratio}).to_excel(writer_ratio, sheet_name=key)
writer_exp.save(), writer.save(), writer_ratio.save()
