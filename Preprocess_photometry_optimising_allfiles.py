##
# Import functions
import os
from pathlib import Path
import numpy as np
from scipy.optimize import curve_fit
from scipy.signal import savgol_filter, butter, filtfilt
from sklearn.impute import SimpleImputer

from StimTimePoints import *

"""
Code optimising 2023/06
"""

##
#Parameter to analyze
animal_folder = "FP114"
File_Name_list = ("20250914_FP114_FC3_day2_4", "20250913_FP114_FC3_day1_4")
Selected_Ch = ['Ch2']


##
for File_Name in File_Name_list:
    # Import data file
    data_folder = Path(f"C:/Users/Kai-Yi.WANG/PycharmProjects/Photometry_analysis/NPZ_file/{animal_folder}/")
    Import_file = data_folder / File_Name
    load_data = np.load(f'{Import_file}.npz', 'r', allow_pickle=True)

    Time, C1, C2, C3, C4 = load_data['Time'], load_data['C1'], load_data['C2'], load_data['C3'], load_data['C4']

    # match the channels
    NamedCH = {"Ch1": C1,
               "Ch2": C2,
               "Ch3": C3,
               "Ch4": C4}
    Ch_toCheck = []
    for CHH in Selected_Ch:
        ext_ch_signals = NamedCH[CHH]
        Ch_toCheck.append(ext_ch_signals)

    # check the numbers of NaN in array
    def CheckNaN(Ch_toCheck):
        NumOfNan = np.isnan(Ch_toCheck).sum()
        print(NumOfNan)

        if NumOfNan != 0:
            imputeNaN = SimpleImputer(missing_values=np.nan, strategy='most_frequent')
            replaceNaN = imputeNaN.fit_transform(Ch_toCheck.reshape(-1, 1))
            new_dataframe = replaceNaN.flatten()
            return new_dataframe
        else:
            return Ch_toCheck

    ChTobeProcessed = []
    for ind_ch in range(len(Selected_Ch)):
        after_repl_NaN = CheckNaN(Ch_toCheck[ind_ch])
        ChTobeProcessed.append(after_repl_NaN)

    # Apply a Savitzky-Golay filter
    win_size = round(len(Time)/1000)  #!!ToDo: test the effect window size on signals
                                      # For now: 1000 for downsample 10 dataset; 5000 for downsample 100 dataset
    if win_size % 2 == 0:
        win_size += 1

    # win_size = 201
    # # check if window is greater than signal shape
    # if win_size > len(ChTobeProcessed[0]):
    # 	win_size = ((win_size+1)/2)+1

    filtered_selected = []
    for ch_to_filter in ChTobeProcessed:
        sg_filter = savgol_filter(ch_to_filter, window_length=win_size, polyorder=2)
        filtered_selected.append(sg_filter)

    namedfilter = {}
    for FCh_into_dic in range(len(Selected_Ch)):
        namedfilter[f'filtered_{Selected_Ch[FCh_into_dic]}'] = filtered_selected[FCh_into_dic]


    # Create the artificial control channel from filtered signal and convert the signal into DFF and z-score

    class CreateCtrl_and_NormToS:
        """
        **Create control channels from filtered signals and normalized the signals to control channels**
        rec_T = Time
        CH = namedfilter (dictionary)
        CHNAME = channel for fitting, key from {namedfilter}
        """

        def __init__(self, rec_T, CH, CHNAME):
            self.rec_T = rec_T
            self.CH = CH
            self.CHNAME = CHNAME

        # def Create_control_from_signal(self):
        #     ChForFit = self.CH[self.CHNAME] # targeted filtered channel
        #
        #     def fit_expo(x, a, b, c):
        #         return a * np.exp(-b * x) + c
        #
        #     try:
        #         popt, pcov = curve_fit(fit_expo, self.rec_T, ChForFit)
        #     except RuntimeError:
        #         print('RuntimeError: Optimal parameters not found! Try Polynomial fit')
        #         Ch_coef = np.polyfit(self.rec_T, ChForFit, deg=2)  # !!ATTENTION to the degree value!!
        #         Ch_polyfit = np.polyval(Ch_coef, self.rec_T)
        #         control_arti = Ch_polyfit
        #         Fit_type = 'Polynomial fit'
        #
        #     else:
        #         print('Works with mono_fit!')
        #         print(popt)
        #         control_mono = fit_expo(self.rec_T, *popt)
        #         control_arti = control_mono
        #         Fit_type = 'Monoexponential fit'
        #
        #     return ChForFit, control_arti, Fit_type
        def Create_control_from_signal(self):
            ChForFit = self.CH[self.CHNAME] # targeted filtered channel


            Ch_coef = np.polyfit(self.rec_T, ChForFit, deg=3)  # !!ATTENTION to the degree value!!
            Ch_polyfit = np.polyval(Ch_coef, self.rec_T)
            control_arti = Ch_polyfit
            Fit_type = 'Polynomial fit'

            return ChForFit, control_arti, Fit_type

        def Convert_into_DFF(self):
            Target_CH, control_CH, Fit_type= self.Create_control_from_signal()

            dF = np.subtract(Target_CH, control_CH)
            dFF = np.divide(dF, control_CH)
            nor_dFF = dFF * 100

            return nor_dFF, Target_CH, control_CH, Fit_type

        def Calculate_zScore(self):
            # DFF_toZ = self.Convert_into_DFF()[0]
            # DFF_toZ = nor_dFF
            dF = np.subtract(Target_CH, control_CH)
            Sub_dff_mean = np.subtract(dF, np.mean(dF))
            zscore = Sub_dff_mean / np.std(dF)
            # zscore = np.divide(Sub_dff_mean, np.std(dF))

            return zscore


    #Compute the DFF and save parameters into npz

    channelmatch = {'Ch1': "filtered_Ch1", 'Ch2': "filtered_Ch2", 'Ch3': "filtered_Ch3", 'Ch4': "filtered_Ch4"}
    output_dir = f"C:/Users/Kai-Yi.WANG/PycharmProjects/Photometry_analysis/NPZ_file/{animal_folder}/"

    for Search_ch in list(channelmatch.keys()):
        for i in range(len(Selected_Ch)):
            if Search_ch == Selected_Ch[i]:
                PrePro = CreateCtrl_and_NormToS(Time, namedfilter, channelmatch[Selected_Ch[i]])
                nor_dFF, Target_CH, control_CH, Fit_type = PrePro.Convert_into_DFF()
                # zscore = PrePro.Calculate_zScore()

                np.savez(os.path.join(output_dir, f'{File_Name[0:-2]}_{Selected_Ch[i]}.npz'),
                         Time = Time, DFF = nor_dFF, Target_CH = Target_CH, control_CH = control_CH, Fit_type = Fit_type)

##
# For FP6-8 and FP1
# def Create_control_from_signal(match_T_Ch):
#     fCh = {key: val for key, val in match_T_Ch.items()} #if key < 640} #!!!only FP6-8 and FP1 need!
#     fCh_k = np.array(list(fCh.keys()))
#     fCh_v = np.array(list(fCh.values()))
#
#     def fit_expo(x, a, b, c):
#         return a * np.exp(-b * x) + c
#
#     try:
#         popt, pcov = curve_fit(fit_expo, fCh_k, fCh_v)
#     except RuntimeError:
#         print('RuntimeError: Optimal parameters not found! Try Polynomial fit')
#         C2_coef = np.polyfit(fCh_k, fCh_v, deg=2)  # !!ATTENTION to the degree value!!
#         C2_polyfit = np.polyval(C2_coef, fCh_k)
#         control_arti = C2_polyfit
#         Fit_type = 'Polynomial fit'
#
#     else:
#         print('Works with mono_fit!')
#         print(popt)
#         control_mono = fit_expo(fCh_k, *popt)
#         control_arti = control_mono
#         Fit_type = 'Monoexponential fit'
#
#     return fCh_k, fCh_v, control_arti, Fit_type


# print(Create_control_from_signal(match_T_C2))
#
# fC1_k, fC1_v, ArtiCtrl_C1, C1_fit_type = Create_control_from_signal(match_T_C1)
# fC2_k, fC2_v, ArtiCtrl_C2, C2_fit_type = Create_control_from_signal(match_T_C2)
# fC3_k, fC3_v, ArtiCtrl_C3, C3_fit_type = Create_control_from_signal(match_T_C3)
# fC4_k, fC4_v, ArtiCtrl_C4, C4_fit_type = Create_control_from_signal(match_T_C4)




