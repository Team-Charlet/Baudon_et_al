##
# Import functions
from pathlib import Path
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
from StimTimePoints import *
# import RFM_STIM
import csv

"""
Code optimising 2023/08
TRY to normalize baseline dF/F to 0, calculate z Score and area under curve
"""
##
para_FC1 = ['202303_Photometry+FC', 'Day2_US']
para_FC2 = ['202306_Photometry+FC2_reexam', 'Day2_US_FC2']
##
#Parameter to analyze
animal_folder = "FP117"
File_Name = "20250914_FP117_FC3_day2_Ch2"
ANYMAZEprotocol = para_FC2[0] #202303_Photometry+FC    202306_Photometry+FC2_reexam
DayForAnalysis = para_FC2[1]
pre_tag = 30
post_tag = 40
T_auc_stim_end = 20 # end time for calculating AUC after footshock (t=0)
T_poststim_start = '-' # specify time for AUC calculation, if not, put '-'
T_poststim_end = '-' # specify time for AUC calculation, if not, put '-'
bin_interval = 0.1 # XXXs interval per bin
plot_start = 0 # Type in the trace to plot if the stimulation points larger than 5 (Ext sessions)
plot_end = 5 # Type in the trace to plot if the stimulation points larger than 5 (Ext sessions)
channel = File_Name[-3:]
Fig_Title = "CeA astrocyte activity_FC3_OTRcKO" #OTR-astrocyte activity
output_dir = f"C:/Users/Kai-Yi.WANG/PycharmProjects/Photometry_analysis/Figure/{animal_folder}/{channel}/"

"""
Option for [DayForAnalysis]: DayX_XS_old; DayX_XS, DayX_XS_FC2
"""

##
# Load data for targeted channel and animal
data_folder = Path(f"C:/Users/Kai-Yi.WANG/PycharmProjects/Photometry_analysis/NPZ_file/{animal_folder}/")
Import_file = data_folder / File_Name
load_data = np.load(f'{Import_file}.npz', 'r', allow_pickle=True)

DFF, Target_CH, control_CH, Time = load_data['DFF'], load_data['Target_CH'], load_data['control_CH'], load_data['Time']

# extract the stimulation time for analysis
search_STIM = StimTimePoints(ANYMAZEprotocol)
STIM = search_STIM[DayForAnalysis]

##
# Extract traces at stimulus time
match_T_CH_DFF = dict(zip(Time, DFF))

all_stim_time = []
CH_stim_response = []

for US_time in STIM:
    eDff = {key: val for key, val in match_T_CH_DFF.items() if US_time-pre_tag <= key < US_time+post_tag}
    eDff_k = list(eDff.keys())
    eDff_v = list(eDff.values())
    all_stim_time.append(eDff_k)
    CH_stim_response.append(eDff_v)

##
# Normalized dF/F to a defined baseline period before the onset of CS
Baseline_ori_value = []
for US_time in STIM:
    eBL = {key: val for key, val in match_T_CH_DFF.items() if US_time-pre_tag <= key < US_time-20} # -20 is the onset of CS time
    eBL_v = list(eBL.values())
    Baseline_ori_value.append(eBL_v)

merge_ori_BL = np.concatenate(Baseline_ori_value)
BL_mean = np.average(merge_ori_BL)

# Subtract mean of defined baseline from each element in the nested list
nor_dFF_toBL = [[x - BL_mean for x in inner_list] for inner_list in CH_stim_response]

# Adjust each stimulus time point to 0
Nor_all_stim_time = []
for i in range(len(all_stim_time)):
    Nor_time = [realT - STIM[i] for realT in all_stim_time[i]]
    Nor_all_stim_time.append(Nor_time)

##
# find the max response and time to max after FC

Max_Resp = []
Max_Resp_time = []
for find_max in range(len(nor_dFF_toBL)):
    max_per_trace = max(nor_dFF_toBL[find_max])
    maxind_per_trace = nor_dFF_toBL[find_max].index(max_per_trace)
    max_time_per_trace = Nor_all_stim_time[find_max][maxind_per_trace]
    Max_Resp.append(max_per_trace)
    Max_Resp_time.append(max_time_per_trace)

avg_Max_Resp = np.mean(Max_Resp)
avg_Max_Resp_time = np.mean(Max_Resp_time)

Max_Resp.append(avg_Max_Resp) #with avg at the end
Max_Resp_time.append(avg_Max_Resp_time) #with avg at the end

# save into csv file
# add label in the first column
Max_Resp.insert(0, 'Peak response (dF/F)')
Max_Resp_time.insert(0, 'Peak latency (s)')

header = ('row name', 'trace_1', 'trace_2', 'trace_3', 'trace_4', 'trace_5', 'average')

with open(f"C:/Users/Kai-Yi.WANG/PycharmProjects/Photometry_analysis/export_csv/{File_Name}_peak_response.csv", 'w', newline ='') as file:
    f = csv.writer(file)
    f.writerow(header)
    f.writerow(Max_Resp)
    f.writerow(Max_Resp_time)

# print(Max_Resp)
# print(Max_Resp_time)
# print(avg_Max_Resp_time)
# print(avg_Max_Resp)

##
# Calculate z score (mean and SD from baseline)
bin_nor_dFF = []
bin_time_zS = []
N_Bin_forzS = (post_tag - (-pre_tag))/0.01
if len(Nor_all_stim_time) <= 5:
    for k in range(len(Nor_all_stim_time)):
        bin_means, bin_edges, binnumber = stats.binned_statistic(Nor_all_stim_time[k], nor_dFF_toBL[k],  statistic='mean', bins=N_Bin_forzS, range=(-pre_tag, post_tag))

        bin_nor_dFF.append(bin_means)
        binTime = np.delete(bin_edges, -1)
        bin_time_zS.append(binTime)

bin_nor_dFF = np.array(bin_nor_dFF)
bin_time_zS = np.array(bin_time_zS)

ext_BL_forzS = bin_nor_dFF[np.where((bin_time_zS >= -pre_tag) & (bin_time_zS < -20))]
mean_BL_forzS = np.mean(ext_BL_forzS)
std_BL_forzS = np.std(ext_BL_forzS)
zScore = [[(x - mean_BL_forzS) / std_BL_forzS for x in inner_list] for inner_list in bin_nor_dFF]

# For plotting AUC with zScore and x(time) value after binning (for zScore)
bin_time_plot_zSocre = np.delete(bin_edges, -1)
avg_bin_zScore = np.mean(zScore, axis=0)


##
# Calculate AUC
AUC_all = [[] for _ in range(4)]

for i in range(len(bin_time_zS)):
    AUC_T = bin_time_zS[i]
    AUC_S = bin_nor_dFF[i]

    # Define the time interval for AUC calculation
    AUC_time_FC = np.where((AUC_T >= 0) & (AUC_T <= T_auc_stim_end)) # **think about the end time
    total_auc_FC = np.trapz(AUC_S[AUC_time_FC], AUC_T[AUC_time_FC])
    AUC_FC = total_auc_FC - (mean_BL_forzS * len(AUC_time_FC))
    AUC_all[2].append(AUC_FC)

    AUC_time_tone = np.where((AUC_T >= -20) & (AUC_T < 0))
    total_auc_tone = np.trapz(AUC_S[AUC_time_tone], AUC_T[AUC_time_tone])
    AUC_tone = total_auc_tone - (mean_BL_forzS * len(AUC_time_tone))
    AUC_all[1].append(AUC_tone)

    AUC_time_BL = np.where((AUC_T >= -pre_tag) & (AUC_T < -20))
    total_auc_BL = np.trapz(AUC_S[AUC_time_BL], AUC_T[AUC_time_BL])
    AUC_BL = total_auc_BL - (mean_BL_forzS * len(AUC_time_BL))
    AUC_all[0].append(AUC_BL)


    if T_poststim_start == '-':
        # continue
        # plot for visualizing individual trace with AUC areas
        plt.plot(AUC_T, AUC_S)
        plt.fill_between(AUC_T[AUC_time_FC], AUC_S[AUC_time_FC], mean_BL_forzS, color='blue', alpha=0.3,
                         label='AUC')
        plt.fill_between(AUC_T[AUC_time_tone], AUC_S[AUC_time_tone], mean_BL_forzS, color='pink', alpha=0.3,
                         label='AUC')
        plt.fill_between(AUC_T[AUC_time_BL], AUC_S[AUC_time_BL], mean_BL_forzS, color='gray', alpha=0.3,
                         label='AUC')
        plt.xlabel('Time')
        plt.ylabel('Fluorescence')
        plt.title('Fiber Photometry Data')
        plt.legend()
        plt.savefig(output_dir + f'{File_Name}_Trace{i + 1}_AUC.png')
        plt.show()
    else:
        AUC_time_FC_after = np.where((AUC_T >= T_poststim_start) & (AUC_T <= T_poststim_end))
        total_auc_FC_after = np.trapz(AUC_S[AUC_time_FC_after], AUC_T[AUC_time_FC_after])
        AUC_FC_after = total_auc_FC_after - (mean_BL_forzS * len(AUC_time_FC_after))
        AUC_all[3].append(AUC_FC_after)

        # plot for visualizing individual trace with AUC areas (including poststim AUC)
        plt.plot(AUC_T, AUC_S)
        plt.fill_between(AUC_T[AUC_time_FC], AUC_S[AUC_time_FC], mean_BL_forzS, color='blue', alpha=0.3,
                         label='AUC')
        plt.fill_between(AUC_T[AUC_time_tone], AUC_S[AUC_time_tone], mean_BL_forzS, color='pink', alpha=0.3,
                         label='AUC')
        plt.fill_between(AUC_T[AUC_time_BL], AUC_S[AUC_time_BL], mean_BL_forzS, color='gray', alpha=0.3,
                         label='AUC')
        plt.fill_between(AUC_T[AUC_time_FC_after], AUC_S[AUC_time_FC_after], mean_BL_forzS, color='green', alpha=0.3,
                         label='AUC')
        plt.xlabel('Time')
        plt.ylabel('Fluorescence')
        plt.title('Fiber Photometry Data')
        plt.legend()
        plt.savefig(output_dir + f'{File_Name}_Trace{i+1}_AUC.png')
        # plt.show()


# print(AUC_all)
# print(np.mean(AUC_all[0]))
# print(np.mean(AUC_all[1]))
# print(np.mean(AUC_all[2]))
# print(np.mean(AUC_all[3]))

# save AUC with different time segments into csv file
AUC_avg_each = []
for avg_each_AUC in AUC_all:
    if len(avg_each_AUC) == 0:
        AUC_avg_each.append(None)
    else:
        avg_subAUC = np.mean(avg_each_AUC)
        AUC_avg_each.append(avg_subAUC)

time_period = [f'BL_-{pre_tag}_to_-20', f'tone_-20_to_0', f'FC_0_to_{T_auc_stim_end}', f'FC_{T_poststim_start}_to_{T_poststim_end}']
for i in range(len(time_period)):
    AUC_all[i].insert(0, time_period[i])
    AUC_all[i].append(AUC_avg_each[i])

header = ('period', 'trace_1', 'trace_2', 'trace_3', 'trace_4', 'trace_5', 'average')

with open(f"C:/Users/Kai-Yi.WANG/PycharmProjects/Photometry_analysis/AUC_export_csv/{File_Name}_AUC.csv", 'w', newline ='') as file:
    f = csv.writer(file)
    f.writerow(header)
    f.writerows(AUC_all)


##
# Binning data into the specific interval
Bin_CH = []
NumOfBin = (post_tag - (-pre_tag))/bin_interval
if len(Nor_all_stim_time) <= 5:
    for k in range(len(Nor_all_stim_time)):
        bin_means, bin_edges, binnumber = stats.binned_statistic(Nor_all_stim_time[k], nor_dFF_toBL[k],  statistic='mean', bins=NumOfBin, range=(-pre_tag, post_tag))
        # bin_means, bin_edges, binnumber = stats.binned_statistic(Nor_all_stim_time[k], zScore_wholetrace[k],
        #                                                          statistic='mean', bins=NumOfBin,
        #                                                          range=(-pre_tag, post_tag))
        Bin_CH.append(bin_means)
else:
    for k in range(plot_start, plot_end):
        bin_means, bin_edges, binnumber = stats.binned_statistic(Nor_all_stim_time[k], nor_dFF_toBL[k],  statistic='mean', bins=NumOfBin, range=(-pre_tag, post_tag))
        # bin_means, bin_edges, binnumber = stats.binned_statistic(Nor_all_stim_time[k], zScore_wholetrace[k],
        #                                                          statistic='mean', bins=NumOfBin,
        #                                                          range=(-pre_tag, post_tag))
        Bin_CH.append(bin_means)

# Calculate average and sem of signals
avg_Bin_CH = np.mean(Bin_CH, axis=0)
sem_Bin_CH = stats.sem(Bin_CH, axis=0)

# For plotting the x(time) value after binning
bin_time_plot = np.delete(bin_edges, -1)

# save the time and avg extracted response into npz
np.savez(f"C:/Users/Kai-Yi.WANG/PycharmProjects/Photometry_analysis/export_csv/{File_Name}_avgTrace", Bin_T = bin_time_plot, Bin_R = avg_Bin_CH)

# write the time and avg extracted response into csv
with open(f"C:/Users/Kai-Yi.WANG/PycharmProjects/Photometry_analysis/export_csv/{File_Name}_avgTrace.csv", 'w', newline ='') as file:
    f = csv.writer(file)
    for a in range(len(avg_Bin_CH)):
        f.writerow([bin_time_plot[a], avg_Bin_CH[a]])



##
# Visualize the fitting result
plt.figure(1)
plt.plot(Time, Target_CH, label='Raw data')
plt.plot(Time, control_CH, 'r-', label='Fitting result')
plt.xlabel("Time (s)")
plt.ylabel("Voltage")
plt.title("Fitting Result")
plt.savefig(output_dir + f'{File_Name}_Fitting_Result.png')

# plot continuous time vs dF/F responses
plt.figure(2)
for j in range(len(STIM)):
    plt.axvline(STIM[j], linestyle='--', color='0.6')
    plt.axvline(STIM[j]-20, linestyle='--', color='r', alpha=0.3)
# plt.axvline(960, linestyle='--', color=(0,0.8,0.2))
plt.plot(Time, DFF)
plt.xlabel("Time (s)")
plt.ylabel("ΔF/F%")
plt.title(Fig_Title)
plt.savefig(output_dir + f'{File_Name}_markSTIM.png')

# # Plot of the extracted responses (dF/F) vs recording time
# plt.figure(3)
# for i in range(len(all_stim_time)):
#     # if i == 4:
#         plt.plot(all_stim_time[i], CH_stim_response[i], label = 'CeA-astrocyte')
#         # plt.plot(all_US_time[i], C3_US_response[i], label = 'PVN-OT')
#         plt.axvline(STIM[i], linestyle='--', color='0.6')
#         plt.axvline(STIM[i]-20, linestyle='--', color=(0.8,0,0,0.2))
#         plt.xlabel("Time (s)")
#         plt.ylabel("ΔF/F%")
#         plt.title(Fig_Title)
# plt.savefig(output_dir + f'{File_Name}_recT_dFF.png')
#
# # Plot of the extracted responses (dF/F) vs normalized time
# # plt.figure(4)
# # for k in range(len(Nor_all_stim_time)):
# #     plt.plot(Nor_all_stim_time[k], CH_stim_response[k], label='CeA-astrocyte')
# #     # plt.plot(all_US_time[i], C3_US_response[i], label = 'PVN-OT')
# # plt.axvline(0, linestyle='--', color='0.6')
# # plt.axvline(-20, linestyle='--', color=(0.8,0,0,0.2))
# # plt.xlabel("Time (s)")
# # plt.ylabel("ΔF/F%")
# # plt.title(Fig_Title)
# # # plt.savefig(output_dir + f'{File_Name}_norT_dFF.png')
#
# plt.figure(5)
# for j in range(len(Bin_CH)):
#     # if j == 4:
#         plt.plot(bin_time_plot, Bin_CH[j], label = j, alpha=0.6)
#     #     plt.plot(all_US_time[i], C3_US_response[i], label = 'PVN-OT')
# plt.axvline(0, linestyle='--', color='0.6')
# plt.axvline(-20, linestyle='--', color=(0.8,0,0,0.2))
# plt.xlabel("Time (s)")
# plt.ylabel("ΔF/F%")
# plt.title(Fig_Title)
# plt.legend()
# plt.savefig(output_dir + f'{File_Name}_norT_ind.png')
#
# plt.figure(6)
# for j in range(len(zScore)):
#     # if i == 0:
#         plt.plot(bin_time_zS[j], zScore[j], label = j, alpha=0.5)
#         # plt.plot(all_US_time[i], C3_US_response[i], label = 'PVN-OT')
# plt.axvline(0, linestyle='--', color='0.6')
# plt.axvline(-20, linestyle='--', color=(0.8,0,0,0.2))
# plt.xlabel("Time (s)")
# plt.ylabel("zScore")
# plt.title(Fig_Title)
# plt.legend()
# plt.savefig(output_dir + f'{File_Name}_norT_zScore.png')

##
plt.figure(7)
plt.plot(bin_time_plot, avg_Bin_CH, label='CeA-astrocyte')
plt.fill_between(bin_time_plot, avg_Bin_CH + sem_Bin_CH, avg_Bin_CH - sem_Bin_CH, alpha=0.4)
# plt.plot(all_US_time[i], C3_US_response[i], label = 'PVN-OT')
plt.axvline(0, linestyle='--', color='0.6')
plt.axvline(-20, linestyle='--', color=(0.8,0,0,0.2))
# plt.ylim(-3, 15) #!!!!
# plt.xlim(-10, 20) #!!!!
plt.axvspan(-20, 0,  color='pink', alpha=0.3)
plt.axvspan(0, 1,  color='yellow', alpha=0.8)
plt.xlabel("Time (s)")
plt.ylabel("ΔF/F%")
plt.title(Fig_Title)
plt.savefig(output_dir + f'{File_Name}_avg.png')
# plt.show()

#
# bin avg plot AUC
plt.figure(8)
plt.plot(bin_time_plot_zSocre, avg_bin_zScore)
plt.fill_between(bin_time_plot_zSocre[AUC_time_FC], avg_bin_zScore[AUC_time_FC], mean_BL_forzS, color='blue', alpha=0.3,
                 label='Baseline_AUC')
plt.fill_between(bin_time_plot_zSocre[AUC_time_tone], avg_bin_zScore[AUC_time_tone], mean_BL_forzS, color='pink', alpha=0.3,
                 label='Tone_AUC')
plt.fill_between(bin_time_plot_zSocre[AUC_time_BL], avg_bin_zScore[AUC_time_BL], mean_BL_forzS, color='gray', alpha=0.3,
                 label='FC_AUC')
# plt.fill_between(bin_time_plot_zSocre[AUC_time_FC_after], avg_bin_zScore[AUC_time_FC_after], mean_BL_forzS, color='green', alpha=0.3,
#                  label='afterFC_AUC')
plt.axvline(0, linestyle='--', color='0.6')
plt.axvline(-20, linestyle='--', color=(0.8,0,0,0.2))
# plt.ylim(-4, 2) #!!!!
plt.xlabel("Time (s)")
plt.ylabel("zScore")
plt.title(Fig_Title)
plt.savefig(output_dir + f'{File_Name}_avg_AUC.png')
plt.show()

##

