# from collections import namedtuple

# Time tag for stimulation (for ANYMAZE protocol 202301_Photometry+FC)
Day1_CS_old = [84, 146, 238, 359, 458]
Day2_CS_old = [40, 132, 269, 351, 462]
Day2_US_old = [60, 152, 289, 371, 482]
Day3_CS_old = [103, 210, 303, 372, 497, 567, 641, 716, 825, 895, 964, 1098, 1190, 1253, 1361, 1469, 1552, 1659,
               1781, 1922]
Day4_CS_old = [87, 176, 293, 395, 491, 607, 714, 783, 893, 966, 1087, 1199, 1322, 1406, 1500, 1615, 1688, 1827,
               1912, 1995]
Day5_CS_old = [96, 200, 291]

# Time tag for stimulation (for ANYMAZE protocol 202303_Photometry+FC)
TimeStamp = []
for New_time in [Day1_CS_old, Day2_CS_old, Day2_US_old, Day3_CS_old, Day4_CS_old, Day5_CS_old]:
    NewTime = [t + 120 for t in New_time]
    TimeStamp.append(NewTime)

Day1_CS, Day2_CS, Day2_US, Day3_CS, Day4_CS, Day5_CS = [TimeStamp[i] for i in range(len(TimeStamp))]

# Time tag for stimulation "FC2" (for ANYMAZE protocol 202306_Photometry+FC2_reexam / 202307_Photometry+RemoteRet [Day1 only])
Day1_CS_FC2 = [196, 275, 374, 493, 554]
Day2_CS_FC2 = [184, 321, 382, 511, 611]
Day2_US_FC2 = [i + 20 for i in Day2_CS_FC2]
Day3_CS_FC2 = [199, 326, 405, 516, 613]

tail_lifting_FP38 = [120, 240, 360, 485, 600, 722, 840, 960, 1080, 1200]
tail_lifting_FP39 = [120, 240, 363, 480, 600, 720, 862, 960, 1082, 1200]

Rat7_FSI2_LCeA_0102 = [141, 152, 211, 227, 236, 241, 289, 301, 321, 336, 368, 379, 393, 409, 499, 525, 548, 568, 591]
Rat8_FSI2_LCeA_0102 = [135, 152, 156, 184, 202, 205, 225, 276, 340, 357, 394, 495, 545]
Rat9_FSI2_LCeA_0102 = [146, 172, 178, 194, 197, 216, 229, 251, 276, 316, 330, 403, 465, 489]
Rat7_FSI3_LCeA_CNO_0104 = [134, 150, 158, 164, 171, 181, 193, 199, 215, 223, 233, 301, 321, 344, 394, 401, 420, 430, 444, 489, 516, 552, 589]
Rat8_FSI3_LCeA_CNO_0104 = [147, 154, 163, 176, 181, 240, 258, 270, 278, 289, 352, 366, 383, 453, 458, 483, 574, 583]
Rat8_FSI1_RaIC_1231 = [135, 142, 160, 170, 187, 223, 228, 239, 268, 279, 346, 360, 372, 427, 474]
Rat9_FSI1_RaIC_1231 = [137, 146, 152, 159, 168, 250, 348, 353, 381, 414, 427, 437, 441, 477, 517, 534, 541]

Rat8_FSI2_LCeA_0102_touch = [141, 146, 225, 497]
Rat8_FSI2_LCeA_0102_facetouch = [184, 214, 266, 276, 341, 356]

Rat7_pinch_LCeA = [83, 135, 180, 233, 263, 290, 339, 407, 566, 618]
Rat8_pinch_LCeA = [85, 114, 154, 231, 268, 288, 348, 513, 587, 612]
Rat4_pinch_LaIC = [76, 125, 207, 269, 359]

def StimTimePoints(ANYMAZEprotocol):
    if ANYMAZEprotocol == '202301_Photometry+FC':
        Stim_dict = {'Day1_CS_old': Day1_CS_old,
                     'Day2_CS_old': Day2_CS_old,
                     'Day2_US_old': Day2_US_old,
                     'Day3_CS_old': Day3_CS_old,
                     'Day4_CS_old': Day4_CS_old,
                     'Day5_CS_old': Day5_CS_old}
        return Stim_dict

    elif ANYMAZEprotocol == '202303_Photometry+FC':
        Stim_dict = {'Day1_CS': Day1_CS,
                     'Day2_CS': Day2_CS,
                     'Day2_US': Day2_US,
                     'Day3_CS': Day3_CS,
                     'Day4_CS': Day4_CS,
                     'Day5_CS': Day5_CS}
        return Stim_dict

    elif ANYMAZEprotocol == '202306_Photometry+FC2_reexam' or ANYMAZEprotocol == '202307_Photometry+RemoteRet':
        Stim_dict = {'Day1_CS_FC2': Day1_CS_FC2,
                     'Day2_CS_FC2': Day2_CS_FC2,
                     'Day2_US_FC2': Day2_US_FC2,
                     'Day3_CS_FC2': Day3_CS_FC2}
        return Stim_dict

    elif ANYMAZEprotocol == '202306_forGRABOTR':
        Stim_dict = {'Day1_CS': Day1_CS,
                     'Day2_CS': Day2_CS,
                     'Day2_US': Day2_US,
                     'Day3_CS': Day3_CS[0:5]}
        return Stim_dict

    elif ANYMAZEprotocol == 'tail_lifting':
        Stim_dict = {'tail_lifting_FP38': tail_lifting_FP38,
                     'tail_lifting_FP39': tail_lifting_FP39}
        return Stim_dict

    elif ANYMAZEprotocol == 'FSI':
        Stim_dict = {'Rat7_FSI2_LCeA_0102': Rat7_FSI2_LCeA_0102,
                     'Rat8_FSI2_LCeA_0102': Rat8_FSI2_LCeA_0102,
                     'Rat9_FSI2_LCeA_0102': Rat9_FSI2_LCeA_0102,
                     'Rat7_FSI3_LCeA_CNO_0104': Rat7_FSI3_LCeA_CNO_0104,
                     'Rat8_FSI3_LCeA_CNO_0104': Rat8_FSI3_LCeA_CNO_0104,
                     'Rat8_FSI1_RaIC_1231': Rat8_FSI1_RaIC_1231,
                     'Rat9_FSI1_RaIC_1231': Rat9_FSI1_RaIC_1231,
                     'Rat8_FSI2_LCeA_0102_touch': Rat8_FSI2_LCeA_0102_touch,
                     'Rat8_FSI2_LCeA_0102_facetouch': Rat8_FSI2_LCeA_0102_facetouch}
        return Stim_dict

    elif ANYMAZEprotocol == 'Pinch':
        Stim_dict = {'Rat7_pinch_LCeA': Rat7_pinch_LCeA,
                     'Rat8_pinch_LCeA': Rat8_pinch_LCeA,
                     'Rat4_pinch_LaIC': Rat4_pinch_LaIC}
        return Stim_dict

##
# # from collections import namedtuple
#
# # Time tag for stimulation (for remote fear memory)
# RFM_Hab = [72, 173, 268, 380, 504]
# RFM_FC_CS = [111, 195, 306, 386, 477]
# RFM_FC_US = [i + 19 for i in Day2_CS_FC2]
# RFM_Ret = [85, 171, 282, 412, 516]
#
#
# tail_lifting_FP38 = [120, 240, 360, 485, 600, 722, 840, 960, 1080, 1200]
# tail_lifting_FP39 = [120, 240, 363, 480, 600, 720, 862, 960, 1082, 1200]
#
#
# def StimTimePoints(ANYMAZEprotocol):
#     if ANYMAZEprotocol == 'RFM_Photometry':
#         Stim_dict = {'RFM_Hab': RFM_Hab,
#                      'RFM_FC_CS': RFM_FC_CS,
#                      'RFM_FC_US': RFM_FC_US,
#                      'RFM_Ret': RFM_Ret}
#         return Stim_dict


##
# def StimTimePoints(ANYMAZEprotocol):
#     if ANYMAZEprotocol == '202301_Photometry+FC':
#         DaySelected = namedtuple('DaySelected',
#                                  ['Day1_CS_old', 'Day2_CS_old', 'Day2_US_old', 'Day3_CS_old', 'Day4_CS_old', 'Day5_CS_old'])
#         return DaySelected(Day1_CS_old, Day2_CS_old, Day2_US_old, Day3_CS_old, Day4_CS_old, Day5_CS_old)
#     elif ANYMAZEprotocol == '202303_Photometry+FC':
#         DaySelected = namedtuple('DaySelected',
#                                  ['Day1_CS', 'Day2_CS', 'Day2_US', 'Day3_CS', 'Day4_CS', 'Day5_CS'])
#         return DaySelected(Day1_CS, Day2_CS, Day2_US, Day3_CS, Day4_CS, Day5_CS)
#     elif ANYMAZEprotocol == '202306_Photometry+FC2_reexam':
#         DaySelected = namedtuple('DaySelected',
#                                  ['Day1_CS_FC2', 'Day2_CS_FC2', 'Day2_US_FC2', 'Day3_CS_FC2'])
#         return DaySelected(Day1_CS_FC2, Day2_CS_FC2, Day2_US_FC2, Day3_CS_FC2)

##
# Day2_CS_FC2 = []
# a = 120
# for j in range(5):
#     b = [79, 107, 59, 91, 77]
#     if j == 0:
#         a += b[j]
#     if j > 0:
#         a += b[j] + 20
#     Day2_CS_FC2.append(a)