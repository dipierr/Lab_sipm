'''
--------------------------------------------------------------------------------
|   Plots_20180725.py                                                          |
|                                                                              |
|   Remember to give permissions to the file:                                  |
|       $ chmod +x *.py                                                        |
|                                                                              |
|   Davide Depaoli 2018                                                        |
--------------------------------------------------------------------------------
'''

'''
From files:
    20180725_HD3-2_01_DARK_AgilentE3641A_31.00_AS_2_100000ev_01.dat_Peaks_8.00.txt
    20180725_HD3-2_01_DARK_AgilentE3641A_32.00_AS_2_100000ev_01.dat_Peaks_8.00.txt
    20180725_HD3-2_01_DARK_AgilentE3641A_33.00_AS_2_100000ev_01.dat_Peaks_8.00.txt
    20180725_HD3-2_01_DARK_AgilentE3641A_34.00_AS_2_100000ev_01.dat_Peaks_8.00.txt
    20180725_HD3-2_01_DARK_AgilentE3641A_35.00_AS_2_100000ev_01.dat_Peaks_8.00.txt
    20180725_HD3-2_01_DARK_AgilentE3641A_36.00_AS_2_100000ev_01.dat_Peaks_8.00.txt
    20180725_HD3-2_01_DARK_AgilentE3641A_37.00_AS_2_100000ev_01.dat_Peaks_8.00.txt

    20180725_HD3-2_02_DARK_AgilentE3641A_31.00_AS_2_100000ev_01.dat_Peaks_8.00.txt
    20180725_HD3-2_02_DARK_AgilentE3641A_32.00_AS_2_100000ev_01.dat_Peaks_8.00.txt
    20180725_HD3-2_02_DARK_AgilentE3641A_33.00_AS_2_100000ev_01.dat_Peaks_8.00.txt
    20180725_HD3-2_02_DARK_AgilentE3641A_34.00_AS_2_100000ev_01.dat_Peaks_8.00.txt
    20180725_HD3-2_02_DARK_AgilentE3641A_35.00_AS_2_100000ev_01.dat_Peaks_8.00.txt
    20180725_HD3-2_02_DARK_AgilentE3641A_36.00_AS_2_100000ev_01.dat_Peaks_8.00.txt
    20180725_HD3-2_02_DARK_AgilentE3641A_37.00_AS_2_100000ev_01.dat_Peaks_8.00.txt

    20180725_HD3-2_03_DARK_AgilentE3641A_31.00_AS_2_100000ev_01.dat_Peaks_8.00.txt
    20180725_HD3-2_03_DARK_AgilentE3641A_32.00_AS_2_100000ev_01.dat_Peaks_8.00.txt
    20180725_HD3-2_03_DARK_AgilentE3641A_33.00_AS_2_100000ev_01.dat_Peaks_8.00.txt
    20180725_HD3-2_03_DARK_AgilentE3641A_34.00_AS_2_100000ev_01.dat_Peaks_8.00.txt
    20180725_HD3-2_03_DARK_AgilentE3641A_35.00_AS_2_100000ev_01.dat_Peaks_8.00.txt
    20180725_HD3-2_03_DARK_AgilentE3641A_36.00_AS_2_100000ev_01.dat_Peaks_8.00.txt
    20180725_HD3-2_03_DARK_AgilentE3641A_37.00_AS_2_100000ev_01.dat_Peaks_8.00.txt
'''


import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
from scipy.optimize import curve_fit
from scipy import asarray as ar,exp
from scipy.odr import *
from pylab import *
import argparse
import ntpath

from numba import jit

# Other files in the same project:
import PlotSettings

__description__ = 'Some Plots'

formatter = argparse.ArgumentDefaultsHelpFormatter
PARSER = argparse.ArgumentParser(description=__description__, formatter_class=formatter)

def main(**kwargs):

    nFilesTot = 3

    # Import Plot Settings from PlotSettings.py:
    PlotSettings.PlotSettings()

    color  = ['#FFD700', 'red', '#FF00FF']
    colorD = ['#FF8C00', '#8B0000', '#8A2BE2']

    HV = [31, 32, 33, 34, 35, 36, 37] # V
    errHV = [0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01] # V

    DCR_Area = np.zeros((nFilesTot, len(HV)))
    errDCR_Area = np.zeros((nFilesTot, len(HV)))
    CT = np.zeros((nFilesTot, len(HV)))
    errCT = np.zeros((nFilesTot, len(HV)))
    DCR_Area_Del = np.zeros((nFilesTot, len(HV)))
    errDCR_Area_Del = np.zeros((nFilesTot, len(HV)))
    CT_Del = np.zeros((nFilesTot, len(HV)))
    errCT_Del = np.zeros((nFilesTot, len(HV)))
    GAIN = np.zeros((nFilesTot, len(HV)))
    errGAIN = np.zeros((nFilesTot, len(HV)))

    Vbd = np.zeros(nFilesTot)

    ###################################
    ########   START VALUES   #########
    ###################################
    DCR_Area[0] = [294.55239308829647, 401.10780094547437, 481.5594748880893, 549.9570029635022, 614.8710647683184, 675.2750223590923, 736.7957314908331] # kHz
    errDCR_Area[0] = [3.0386181141030875, 1.78966237339165, 3.1525066032232303, 4.799358061076759, 8.803392168289065, 10.27524087267318, 8.725943178186299] # kHz
    CT[0] = [0.2148099778252095, 0.23307444372324945, 0.2776512795526091, 0.3147572595004262, 0.3457156244956133, 0.3885942546231689, 0.4256236428823916]
    errCT[0] = [0.01336197231057995, 0.012500654629939656, 0.004584626037726358, 0.004243450880588029, 0.00720313669766709, 0.005928814192721454, 0.00760855024417878]
    DCR_Area_Del[0] = [272.74344113795297, 380.42720849184417, 464.98887694524484, 532.499853999949, 594.6706238660611, 651.9969510665262, 708.1157090229007] # kHz
    errDCR_Area_Del[0] = [2.966118393524823, 1.4227868722359176, 2.7187353104793033, 6.786212874018361, 14.801194062299032, 17.444852818896834, 14.59177247571381] # kHz
    CT_Del[0] = [0.2394530749504907, 0.2548734164832318, 0.29557237944209974, 0.3250611924543762, 0.34074501014823616, 0.3776364339716225, 0.4076174777541066]
    errCT_Del[0] = [0.009779508110357521, 0.01568543594482169, 0.00720633927118175, 0.005665735851918929, 0.011731655476765412, 0.008747482092274605, 0.01139913637176182]


    ###############
    DCR_Area[1] = [296.847346199881, 410.3793891501269, 511.5150706459789, 588.1319627964064, 658.3614664238403, 724.9408376814074, 787.4768893407528] # kHz
    errDCR_Area[1] = [5.827023412838685, 3.331583230179433, 10.307706088634177, 14.669925831273758, 19.67312417043422, 19.870341858409688, 12.706634076664841] # kHz
    CT[1] = [0.19544637223324274, 0.21687582393124347, 0.2640399128444377, 0.30081257044216403, 0.3297298127387233, 0.3772110863550284, 0.4164885625996709]
    errCT[1] = [0.01493337512423859, 0.016207807327326573, 0.0058997847487703425, 0.006264431061789644, 0.012122407433773963, 0.008625180870549642, 0.010616284033044099]
    DCR_Area_Del[1] = [271.89514373346725, 380.7215354566125, 475.26944458879683, 551.761475214314, 618.982823785793, 679.2082437319115, 738.4057332572467] # kHz
    errDCR_Area_Del[1] = [5.743688454710195, 3.02451513443998, 8.410867738606385, 15.47139835036296, 26.051619657540982, 28.026725056602118, 20.087153944710508] # kHz
    CT_Del[1] = [0.2104209148658025, 0.23358999904512381, 0.2786444423452531, 0.3015271098574084, 0.3147558681764255, 0.35969474861028217, 0.39239692034287454]
    errCT_Del[1] = [0.01543589905375578, 0.023426770789291923, 0.013644200632375147, 0.011135787396807484, 0.022407434201921095, 0.015153625802228954, 0.01796038354991869]

    ###############
    DCR_Area[2] = [359.39516497231796, 478.429222571692, 575.8075759771831, 659.0522837055855, 738.3088370081614, 813.3034045282803, 886.8962105508891] # kHz
    errDCR_Area[2] = [3.2858639099752054, 1.5244063552482316, 3.1571332137751824, 6.1651820899828635, 12.589856067517985, 16.0603767746353, 14.004122229816517] # kHz
    CT[2] = [0.22607862338017926, 0.24725468704458595, 0.2929115607038434, 0.3299667217425706, 0.36189307895050604, 0.40501236954279574, 0.4432044612308677]
    errCT[2] = [0.013763407856059728, 0.012545310492329514, 0.005326304358050149, 0.004653347640773742, 0.008293790805483536, 0.006870571792700342, 0.009075423979364006]
    DCR_Area_Del[2] = [325.9481181066461, 443.14670604434417, 540.8731063431654, 621.0070898550982, 691.9363825667543, 758.2650004791124, 819.4434763623906] # kHz
    errDCR_Area_Del[2] = [3.4768138692421076, 1.3514164602023584, 3.2163775937924584, 9.799861727669281, 21.91140650345642, 27.76688217664207, 23.708674856866764] # kHz
    CT_Del[2] = [0.23957521625430883, 0.26220029508694576, 0.29946130529605874, 0.3261904776594619, 0.34484267895180104, 0.3839084101720149, 0.4159578949914569]
    errCT_Del[2] = [0.010104300956073814, 0.013881900333676256, 0.007282994241000484, 0.006671619225388115, 0.01329724513182734, 0.01044506456362615, 0.013505550901086238]


    ###############
    GAIN[0] = [9.532068027089933, 13.139574896383952, 16.314733573439213, 19.24271583423075, 21.996579499060108, 24.86324611982192, 27.576764941287024]
    Vbd[0] = 27.62668949381624

    GAIN[1] = [7.105609603031427, 11.882549769831478, 14.851595299041708, 17.637796668693014, 20.265761352780117, 22.934266147655894, 25.514368833877757]
    Vbd[1] = 28.18968736995698

    GAIN[2] = [9.718378350080377, 13.420000150043908, 16.533649811619675, 19.572265182161193, 22.434070408823594, 25.273139971643666, 27.98337237545008]
    Vbd[2] = 27.605109633725824


    ###################################
    #########   END VALUES   ##########
    ###################################


    titleHV = "HV (V)"
    titleDCR = "$\\frac{\si{DCR}}{\si{Area}} \, (\\frac{\si{\kilo\hertz}}{\si{mm^2}})$"
    titleCT = "Cross Talk"
    titleGAIN = "GAIN (mV)"


    # Plot DCR
    plt.figure(figsize=(10, 6))
    for i in range(nFilesTot):
        label = "FBK NUV HD3-2 (" + str(i+1) + ")"
        plt.errorbar(HV, DCR_Area[i], xerr=errHV, yerr=errDCR_Area[i], color=color[i], fmt='o', markersize=4, label=label)
    plt.xlabel(titleHV)
    plt.ylabel(titleDCR)
    plt.legend()

    # Plot CT
    plt.figure(figsize=(10, 6))
    for i in range(nFilesTot):
        label = "FBK NUV HD3-2 (" + str(i+1) + ")"
        plt.errorbar(HV, CT[i], xerr=errHV, yerr=errCT[i], color=color[i], fmt='o', markersize=4, label=label)
    plt.xlabel(titleHV)
    plt.ylabel(titleCT)
    plt.legend()

    # Plot GAIN
    plt.figure(figsize=(10, 6))
    for i in range(nFilesTot):
        label = "FBK NUV HD3-2 (" + str(i+1) + ")"
        plt.errorbar(HV, GAIN[i], xerr=errHV, yerr=errGAIN[i], color=color[i], fmt='o', markersize=4, label=label)
    plt.xlabel(titleHV)
    plt.ylabel(titleGAIN)
    plt.legend()

    ####################################

    # Plot DCR CNT and Del
    plt.figure(figsize=(10, 6))
    for i in range(nFilesTot):
        label = "FBK NUV HD3-2 (" + str(i+1) + ")" + " CNT"
        plt.errorbar(HV, DCR_Area[i], xerr=errHV, yerr=errDCR_Area[i], color=color[i], fmt='o', markersize=4, label=label)
        label = "FBK NUV HD3-2 (" + str(i+1) + ")" + " Del"
        plt.errorbar(HV, DCR_Area_Del[i], xerr=errHV, yerr=errDCR_Area_Del[i], color=colorD[i], fmt='s', markersize=4, label=label)
    plt.xlabel(titleHV)
    plt.ylabel(titleDCR)
    plt.legend()

    # Plot CT CNT and Del
    plt.figure(figsize=(10, 6))
    for i in range(nFilesTot):
        label = "FBK NUV HD3-2 (" + str(i+1) + ")" + " CNT"
        plt.errorbar(HV, CT[i], xerr=errHV, yerr=errCT[i], color=color[i], fmt='o', markersize=4, label=label)
        label = "FBK NUV HD3-2 (" + str(i+1) + ")" + " Del"
        plt.errorbar(HV, CT_Del[i], xerr=errHV, yerr=errCT_Del[i], color=colorD[i], fmt='s', markersize=4, label=label)
    plt.xlabel(titleHV)
    plt.ylabel(titleCT)
    plt.legend()


    # Vbd
    plt.figure(figsize=(10, 6))
    plt.hist(np.array(Vbd), 10)
    plt.xlabel("$V_{bd}$")


    plt.show(block=False)
    input("Press Enter to continue... ")
    plt.close()

if __name__ == '__main__':
    args = PARSER.parse_args()
    main(**args.__dict__)
