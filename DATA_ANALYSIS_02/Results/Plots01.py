#!/usr/bin/env python3

'''
--------------------------------------------------------------------------------
|   Plots.py                                                                   |
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

    color = ['#FFD700', 'red', '#FF00FF']

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

    ###################################
    ########   START VALUES   #########
    ###################################
    DCR_Area[0] = [294.55239308829647, 401.10780094547437, 481.5594748880893, 549.9570029635022, 614.8710647683184, 675.2750223590923, 736.7957314908331] # kHz
    errDCR_Area[0] = [3.0386181141030875, 1.78966237339165, 3.1525066032232303, 4.799358061076759, 8.803392168289065, 10.27524087267318, 8.725943178186299] # kHz
    CT[0] = [0.2148099778252095, 0.23307444372324945, 0.2776512795526091, 0.3147572595004262, 0.3457156244956133, 0.3885942546231689, 0.4256236428823916]
    errCT[0] = [0.01336197231057995, 0.012500654629939656, 0.004584626037726358, 0.004243450880588029, 0.00720313669766709, 0.005928814192721454, 0.00760855024417878]
    DCR_Area_Del[0] = [293.9010277781799, 409.0305565036102, 497.9039056028509, 570.5985435567897, 640.553507256302, 702.3456667860879, 766.2647799571932] # kHz
    errDCR_Area_Del[0] = [3.2654540272793042, 1.6575884897989113, 3.0528927577615264, 7.221939717736518, 16.103445331254193, 19.09845410264745, 16.443833933366932] # kHz
    CT_Del[0] = [0.24029019676947436, 0.25055339352501815, 0.28563870711662065, 0.3100126528602453, 0.3244324074330121, 0.35800378945745315, 0.3857695160986809]
    errCT_Del[0] = [0.008470631175901994, 0.012963429199379928, 0.006327938325368099, 0.004952056957554407, 0.010142416561450651, 0.007724756694236279, 0.010001186921550997]

    DCR_Area[1] = [296.847346199881, 410.3793891501269, 511.5150706459789, 588.1319627964064, 658.3614664238403, 724.9408376814074, 787.4768893407528] # kHz
    errDCR_Area[1] = [5.827023412838685, 3.331583230179433, 10.307706088634177, 14.669925831273758, 19.67312417043422, 19.870341858409688, 12.706634076664841] # kHz
    CT[1] = [0.19544637223324274, 0.21687582393124347, 0.2640399128444377, 0.30081257044216403, 0.3297298127387233, 0.3772110863550284, 0.4164885625996709]
    errCT[1] = [0.01493337512423859, 0.016207807327326573, 0.0058997847487703425, 0.006264431061789644, 0.012122407433773963, 0.008625180870549642, 0.010616284033044099]
    DCR_Area_Del[1] = [299.48309069103595, 412.97481959117647, 518.8118383704758, 601.2958378800636, 675.0196988230734, 743.1854271735966, 809.9829646254701] # kHz
    errDCR_Area_Del[1] = [6.437385222592184, 3.5864450065351434, 10.69461070533356, 18.248642955495598, 29.72675914910269, 32.256706847008445, 23.22202791326481] # kHz
    CT_Del[1] = [0.21841137185090037, 0.2323226305404412, 0.2682587271920454, 0.2874103940747575, 0.3008703524470603, 0.33909462301652626, 0.369369676588554]
    errCT_Del[1] = [0.014686122719732131, 0.020929838570120124, 0.010495854557787787, 0.009110002648859583, 0.01887475652747328, 0.012745070173772965, 0.01542339576621704]

    DCR_Area[2] = [359.39516497231796, 478.429222571692, 575.8075759771831, 659.0522837055855, 738.3088370081614, 813.3034045282803, 886.8962105508891] # kHz
    errDCR_Area[2] = [3.2858639099752054, 1.5244063552482316, 3.1571332137751824, 6.1651820899828635, 12.589856067517985, 16.0603767746353, 14.004122229816517] # kHz
    CT[2] = [0.22607862338017926, 0.24725468704458595, 0.2929115607038434, 0.3299667217425706, 0.36189307895050604, 0.40501236954279574, 0.4432044612308677]
    errCT[2] = [0.013763407856059728, 0.012545310492329514, 0.005326304358050149, 0.004653347640773742, 0.008293790805483536, 0.006870571792700342, 0.009075423979364006]
    DCR_Area_Del[2] = [353.6160663376103, 480.1130440684083, 586.7492370182435, 676.4479521438074, 756.8861052670111, 834.2502386030953, 908.3138153793002] # kHz
    errDCR_Area_Del[2] = [3.8554800465604444, 1.4890912435445784, 3.6008060604840466, 10.862716082820043, 24.97204843786767, 32.1317446357408, 28.360546512396013] # kHz
    CT_Del[2] = [0.23886716562925067, 0.25320747969843993, 0.2843428246541058, 0.30735502775605583, 0.3236660230877801, 0.3573989363934111, 0.38571488533575576]
    errCT_Del[2] = [0.009806586179418736, 0.011612483837787901, 0.006114319757976994, 0.005988202324395953, 0.011787863820013122, 0.009440293587447024, 0.012159606241237852]

    GAIN[0]  = np.array([9.532071122414779, 13.139574248904314, 16.314732832489092, 19.242722621276968, 21.996580651998517, 24.86324611982192, 27.576760517468827])
    GAIN[1]  = np.array([7.1056072764475875, 11.882550134545795, 14.851595435912348, 17.63779762162774, 20.265761421116025, 22.934266147655894, 25.514369755498596])
    GAIN[2]  = np.array([9.718380716038025, 13.420000142430673, 16.533648866707814, 19.572260588104687, 22.434069524835564, 25.273139971643666, 27.983371569666005])
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
        plt.errorbar(HV, DCR_Area_Del[i], xerr=errHV, yerr=errDCR_Area_Del[i], color=color[i], fmt='s', markersize=4, label=label)
    plt.xlabel(titleHV)
    plt.ylabel(titleDCR)
    plt.legend()

    # Plot CT CNT and Del
    plt.figure(figsize=(10, 6))
    for i in range(nFilesTot):
        label = "FBK NUV HD3-2 (" + str(i+1) + ")" + " CNT"
        plt.errorbar(HV, CT[i], xerr=errHV, yerr=errCT[i], color=color[i], fmt='o', markersize=4, label=label)
        label = "FBK NUV HD3-2 (" + str(i+1) + ")" + " Del"
        plt.errorbar(HV, CT_Del[i], xerr=errHV, yerr=errCT_Del[i], color=color[i], fmt='s', markersize=4, label=label)
    plt.xlabel(titleHV)
    plt.ylabel(titleCT)
    plt.legend()



    plt.show(block=False)
    input("Press Enter to continue... ")
    plt.close()

if __name__ == '__main__':
    args = PARSER.parse_args()
    main(**args.__dict__)
