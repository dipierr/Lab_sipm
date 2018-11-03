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
    GAIN = np.zeros((nFilesTot, len(HV)))
    errGAIN = np.zeros((nFilesTot, len(HV)))

    ###################################
    ########   START VALUES   #########
    ###################################
    DCR_Area[0] = [294.55239308829647, 401.10780094547437, 481.5594748880893, 549.9570029635022, 614.8710647683184, 675.2750223590923, 736.7957314908331] # kHz
    errDCR_Area[0] = [3.0386181141030875, 1.78966237339165, 3.1525066032232303, 4.799358061076759, 8.803392168289065, 10.27524087267318, 8.725943178186299] # kHz
    CT[0] = [0.2148099778252095, 0.23307444372324945, 0.2776512795526091, 0.3147572595004262, 0.3457156244956133, 0.3885942546231689, 0.4256236428823916]
    errCT[0] = [0.01336197231057995, 0.012500654629939656, 0.004584626037726358, 0.004243450880588029, 0.00720313669766709, 0.005928814192721454, 0.00760855024417878]

    DCR_Area[1] = [296.847346199881, 410.3793891501269, 511.5150706459789, 588.1319627964064, 658.3614664238403, 724.9408376814074, 787.4768893407528] # kHz
    errDCR_Area[1] = [5.827023412838685, 3.331583230179433, 10.307706088634177, 14.669925831273758, 19.67312417043422, 19.870341858409688, 12.706634076664841] # kHz
    CT[1] = [0.19544637223324274, 0.21687582393124347, 0.2640399128444377, 0.30081257044216403, 0.3297298127387233, 0.3772110863550284, 0.4164885625996709]
    errCT[1] = [0.01493337512423859, 0.016207807327326573, 0.0058997847487703425, 0.006264431061789644, 0.012122407433773963, 0.008625180870549642, 0.010616284033044099]

    DCR_Area[2] = [359.39516497231796, 478.429222571692, 575.8075759771831, 659.0522837055855, 738.3088370081614, 813.3034045282803, 886.8962105508891] # kHz
    errDCR_Area[2] = [3.2858639099752054, 1.5244063552482316, 3.1571332137751824, 6.1651820899828635, 12.589856067517985, 16.0603767746353, 14.004122229816517] # kHz
    CT[2] = [0.22607862338017926, 0.24725468704458595, 0.2929115607038434, 0.3299667217425706, 0.36189307895050604, 0.40501236954279574, 0.4432044612308677]
    errCT[2] = [0.013763407856059728, 0.012545310492329514, 0.005326304358050149, 0.004653347640773742, 0.008293790805483536, 0.006870571792700342, 0.009075423979364006]

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
        plt.errorbar(HV, DCR_Area[i], xerr=errHV, yerr=errDCR_Area[i], color=color[i], fmt='o', markersize=4)
    plt.xlabel(titleHV)
    plt.ylabel(titleDCR)

    # Plot CT
    plt.figure(figsize=(10, 6))
    for i in range(nFilesTot):
        plt.errorbar(HV, CT[i], xerr=errHV, yerr=errCT[i], color=color[i], fmt='o', markersize=4)
    plt.xlabel(titleHV)
    plt.ylabel(titleCT)

    # Plot GAIN
    plt.figure(figsize=(10, 6))
    for i in range(nFilesTot):
        plt.errorbar(HV, GAIN[i], xerr=errHV, yerr=errGAIN[i], color=color[i], fmt='o', markersize=4)
    plt.xlabel(titleHV)
    plt.ylabel(titleGAIN)


    plt.show(block=False)
    input("Press Enter to continue... ")
    plt.close()

if __name__ == '__main__':
    args = PARSER.parse_args()
    main(**args.__dict__)
