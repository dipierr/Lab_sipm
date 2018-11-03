#!/usr/bin/env python3

'''
--------------------------------------------------------------------------------
|   FindDCRCT.py                                                               |
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

__description__ = 'Find DCR and CT'

formatter = argparse.ArgumentDefaultsHelpFormatter
PARSER = argparse.ArgumentParser(description=__description__, formatter_class=formatter)
PARSER.add_argument('-flist', '--input_filelist', type=str, required=False, default='file.txt', help='File to be analyzed.')

def main(**kwargs):
    file = kwargs['input_filelist']
    ResultsFileName = str(file) + "_DCRCT.txt"
    Area = 36 # mm^2
    TimeWindow = 1024 # ns
    dleddt = 0 # ns
    nFile = 0

    HV    = np.array([31, 32, 33, 34, 35, 36, 37])
    errHV = np.full(len(HV), 0.01)

    thr_0_5  = np.zeros((3, len(HV)))
    thr_1_5  = np.zeros((3, len(HV)))
    DCR_Area = np.zeros((3, len(HV)))
    CT       = np.zeros((3, len(HV)))

    thr_0_5[0]  = [ 8,  8,  8,  8,  8,  8,  8]
    thr_0_5[1]  = [ 8,  8,  8,  9,  9,  9, 10]
    thr_0_5[2]  = [ 9,  9, 10, 12, 14, 15, 15]

    thr_1_5[0]  = [18.5, 22, 26, 29, 33, 35, 37]
    thr_1_5[1]  = [19,   23, 26, 30, 35, 37, 40]
    thr_1_5[2]  = [19.5, 24, 27, 31, 36, 39, 42]

    thrs     = np.arange(8,70,1)
    DCR_thr  = np.zeros((len(HV), len(thrs)))
    Files    = []

    titleHV = "HV (V)"
    titleDCR = "$\\frac{\si{DCR}}{\si{Area}} \, (\\frac{\si{\kilo\hertz}}{\si{mm^2}})$"
    titleCT = "Cross Talk"
    titleThrs = "Thresholds (mV)"
    titleDCR_multi = "$\\frac{\si{DCR}}{\si{Area}} \, (\\frac{\si{\hertz}}{\si{mm^2}})$"


    with open(file, "r") as infilelist:
        for f in infilelist:
            f = f.rstrip('\n')
            print("\nReading file ", FindNameAndPath(f))
            Files.append(f)
            with open(f, "r") as file:
                lines = file.read().split("\n")

                nTraks  = 0
                nEvents = 0
                cnt_0_5_pe = np.zeros(3)
                cnt_1_5_pe = np.zeros(3)
                introduction_bool = True

                for l in lines:
                    if l == "END_INTRODUCTION":
                        introduction_bool = False
                        continue
                    if introduction_bool:
                        if "dleddt" in l:
                            dleddt = double(re.findall(r"[-+]?\d*\.\d+|\d+", l)[0])
                        print(l)
                    else:
                        if l == "N":
                            if (nTraks%10000 == 0):
                                print("Read trace " + str(nTraks))
                            nTraks = nTraks + 1
                        else:
                            # temp = re.findall(r"[-+]?\d*\.\d+|\d+", l)
                            temp = l.split("\t")
                            if len(temp) > 1:
                                nEvents = nEvents + 1
                                peak = float(temp[1])
                                FindDCRthrs(peak, thr_0_5, thr_1_5, cnt_0_5_pe, cnt_1_5_pe, thrs, DCR_thr, nFile)


                TimeWindowCorrected = TimeWindow*nTraks - nEvents*2*dleddt
                DCR = cnt_0_5_pe/(TimeWindowCorrected)*1e3
                for i in range(3):
                    DCR_Area[i][nFile] = float(DCR[i]/Area*1e3)
                    CT[i][nFile] = float(cnt_1_5_pe[i]/cnt_0_5_pe[i])
                DCR_thr[nFile] =  DCR_thr[nFile]/(TimeWindowCorrected*1e-9*Area)
                nFile = nFile + 1


    errDCR_Area = np.maximum(np.abs(DCR_Area[1] - DCR_Area[0]), np.abs(DCR_Area[1] - DCR_Area[2]))
    errCT = np.maximum(np.abs(CT[1] - CT[0]), np.abs(CT[1] - CT[2]))

    # Import Plot Settings from PlotSettings.py:
    PlotSettings.PlotSettings()

    # Plot DCR
    plt.figure(figsize=(10, 6))
    plt.errorbar(HV, DCR_Area[1], xerr=errHV, yerr=errDCR_Area, color='blue', fmt='o', markersize=4)
    plt.xlabel(titleHV)
    plt.ylabel(titleDCR)

    # Plot CT
    plt.figure(figsize=(10, 6))
    plt.errorbar(HV, CT[1], xerr=errHV, yerr=errCT, color='orange', fmt='o', markersize=4)
    plt.xlabel(titleHV)
    plt.ylabel(titleCT)

    # Plot DCR vs thrs
    color = ['black', '#964B00', 'red', '#FFD700', 'green', 'cyan', 'blue']
    plt.figure(figsize=(10, 6))
    for i in range(len(HV)):
        label = str(HV[i]) + " V"
        plt.semilogy(thrs, DCR_thr[i], color=color[i], linestyle='-', label=label)
    plt.xlabel(titleThrs)
    plt.ylabel(titleDCR_multi)
    plt.legend()

    # Print on File
    ResultsFile = open(ResultsFileName, "w")
    ResultsFile.write("Files Analized:\n")
    for i in range(len(Files)):
        ResultsFile.write(FindNameAndPath(Files[i]) + "\n")
    ResultsFile.write("\n\nParameters:\n")
    ResultsFile.write("HV = " + str(list(HV)) + " # V \n")
    ResultsFile.write("errHV = " + str(list(errHV)) + " # V \n")
    for i in range(3):
        ResultsFile.write("thr_0_5[" + str(i) + "] = " + str(list(thr_0_5[i])) + "\n")
    for i in range(3):
        ResultsFile.write("thr_1_5[" + str(i) + "] = " + str(list(thr_1_5[i])) + "\n")
    ResultsFile.write("\n\nResults:\n")
    ResultsFile.write("DCR_Area = " + str(list(DCR_Area[1])) + " # kHz \n")
    ResultsFile.write("errDCR_Area = " + str(list(errDCR_Area)) + " # kHz \n")
    ResultsFile.write("CT = " + str(list(CT[1])) + "\n")
    ResultsFile.write("errCT = " + str(list(errCT)) + "\n")
    ResultsFile.close()

    plt.show(block=False)
    input("Press Enter to continue... ")
    plt.close()

@jit(nopython=True)
def FindDCRthrs(peak, thr_0_5, thr_1_5, cnt_0_5_pe, cnt_1_5_pe, thrs, DCR_thr, nFile):
    # DCR at 0.5 pe:
    for i in range(3):
        if (peak > thr_0_5[i][nFile]):
            cnt_0_5_pe[i] = cnt_0_5_pe[i] + 1
        else:
            break
    # DCR at 1.5 pe:
    for i in range(3):
        if (peak > thr_1_5[i][nFile]):
            cnt_1_5_pe[i] = cnt_1_5_pe[i] + 1
        else:
            break
    # DCR at different thresholds:
    for i in range(len(thrs)):
        if (peak > thrs[i]):
            DCR_thr[nFile][i] = DCR_thr[nFile][i] + 1
        else:
            break


@jit(nopython=True)
def FillHistPeaks(value, vector, min_value, max_value, binw):
    for i in range(len(vector)):
        if(value < min_value + i*binw):
            vector[i] = vector[i] + 1
            break

def FindNameAndPath(path):
    head, tail = ntpath.split(path)
    return tail or ntpath.basename(head)

if __name__ == '__main__':
    args = PARSER.parse_args()
    main(**args.__dict__)
