# FindDCRCT.py

import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
from scipy.optimize import curve_fit
from scipy import asarray as ar,exp
from scipy.odr import *
from pylab import *
import argparse

from numba import jit

# Other files in the same project:
import PlotSettings

__description__ = 'Find DCR and CT'

formatter = argparse.ArgumentDefaultsHelpFormatter
PARSER = argparse.ArgumentParser(description=__description__, formatter_class=formatter)
PARSER.add_argument('-flist', '--input_filelist', type=str, required=False, default='file.txt', help='File to be analyzed.')

def main(**kwargs):
    file = kwargs['input_filelist']
    Area = 36 # mm^2
    TimeWindow = 1024 # ns
    dleddt = 0 # ns
    nFile = 0

    HV       = np.array([31, 32, 33, 34, 35, 36, 37])
    thr_0_5  = np.array([8,  8,  8,  9,  9,  9,  10])
    thr_1_5  = np.array([19, 23, 26, 30, 35, 37, 40])
    DCR_Area = np.array([])
    CT       = np.array([])
    thrs     = np.arange(8,50,1)
    DCR_thr  = np.zeros((len(HV), len(thrs)))
    peaks    = np.array([])


    titleHV = "HV (V)"
    titleDCR = "$\\frac{\si{DCR}}{\si{Area}} \, (\\frac{\si{\kilo\hertz}}{\si{mm^2}})$"
    titleCT = "Cross Talk"
    titleThrs = "Thresholds (mV)"
    titleDCR_multi = "$\\frac{\si{DCR}}{\si{Area}} \, (\\frac{\si{\hertz}}{\si{mm^2}})$"

    with open(file, "r") as infilelist:
        for f in infilelist:
            f = f.rstrip('\n')
            print("Reading file ", f)
            with open(f, "r") as file:
                lines = file.read().split("\n")

                nTraks = 0
                nEvents = 0
                cnt_0_5_pe = 0
                cnt_1_5_pe = 0
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
                                cnt_0_5_pe, cnt_1_5_pe = FindDCRthrs(peak, thr_0_5, thr_1_5, cnt_0_5_pe, cnt_1_5_pe, thrs, len(thrs), DCR_thr, nFile)

                TimeWindowCorrected = TimeWindow*nTraks - nEvents*2*dleddt
                DCR = cnt_0_5_pe/(TimeWindowCorrected)*1e3
                DCR_Area = np.append(DCR_Area, float(DCR/Area*1e3))
                CT = np.append(CT, float(cnt_1_5_pe/cnt_0_5_pe))
                DCR_thr[nFile] =  DCR_thr[nFile]/(TimeWindowCorrected*1e-9*Area)
                nFile = nFile + 1



    # Import Plot Settings from PlotSettings.py:
    PlotSettings.PlotSettings()

    # Plot DCR
    plt.figure(figsize=(10, 6))
    plt.plot(HV, DCR_Area, color='blue', marker='o', linestyle='None', markersize=4)
    plt.xlabel(titleHV)
    plt.ylabel(titleDCR)
    plt.grid(True)

    # Plot CT
    plt.figure(figsize=(10, 6))
    plt.plot(HV, CT, color='orange', marker='o', linestyle='None', markersize=4)
    plt.xlabel(titleHV)
    plt.ylabel(titleCT)
    plt.grid(True)

    # Plot DCR vs thrs
    plt.figure(figsize=(10, 6))
    plt.semilogy(thrs, DCR_thr[0], color='black',  linestyle='-')
    plt.semilogy(thrs, DCR_thr[1], color='#964B00',   linestyle='-')
    plt.semilogy(thrs, DCR_thr[2], color='red',    linestyle='-')
    plt.semilogy(thrs, DCR_thr[3], color='yellow',  linestyle='-')
    plt.semilogy(thrs, DCR_thr[4], color='green',   linestyle='-')
    plt.semilogy(thrs, DCR_thr[5], color='cyan',    linestyle='-')
    plt.semilogy(thrs, DCR_thr[6], color='blue',    linestyle='-')
    plt.xlabel(titleThrs)
    plt.ylabel(titleDCR_multi)
    plt.grid(True)

    plt.show(block=False)
    input("Press Enter to continue... ")
    plt.close()

@jit
def FindDCRthrs(peak, thr_0_5, thr_1_5, cnt_0_5_pe, cnt_1_5_pe, thrs, len, DCR_thr, nFile):
    if(peak > thr_0_5[nFile]):
        cnt_0_5_pe = cnt_0_5_pe + 1
        if(peak > thr_1_5[nFile]):
            cnt_1_5_pe = cnt_1_5_pe + 1
    for i in range(len):
        if (peak > thrs[i]):
            DCR_thr[nFile][i] = DCR_thr[nFile][i] + 1
        else:
            break
    return cnt_0_5_pe, cnt_1_5_pe

if __name__ == '__main__':
    args = PARSER.parse_args()
    main(**args.__dict__)
