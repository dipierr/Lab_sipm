#!/usr/bin/env python3

'''
--------------------------------------------------------------------------------
|   FindGAIN.py                                                                |
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

__description__ = 'Find GAIN'

formatter = argparse.ArgumentDefaultsHelpFormatter
PARSER = argparse.ArgumentParser(description=__description__, formatter_class=formatter)
PARSER.add_argument('-flist', '--input_filelist', type=str, required=False, default='file.txt', help='File to be analyzed.')

def main(**kwargs):
    file = kwargs['input_filelist']
    ResultsFileName = str(file) + "_GAIN.txt"
    Area = 36 # mm^2
    TimeWindow = 1024 # ns
    dleddt = 0 # ns
    nFile = 0

    HV       = np.array([31, 32, 33, 34, 35, 36, 37])
    errHV = np.full(len(HV), 0.01)
    Files    = []

    min_peak = 0
    max_peak = 200
    nbins    = 200
    binw     = (max_peak - min_peak)/nbins
    peaksHy  = np.zeros((len(HV), nbins))
    peaksHx  = np.arange(min_peak, max_peak, binw)

    titleHx = "Peaks (mV)"
    titleHy = "Normalized Counts"



    with open(file, "r") as infilelist:
        for f in infilelist:
            f = f.rstrip('\n')
            print("\nReading file:\n" + FindNameAndPath(f))
            Files.append(f)
            with open(f, "r") as file:
                lines = file.read().split("\n")

                nTraks = 0
                introduction_bool = True

                for l in lines:
                    if l == "END_INTRODUCTION":
                        introduction_bool = False
                        continue
                    if introduction_bool:
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
                                FillHistPeaks(float(temp[1]), peaksHy[nFile], min_peak, max_peak, binw)

                nFile = nFile + 1


    for i in range(len(HV)):
        peaksHy[i] = peaksHy[i]/sum(peaksHy[i])

    # Import Plot Settings from PlotSettings.py:
    PlotSettings.PlotSettings()

    color = ['black', '#964B00', 'red', '#FFD700', 'green', 'cyan', 'blue']

    # Plot Hists
    for i in range(len(HV)):
        plt.figure(figsize=(10, 6))
        plt.step(peaksHx, peaksHy[i], color=color[i],   linestyle='-')
        plt.xlabel(titleHx)
        plt.ylabel(titleHy)

    # Plot Hist at different HVs
    plt.figure(figsize=(10, 6))
    for i in range(len(HV)):
        plt.step(peaksHx, peaksHy[i], color=color[i],   linestyle='-')
    plt.yscale("log")
    plt.xlabel(titleHx)
    plt.ylabel(titleHy)

    # Print on File
    ResultsFile = open(ResultsFileName, "w")
    ResultsFile.write("Files Analized:\n")
    for i in range(len(Files)):
        ResultsFile.write(FindNameAndPath(Files[i]) + "\n")
    ResultsFile.write("\n\nParameters:\n")
    ResultsFile.write("HV = " + str(list(HV)) + " # V \n")
    ResultsFile.write("errHV = " + str(list(errHV)) + " # V \n")
    ResultsFile.write("\n\nResults:\n")
    ResultsFile.close()

    plt.show(block=False)
    input("Press Enter to continue... ")
    plt.close()


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
