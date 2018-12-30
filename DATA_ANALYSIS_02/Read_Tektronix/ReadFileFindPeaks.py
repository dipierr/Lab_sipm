'''
--------------------------------------------------------------------------------
|   ReadFileFindPeaks.py                                                               |
|                                                                              |
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
from scipy.signal import find_peaks
import argparse
from numba import jit
import datetime
import time
import ntpath

# Other files in the same project:
import tekwfm
import PlotSettings


__description__ = 'ReadFileFindPeaks'

formatter = argparse.ArgumentDefaultsHelpFormatter
PARSER = argparse.ArgumentParser(description=__description__, formatter_class=formatter)
PARSER.add_argument('-f', '--input_file', type=str, required=False, default='file.wfm', help='File Name to be analyzed (.wfm)')


def main(**kwargs):

    File = kwargs['input_file']
    FilePath, FileName = ntpath.split(File)
    FilePeaksName = FilePath + "/Peaks/" + FileName + "_Peaks.txt"
    dleddt = 2*12

    # Read File
    volts, tstart, tscale, tfrac, tdatefrac, tdate = tekwfm.read_wfm(File)
    # Create time vector
    samples, frames = volts.shape
    tstop = samples * tscale + tstart
    t = np.linspace(tstart, tstop, num=samples, endpoint=False)
    times = np.zeros(volts.shape)
    for frame, subsample in enumerate(tfrac):
        toff = subsample * tscale
        times[:,frame] = t + toff
    # Conversions
    k_times = 1e9
    k_volts = 1e3
    times *= k_times  # times in ns
    volts *= k_volts  # volts in mV
    # DLED
    times_dled, volts_dled = DLED(times, volts, dleddt)
    # Open File Peaks
    FilePeaks = open(FilePeaksName, "w")
    # Print Introduction
    FilePeaks.write("Peaks from file:\n")
    FilePeaks.write(FileName + "\n")
    FilePeaks.write("Sampling = " + str(1/(tscale*k_times)) + " GHz     (Time Resolution = " + str(tscale*k_times*1e3) + " ps)\n")
    FilePeaks.write("Peaks found with find_peaks(volts_dled[:,i], height=0)\n")
    FilePeaks.write("NO trace smoothing\n")
    FilePeaks.write("TimeWindow = [" + str(tstart*k_times) + " , " + str(tstop*k_times) + "] ns,  Dt = " + str((tstop-tstart)*k_times) + " ns\n" )
    FilePeaks.write("DLED related: dleddt = " + str(dleddt) + " indices\n")
    FilePeaks.write("END_INTRODUCTION\n")

    # Find Peaks
    for i in range(frame + 1):
    # for i in range(3):
        indices, _ = find_peaks(volts_dled[:,i], height=0)
        for j in range(len(indices)):
            FilePeaks.write(str(times_dled[indices[j],i]) + "\t" + str(volts_dled[indices[j],i]) + "\n")
        FilePeaks.write("N\n")
        # Show the Plots
        # ShowPlot(times[:,i], volts[:,i], times_dled[:,i], volts_dled[:,i], indices)



def DLED(times, volts, dleddt):
    times_dled = times[dleddt:,:]
    volts_dled = volts[dleddt:,:] - volts[:-dleddt,:]
    return times_dled, volts_dled



def ShowPlot(x1, y1, x2, y2, indices):
    plt.figure()
    plt.plot(x1, y1, "b")
    plt.plot(x2, y2, "r")
    plt.plot(x2[indices], y2[indices], "x")
    plt.xlabel("Time (ns)")
    plt.ylabel("Amplitude (mV)")
    plt.show(block=False)
    input("Press Enter to continue... ")
    plt.close()





if __name__ == '__main__':
    args = PARSER.parse_args()
    main(**args.__dict__)
