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

# Other files in the same project:
import tekwfm
import PlotSettings


__description__ = 'ReadFileFindPeaks'

formatter = argparse.ArgumentDefaultsHelpFormatter
PARSER = argparse.ArgumentParser(description=__description__, formatter_class=formatter)
PARSER.add_argument('-f', '--input_file', type=str, required=False, default='file.wfm', help='File Name to be analyzed (.wfm)')


def main(**kwargs):

    file = kwargs['input_file']
    dleddt = 2*12



    # Read File
    volts, tstart, tscale, tfrac, tdatefrac, tdate = tekwfm.read_wfm(file)
    # Create time vector
    samples, frames = volts.shape
    tstop = samples * tscale + tstart
    t = np.linspace(tstart, tstop, num=samples, endpoint=False)
    times = np.zeros(volts.shape)
    for frame, subsample in enumerate(tfrac):
        toff = subsample * tscale
        times[:,frame] = t + toff
    # Conversions
    times *= 1e9  # times in ns
    volts *= 1e3  # volts in mV
    # DLED
    times_dled, volts_dled = DLED(times, volts, dleddt)
    # Print Introduction
    print(file)
    print("Peaks found with find_peaks(volts_dled[:,i], height=0)")
    print("NO trace smoothing")
    print("TimeWindow = [" + str(times[0,0]) + " , " + str(times[-1,0]) + "] ns,  Dt = " + str(times[-1,0]-times[0,0]) + " ns" )
    print("DLED related: dleddt = " + str(dleddt))
    print("END_INTRODUCTION")
    input("Press Enter to continue... ")

    # Find Peaks
    for i in range(frame + 1):
    # for i in range(3):
        indices, _ = find_peaks(volts_dled[:,i], height=0)
        for j in range(len(indices)):
            print(str(times_dled[indices[j],i]) + "\t" + str(volts_dled[indices[j],i]))
        print("N")
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
