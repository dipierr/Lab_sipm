# FindDCR.py

import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
from scipy.optimize import curve_fit
from scipy import asarray as ar,exp
from scipy.odr import *
from pylab import *
import argparse

# Other files in the same project:
import PlotSettings

__description__ = 'Plot X Y'

formatter = argparse.ArgumentDefaultsHelpFormatter
PARSER = argparse.ArgumentParser(description=__description__, formatter_class=formatter)
PARSER.add_argument('-f', '--input_file', type=str, required=False, default='file.txt', help='File to be analyzed.')


def main(**kwargs):

    Area = 36 # mm^2
    TimeWindow = 1024 # ns
    dleddt = 0 # ns

    # Import Plot Settings from PlotSettings.py:
    PlotSettings.PlotSettings()

    titleX = "$X_{title}$"
    titleY = "$Y_{title}$"

    # Open File for XY plot
    file = open(kwargs['input_file'], "r")
    lines = file.read().split("\n")

    print("\nReading file:")
    print(kwargs['input_file'] + '\n')

    time = np.array([])
    peaks = np.array([])

    nTraks = 0

    introduction_bool = True

    for l in lines:
        if l == "END_INTRODUCTION":
            introduction_bool = False
        if introduction_bool:
            if "dleddt" in l:
                dleddt = double(re.findall(r"[-+]?\d*\.\d+|\d+", l)[0])
            print(l)
        else:
            if l == "N":
                nTraks = nTraks + 1
            else:
                temp = re.findall(r"[-+]?\d*\.\d+|\d+", l)
                if len(temp) > 0:
                    time = np.append(time,float(temp[0]))
                    peaks = np.append(peaks,float(temp[1]))

    DCR = len(peaks)/(TimeWindow*nTraks-len(peaks)*2*dleddt)*1e3
    DCR_Area = DCR/Area*1e3
    print("DCR      = " + str(DCR) + "\t MHz")
    print("DCR/Area = " + str(DCR_Area) + "\t kHz")

    # Plot
    plt.figure(figsize=(10, 6))
    plt.plot(time, peaks, color='blue', marker='o', linestyle='None', markersize=4)
    plt.xlabel(titleX)
    plt.ylabel(titleY)
    plt.grid(True)

    plt.show(block=False)
    input("Press Enter to continue... ")
    plt.close()

if __name__ == '__main__':
    args = PARSER.parse_args()
    main(**args.__dict__)
