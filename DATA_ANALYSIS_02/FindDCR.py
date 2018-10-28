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
PARSER.add_argument('-flist', '--input_filelist', type=str, required=False, default='file.txt', help='File to be analyzed.')


def main(**kwargs):

    Area = 36 # mm^2
    TimeWindow = 1024 # ns
    dleddt = 0 # ns

    HV = np.array([36,37])
    DCR_Area = np.array([])

    titleX = "HV (V)"
    titleY = "DCR $\\frac{\si{\kilo\hertz}}{\si{mm^2}}$"

    with open(kwargs['input_filelist'], "r") as infilelist:
        for f in infilelist:
            f = f.rstrip('\n')
            print("Reading file ", f)
            with open(f, "r") as file:
                lines = file.read().split("\n")

                time = np.array([])
                peaks = np.array([])
                nTraks = 0
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
                            if (nTraks%1000 == 0):
                                print("Read trace " + str(nTraks))
                            nTraks = nTraks + 1
                        else:
                            # temp = re.findall(r"[-+]?\d*\.\d+|\d+", l)
                            temp = l.split("\t")
                            if len(temp) > 1:
                                time = np.append(time,float(temp[0]))
                                peaks = np.append(peaks,float(temp[1]))

                DCR = len(peaks)/(TimeWindow*nTraks-len(peaks)*2*dleddt)*1e3
                DCR_Area = np.append(DCR_Area, float(DCR/Area*1e3))
                print("DCR      = " + str(DCR) + "\t MHz")


    # Import Plot Settings from PlotSettings.py:
    PlotSettings.PlotSettings()

    # Plot
    plt.figure(figsize=(10, 6))
    plt.plot(HV, DCR_Area, color='blue', marker='o', linestyle='None', markersize=4)
    plt.xlabel(titleX)
    plt.ylabel(titleY)
    plt.grid(True)

    plt.show(block=False)
    input("Press Enter to continue... ")
    plt.close()

if __name__ == '__main__':
    args = PARSER.parse_args()
    main(**args.__dict__)
