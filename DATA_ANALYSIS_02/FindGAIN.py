# FindGAIN.py

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
    Files    = []

    min_peak = 0
    max_peak = 200
    nbins    = 200
    binw     = (max_peak - min_peak)/nbins
    peaksHy  = np.zeros((len(HV), nbins))
    peaksHx  = np.arange(min_peak, max_peak, binw)

    titleHx = "Peak (mV)"
    titleHy = "Normalized Counts"



    with open(file, "r") as infilelist:
        for f in infilelist:
            f = f.rstrip('\n')
            print("Reading file ", f)
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


    # Plot Hist
    plt.figure(figsize=(10, 6))
    plt.semilogy(peaksHx, peaksHy[0], color='black',   linestyle='-')
    plt.semilogy(peaksHx, peaksHy[1], color='#964B00', linestyle='-')
    plt.semilogy(peaksHx, peaksHy[2], color='red',     linestyle='-')
    plt.semilogy(peaksHx, peaksHy[3], color='#FFD700', linestyle='-')
    plt.semilogy(peaksHx, peaksHy[4], color='green',   linestyle='-')
    plt.semilogy(peaksHx, peaksHy[5], color='cyan',    linestyle='-')
    plt.semilogy(peaksHx, peaksHy[6], color='blue',    linestyle='-')
    plt.xlabel(titleHx)
    plt.ylabel(titleHy)
    plt.grid(True)

    # Print on File
    print("\n\nFiles Analized:\n")
    for i in range(len(Files)):
        print(Files[i])

    plt.show(block=False)
    input("Press Enter to continue... ")
    plt.close()


@jit(nopython=True)
def FillHistPeaks(value, vector, min_value, max_value, binw):
    for i in range(len(vector)):
        if(value < min_value + i*binw):
            vector[i] = vector[i] + 1
            break

if __name__ == '__main__':
    args = PARSER.parse_args()
    main(**args.__dict__)
