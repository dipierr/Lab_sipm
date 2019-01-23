'''
--------------------------------------------------------------------------------
|   FindGAIN.py                                                                |
|                                                                              |
|   Remember to give permissions to the file:                                  |
|       $ chmod +x *.py                                                        |
|                                                                              |
|   Davide Depaoli 2018 -2019                                                  |
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
from numba import jit
import ntpath
import datetime

# Other files in the same project:
import PlotSettings

__description__ = 'Find GAIN'

formatter = argparse.ArgumentDefaultsHelpFormatter
PARSER = argparse.ArgumentParser(description=__description__, formatter_class=formatter)
PARSER.add_argument('-flist', '--input_filelist', type=str, required=False, default='file.txt', help='File to be analyzed')
PARSER.add_argument('-NoFit', '--no_fit_gain', action='store_true', required=False, default=False, help='Disables histogram fitting')


def main(**kwargs):
    file = kwargs['input_filelist']
    ResultsFileName = str(file) + "_GAIN.txt"
    Area = 36 # mm^2
    TimeWindow = 1024 # ns
    dleddt = 0 # ns
    nFile = 0

    fit_bool = True

    HV    = np.array([30, 31, 32, 33, 34, 35, 36, 37])
    errHV = np.full(len(HV), 0.01)
    Files = []

    min_peak = 0
    max_peak = 200
    nbins    = 400
    binw     = (max_peak - min_peak)/nbins
    peaksHy  = np.zeros((len(HV), nbins))
    peaksHx  = np.arange(min_peak, max_peak, binw)

    titleHx   = "Peaks (mV)"
    titleHy   = "Normalized Counts"
    titleHV   = "HV (V)"
    titleGAIN = "GAIN (mV)"


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

    # Normalize Hists
    for i in range(len(HV)):
        peaksHy[i] = peaksHy[i]/sum(peaksHy[i])


    # Import Plot Settings from PlotSettings.py:
    PlotSettings.PlotSettings()

    color = ['black', '#964B00', 'red', '#FFD700', 'green', 'cyan', 'blue', '#5000FF']

    gaus1low   = [14,   14,  13,  13,  13,  13,  14, 40]
    gaus1high  = [24,  24,  25,  27,  32,  33,  35, 40]
    gaus2low   = [30.0,  30,  27,  30,  34,  40,  42, 40]
    gaus2high  = [40.0,  40,  42,  47,  52,  60,  65, 40]
    a1Init     = np.full(len(HV), 1e-3)
    mean1Init  = (np.array(gaus1high) + np.array(gaus1low)) / 2
    sigma1Init = (np.array(gaus1high) - np.array(gaus1low)) / 2
    a2Init     = np.full(len(HV), 1e-3)
    mean2Init  = (np.array(gaus2high) + np.array(gaus2low)) / 2
    sigma2Init = (np.array(gaus2high) - np.array(gaus2low)) / 2

    GAIN       = np.zeros(len(HV))
    errGAIN    = np.zeros(len(HV))

    # Plot Hists
    for i in range(len(HV)):
        plt.figure(figsize=(10, 6))
        plt.plot(peaksHx, peaksHy[i], color=color[i],   linestyle='-')
        plt.xlabel(titleHx)
        plt.ylabel(titleHy)
        # Fit 1 Init
        yfit1 = peaksHy[i][peaksHx>gaus1low[i]]
        xfit1 = peaksHx[peaksHx>gaus1low[i]]
        yfit1 = yfit1[xfit1<gaus1high[i]]
        xfit1 = xfit1[xfit1<gaus1high[i]]
        # Fit 2 Init
        yfit2 = peaksHy[i][peaksHx>gaus2low[i]]
        xfit2 = peaksHx[peaksHx>gaus2low[i]]
        yfit2 = yfit2[xfit2<gaus2high[i]]
        xfit2 = xfit2[xfit2<gaus2high[i]]
        if(kwargs["no_fit_gain"] == False):
            popt1, pcov1 = curve_fit(gaus, xfit1, yfit1, p0=[a1Init[i], mean1Init[i], sigma1Init[i]])
            plt.plot(peaksHx, gaus(peaksHx, *popt1), color='red', linestyle='--')
            popt2, pcov2 = curve_fit(gaus, xfit2, yfit2, p0=[a2Init[i], mean2Init[i], sigma2Init[i]])
            plt.plot(peaksHx, gaus(peaksHx, *popt2), color='red', linestyle='--')
            GAIN[i] = popt2[1] - popt1[1]
            errGAIN[i] = pcov1[1][1] + pcov2[1][1]


    # Find Vbd
    poptVbd, pcovVbd = curve_fit(line, HV, GAIN)
    Vbd = - poptVbd[1] / poptVbd[0]
    print("Vbd = " + str(Vbd))

    # Plot GAIN
    plt.figure(figsize=(10, 6))
    plt.errorbar(HV, GAIN, xerr=errHV, yerr=errGAIN, color='blue', fmt='o', markersize=4)
    plt.xlabel(titleHV)
    plt.ylabel(titleGAIN)

    # Plot Hist at different HVs
    plt.figure(figsize=(10, 6))
    for i in range(len(HV)):
        label = str(HV[i]) + " V"
        plt.step(peaksHx, peaksHy[i], color=color[i], linestyle='-', label=label)
    plt.yscale("log")
    plt.xlabel(titleHx)
    plt.ylabel(titleHy)
    plt.legend()

    # Print on File
    ResultsFile = open(ResultsFileName, "w")
    ResultsFile.write(str(datetime.datetime.today().strftime('%Y-%m-%d')) + "    " + str(datetime.datetime.now().time()) + "\n")
    ResultsFile.write("\n\nFiles Analized:\n")
    for i in range(len(Files)):
        ResultsFile.write(FindNameAndPath(Files[i]) + "\n")
    ResultsFile.write("\n\nParameters:\n")
    ResultsFile.write("HV    = " + str(list(HV)) + " # V \n")
    ResultsFile.write("errHV = " + str(list(errHV)) + " # V \n")
    ResultsFile.write("gaus1low   = " + str(list(gaus1low)) + " \n")
    ResultsFile.write("gaus1high  = " + str(list(gaus1high)) + " \n")
    ResultsFile.write("gaus2low   = " + str(list(gaus2low)) + " \n")
    ResultsFile.write("gaus2high  = " + str(list(gaus2high)) + " \n")
    ResultsFile.write("min_peak   = " + str(min_peak) + "\n")
    ResultsFile.write("max_peak   = " + str(max_peak) + "\n")
    ResultsFile.write("nbins      = " + str(nbins) + "\n")
    ResultsFile.write("\n\nResults:\n")
    ResultsFile.write("GAIN = " + str(list(GAIN)) + " \n")
    ResultsFile.write("Vbd = "  + str(Vbd) + " \n")
    ResultsFile.close()

    # plt.show(block=False)
    plt.show()
    input("Press Enter to continue... ")
    plt.close()


def gaus(x, a, x0, sigma):
    return a*exp(-(x-x0)**2/(2*sigma**2))

def line(x, m, q):
    return m*x + q

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
