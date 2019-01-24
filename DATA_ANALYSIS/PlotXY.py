# PlotXY.py

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

    # Import Plot Settings from PlotSettings.py:
    PlotSettings.PlotSettings()

    titleX = "$X_{title}$"
    titleY = "$Y_{title}$"

    # Open File for XY plot
    file = open(kwargs['input_file'], "r")
    lines = file.read().split("\n")

    print("\nReading file:")
    print(kwargs['input_file'] + '\n')

    x = np.array([])
    y = np.array([])

    introduction_bool = True

    for l in lines:
        if l == "END_INTRODUCTION":
            introduction_bool = False
        if introduction_bool:
            print(l)
        else:
            temp = re.findall(r"[-+]?\d*\.\d+|\d+", l)
            if temp:
                x = np.append(x,float(temp[0]))
                y = np.append(y,float(temp[1]))

    # Plot
    plt.figure(figsize=(10, 6))
    plt.plot(x, y, color='blue', marker='o', linestyle='None', markersize=4)
    plt.xlabel(titleX)
    plt.ylabel(titleY)
    plt.grid(True)

    plt.show(block=False)
    input("Press Enter to continue... ")
    plt.close()

if __name__ == '__main__':
    args = PARSER.parse_args()
    main(**args.__dict__)
