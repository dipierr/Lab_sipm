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

__description__ = 'Plot I-V'

formatter = argparse.ArgumentDefaultsHelpFormatter
PARSER = argparse.ArgumentParser(description=__description__, formatter_class=formatter)
PARSER.add_argument('-flist', '--input_filelist', type=str, required=False, default='file.txt', help='File to be analyzed')

def main(**kwargs):
    file = kwargs['input_filelist']
    nFile = 0
    
    T    = np.array([0, 5])
    errT = np.full(len(T), 0.5)
    Files = []
    
    min_I = 23
    max_I = 31
    
    titleI="I (mA)"
    titleV="V (V)"
     
    with open(file, "r") as infilelist:
        for f in infilelist:
            f = f.rstrip('\n')
            print("\nReading file:\n" + FindNameAndPath(f))
            Files.append(f)
            with open(f, "r") as file:
                 lines = file.read().split("\n")
                 I = 0
                 V = 0
                 nPoint = 0
                 introduction_bool = True
                
                 for l in lines:
                     if l == "END_INTRODUCTION":
                        introduction_bool = False
                        continue
                     if introduction_bool:
                        print(l)
                     else:
                        values = l.split("\t")
                     if(len(values) > 1):
                        temp = float(values[0])

                     if(temp < max_I):
                        nPoint = nPoint +1

            nFile = nFile + 1       
     
     # define length of columns
    I = ["" for x in range(nPoint)]
    V = ["" for x in range(nPoint)]
    
    j = 0
    k = 0
    
    
    introduction_bool = True

    for l in lines:
        if l == "END_INTRODUCTION":
            introduction_bool = False
            continue
        if not introduction_bool:
            values = l.split("\t")
            if(len(values) > 1):
                temp=float(values[0])
                if temp < up_lim_ALL:
                    V[j]=float(values[0])
                    I[k]=float(values[1])
# Import Plot Settings from PlotSettings.py:
    PlotSettings.PlotSettings()

    color = ['black', '#964B00', 'red', '#FFD700', 'green', 'cyan', 'blue', '#5000FF']
    
# Plot I-V at different T
    plt.figure(figsize=(10, 6))
    for i in range(len(T)):
        label = str(T[i]) + " °C"
        plt.step(V, I[i], color=color[i], linestyle='-', label=label)
    plt.xlabel(titleV)
    plt.ylabel(titleI)
    plt.legend()
    
    # plt.show(block=False)
    plt.show()
    input("Press Enter to continue... ")
    plt.close()
    
def FindNameAndPath(path):
    head, tail = ntpath.split(path)
    return tail or ntpath.basename(head)

if __name__ == '__main__':
    args = PARSER.parse_args()
    main(**args.__dict__)