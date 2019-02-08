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
    
    T    = np.array([-10, -5, 0, 5, 10, 15, 20, 25, 30])
    errT = np.full(len(T), 0.5)
    V = np.array([])
    I = np.array([])
    Files = []
    
    titleI="I (A)"
    titleV="V (V)"
     
    with open(file, "r") as infilelist:
        for f in infilelist:
            f = f.rstrip('\n')
            print("\nReading file:\n" + FindNameAndPath(f))
            Files.append(f)
            with open(f, "r") as file:
                lines = file.read().split("\n")

                V_temp = np.array([])
                I_temp = np.array([])

                introduction_bool = True

                for l in lines:
                    if l == "END_INTRODUCTION":
                        introduction_bool = False
                    if introduction_bool:
                        print(l)
                    else:
                        temp = l.split("\t")
                        if len(temp) > 1:
                            if float(temp[0])>23 and float(temp[0])<31:
                                V_temp = np.append(V_temp,float(temp[0]))
                                I_temp = np.append(I_temp,float(temp[1]))
                
                if nFile == 0:
                    V = V_temp
                    I = I_temp
                else:
                    V = np.vstack([V, V_temp])
                    I = np.vstack([I, I_temp])
                
                nFile = nFile + 1
     
    
# Import Plot Settings from PlotSettings.py:
    PlotSettings.PlotSettings()

    color = ['#c6ecd7', '#9fdfbc', '#79d2a1', '#53c687', '#3cb371', '#339961', '#267349', '#194d30', '#0d2618']
    
# Plot I-V at different T
    plt.figure(figsize=(10, 6))
    for i in range(len(T)):
        label = str(T[i]) + " °C"
        plt.plot(V[i], I[i], color=color[i], linestyle='-', label=label)
    #plt.yscale("log")
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