# Plot_value.py

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

def line_fit(p,x):
    m, q = p
    return m*x + q

def main(**kwargs):

    # Import Plot Settings from PlotSettings.py:
    PlotSettings.PlotSettings()

    titleT = "$T$ (°C)"
    titleVbd_retta_parabola = "$V$ (V)"
    titleVbd_sqrt = "$V$ (V)"
    
    # Open File for XY plot
    file = open(kwargs['input_file'], "r")
    lines = file.read().split("\n")

    print("\nReading file:")
    print(kwargs['input_file'] + '\n')

    T = np.array([])
    Vbd_retta_parabola = np.array([])
    Vbd_sqrt = np.array([])
    errVbd_retta_parabola = np.array([])
    errVbd_sqrt = np.array([])

    introduction_bool = True

    for l in lines:
        if l == "END_INTRODUCTION":
            introduction_bool = False
        if introduction_bool:
            print(l)
        else:
            temp = l.split("\t")
            if len(temp) > 1:
                T = np.append(T,float(temp[0]))
                Vbd_retta_parabola = np.append(Vbd_retta_parabola,float(temp[1])-float(temp[2]))
                Vbd_sqrt = np.append(Vbd_sqrt,float(temp[3]))
                errVbd_retta_parabola = np.append(errVbd_retta_parabola,float(temp[2]))
                errVbd_sqrt = np.append(errVbd_sqrt,float(temp[4]))

    # Plot retta parabola
    plt.figure(figsize=(10, 6))
    plt.errorbar(T, Vbd_retta_parabola, yerr=errVbd_retta_parabola, xerr=0.5, color='red', marker='o', linestyle='None', markersize=4)
    plt.xlabel(titleT)
    plt.ylabel(titleVbd_retta_parabola)
    plt.grid(True)
    axes = plt.gca()
    axes.set_xlim([-11,31])
    # axes.set_ylim([ymin,ymax])
    
     # FIT LINEARE RETTA PARABOLA 
    fit1_model = Model(line_fit)
    data1 = RealData(T, Vbd_retta_parabola, sx=0.5, sy=errVbd_retta_parabola)
    odr1 = ODR(data1, fit1_model, beta0=[1.,1.])
    out = odr1.run()
    out.pprint()
    pcov1ODR = odr1.output.cov_beta # Covariance Matrix
    optimizedParameters1 = odr1.output.beta # fit parameters
    
    #PLOT
    xfit1 = np.arange(-10,30,0.01)
    yfit1 = line_fit(optimizedParameters1, xfit1)
    
    plt.plot(xfit1, yfit1, color = 'black')
    
     # Plot sqrt
    plt.figure(figsize=(10, 6))
    plt.errorbar(T, Vbd_sqrt, xerr=0.5, yerr=errVbd_sqrt, color='blue', marker='o', linestyle='None', markersize=4)
    plt.xlabel(titleT)
    plt.ylabel(titleVbd_sqrt)
    plt.grid(True)
    axes = plt.gca()
    axes.set_xlim([-11,31])
    # axes.set_ylim([ymin,ymax])
        
     # FIT LINEARE SQRT 
    fit2_model = Model(line_fit)
    data2 = RealData(T, Vbd_sqrt, sx=0.5, sy=errVbd_sqrt)
    odr1 = ODR(data2, fit2_model, beta0=[1.,1.])
    out = odr1.run()
    out.pprint()
    pcov1ODR = odr1.output.cov_beta # Covariance Matrix
    optimizedParameters1 = odr1.output.beta # fit parameters
    
    #PLOT
    xfit2 = np.arange(-11,30,0.01)
    yfit2 = line_fit(optimizedParameters1, xfit2)
    
    plt.plot(xfit2, yfit2, color = 'black')
    #plt.plot(xfit2, yfit2)
    
    # Plot Vbd retta parabole e Vbd sqrt
    #plt.figure(figsize=(10, 6))
    #lines=plt.errorbar(T, Vbd_retta_parabola, xerr=0.5, T, Vbd_sqrt)
    #l1,l2=lines
    #plt.setp(l1, color='red',marker='o', linestyle='None')
    #plt.setp(l2,color='blue',marker='o', linestyle='None')
    #plt.yscale("log")
    #plt.xlabel(titleT)
    #plt.ylabel('Vbd')
    # plt.legend()
    
    Vbd = [Vbd_retta_parabola,Vbd_sqrt]
    errorVbd = [errVbd_retta_parabola, errVbd_sqrt]
    label = ['Vbd,1', 'Vbd,2']
    
    color = ['red', 'blue']
    
    plt.figure(figsize=(10, 6))
    for i in range(len(Vbd)):
        # label = Vbds[i]
        plt.errorbar(T, Vbd[i], xerr=0.5, yerr= errorVbd[i], color=color[i], linestyle='None', marker='o', label=label[i])
    #plt.yscale("log")
    plt.xlabel(titleT)
    plt.ylabel('V (V)')
    plt.legend()

    plt.show(block=False)
    input("Press Enter to continue... ")
    plt.close()

if __name__ == '__main__':
    args = PARSER.parse_args()
    main(**args.__dict__)
