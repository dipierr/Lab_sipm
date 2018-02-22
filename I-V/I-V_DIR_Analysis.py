#I-V_DIR_Analysis.py

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++       ANALYSIS I-V CURVE for SiPM (DIRECT)      +++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

'''
File to be analyzed: tab separated, written as:
V (V)	I (muA) errI (muA)
with NO header
'''

import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
from scipy.optimize import curve_fit
from scipy import asarray as ar,exp
import argparse


__description__ = 'ANALYSIS I-V CURVE for SiPM (DIRECT)'

formatter = argparse.ArgumentDefaultsHelpFormatter
PARSER = argparse.ArgumentParser(description=__description__, formatter_class=formatter)
PARSER.add_argument('-f', '--input_file', type=str, required=False, default='file.txt', help='File to be analyzed.')


def general_exp(x, a, b,c, x0):
    return c*(np.exp(a*(x-x0))) + b

def general_line(x,m,q):
    return m*x + q

def main(**kwargs):

    low_lim = 0.64
    up_lim = 0.96


    # Read Text File in Python
    #input_file ="I-V_Dark_Dir"
    file = open(kwargs['input_file'], "r")
    lines = file.read().split("\n")
       
    print kwargs['input_file'], '\n'
    
    # define length of columns
    I = ["" for x in range(len(lines))]
    V = ["" for x in range(len(lines))]
    errI = ["" for x in range(len(lines))]


    #separate vectors before and after breakdown
    j=0
    index=0
    after_Vbd=0
    for l in lines:
        values = l.split("\t")
        V[j]=float(values[0])
        I[j]=float(values[1])*10**(-6)
        errI[j]=float(values[2])*10**(-6)
        if low_lim < V[j] < up_lim:
            index=index+1
        j=j+1
        

    V1 = np.linspace(0,1,index)
    I1 = np.linspace(0,1,index)
    errI1 = np.linspace(0,1,index)

    j=0
    k=0
    for l in lines:
        values = l.split("\t")
        temp=float(values[0])
        if low_lim < temp < up_lim:
            V1[j]=float(values[0])
            I1[j]=float(values[1])*10**(-6) 
            j=j+1

    #plot I-V
    plt.plot(V, I, color='blue', marker='o', linestyle='None', markersize=1)
    plt.xlabel('V (V)')
    plt.ylabel('I ('+'A)')
    plt.grid(True)
    plt.errorbar(V, I, xerr=0.01, yerr=errI, linestyle='None')

    #plt.yscale('log')


    #fit
    optimizedParameters1, pcov1 = curve_fit(general_line, V1, I1)
    xfit1 = np.arange(low_lim, up_lim, 0.001)
    yfit1 = general_line(xfit1, *optimizedParameters1)
    plt.plot(xfit1, yfit1, label="fit_after")
    m = optimizedParameters1[0]
    q = optimizedParameters1[1]
    errm = pcov1[0][0]**(0.5)
    Req = optimizedParameters1[0]**(-1)
    errReq= errm/m*Req
    Rq=22500*Req
    errRq=errReq/Req*Rq
    print('[m, q] = {}'.format(optimizedParameters1))
    print('errm = {}'.format(errm))
    print('Req = {}'.format(Req))
    print('errReq = {}'.format(errReq))
    print('Rq = {}'.format(Rq))
    print('errRq = {}'.format(errRq))

    print stats.chisquare(I1, f_exp = general_line(V1, *optimizedParameters1))

    plt.show()
    
if __name__ == '__main__':
	args = PARSER.parse_args()
	main(**args.__dict__)
