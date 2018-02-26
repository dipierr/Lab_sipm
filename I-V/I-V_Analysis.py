#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++       ANALYSIS I-V CURVE for SiPM       +++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

'''
File to be analyzed: tab separated, written as:
V (V)   I (muA) errI (muA)
with NO header

Before Vdb: linear fit
After Vbd: parabolic fit


Reference:
    Nagy Ferenc & al. - A model based DC analysis of SiPM breakdown voltages

'''

import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
from scipy.optimize import curve_fit
from scipy import asarray as ar,exp
from scipy.odr import *
from pylab import *
import argparse


__description__ = 'ANALYSIS I-V CURVE for SiPM'

formatter = argparse.ArgumentDefaultsHelpFormatter
PARSER = argparse.ArgumentParser(description=__description__, formatter_class=formatter)
PARSER.add_argument('-f', '--input_file', type=str, required=False, default='file.txt', help='File to be analyzed.')



def line_fit(p,x):
    m, q = p
    return m*x + q

def parabola_fit(p,x):
    a, b, c = p
    return a*(x)**2+b*x+c
    
def main(**kwargs):
    
    low_lim = 24.2
    central_lim_low=26
    central_lim_high=26.5
    central_lim_high_2=28
    #central_lim_high_2=27.5 # non va il test z con il confronto con GAIN-V
    up_lim = 30

    intersect_low = 26.2;
    intersect_up = 27;

    void=1.5

    errV = 0.01


    # Read Text File in Python
    #input_file ="I-V_Dark"
    #input_file = "I-V_Dark_02"
    file = open(kwargs['input_file'], "r")
    lines = file.read().split("\n")
    
    print(kwargs['input_file']+ '\n')
    
    # define length of columns
    I = ["" for x in range(len(lines))]
    V = ["" for x in range(len(lines))]
    errI = ["" for x in range(len(lines))]


    #separate vectors before and after breakdown
    j=0
    before_Vbd=0 # lenght of vector V < Vbd
    after_Vbd=0 # lenght of vector V > Vbd, 1st case
    after_Vbd_2=0 # lenght of vector V > Vbd, 2nd case
    for l in lines:
        values = l.split("\t")
        V[j]=float(values[0])
        I[j]=float(values[1])
        errI[j]=float(values[2])
        
        if low_lim < V[j] < central_lim_low:
            before_Vbd=before_Vbd+1
        if central_lim_high < V[j] < up_lim:
            after_Vbd=after_Vbd+1
        if central_lim_high_2 < V[j] < up_lim:
            after_Vbd_2=after_Vbd_2+1
        j=j+1
        

    V1 = ["" for x in range(before_Vbd)]
    I1 = ["" for x in range(before_Vbd)]
    errI1 = ["" for x in range(before_Vbd)]

    V2 = ["" for x in range(after_Vbd)]
    I2 = ["" for x in range(after_Vbd)]
    errI2 = ["" for x in range(after_Vbd)]

    V2_2 = ["" for x in range(after_Vbd_2)]
    I2_2 = ["" for x in range(after_Vbd_2)]
    errI2_2 = ["" for x in range(after_Vbd_2)]

    j=0
    k=0
    h=0
    for l in lines:
        values = l.split("\t")
        temp=float(values[0])
        if low_lim < temp < central_lim_low:
            V1[j]=float(values[0])
            I1[j]=float(values[1])
            errI1[j]=float(values[2])
            j=j+1
        if central_lim_high < temp < up_lim:
            V2[k]=float(values[0])
            I2[k]=float(values[1])
            errI2[k]=float(values[2])
            k=k+1
        if central_lim_high_2 < temp < up_lim:
            V2_2[h]=float(values[0])
            I2_2[h]=float(values[1])
            errI2_2[h]=float(values[2])
            h=h+1

    #plot I-V
    plt.plot(V, I, color='blue', marker='o', linestyle='None', markersize=1)
    plt.xlabel('V (V)', fontsize = 18)
    plt.ylabel('I ('+'$\mu$'+'A)', fontsize = 18)
    plt.grid(True)
    plt.errorbar(V, I, xerr=errV, yerr=errI, linestyle='None')
    plt.ylim(10**(-3), 10**(2))
    plt.yscale('log')


    #---------
    #  ODR FIT
    #---------

    #FIT BEFORE Vbd: LINE
    print( 'BEFORE')
    fit1_model = Model(line_fit)
    data1 = RealData(V1, I1, sx=errV, sy=errI1)
    odr1 = ODR(data1, fit1_model, beta0=[1.,1.])
    out = odr1.run()
    out.pprint() 
    pcov1ODR = odr1.output.cov_beta # Covariance Matrix
    optimizedParameters1 = odr1.output.beta # fit parameters

    #FIT AFTER Vbd: PARABOLA 1
    print( '\n\nAFTER_1')
    fit2_model = Model(parabola_fit)
    data2 = RealData(V2, I2, sx=errV, sy=errI2)
    odr2 = ODR(data2, fit2_model, beta0=[1.,1.,1.])
    out = odr2.run()
    out.pprint()
    pcov2ODR = odr2.output.cov_beta # Covariance Matrix
    optimizedParameters2 = odr2.output.beta # fit parameters

    #FIT AFTER Vbd: PARABOLA 2
    print( '\n\nAFTER_2')
    fit2_2_model = Model(parabola_fit)
    data2_2 = RealData(V2_2, I2_2, sx=errV, sy=errI2_2)
    odr2_2 = ODR(data2_2, fit2_2_model, beta0=[1.,1.,1.])
    out = odr2_2.run()
    out.pprint()
    pcov2ODR_2 = odr2_2.output.cov_beta # Covariance Matrix
    optimizedParameters2_2 = odr2_2.output.beta # fit parameters

    #PLOT
    xfit1 = np.arange(low_lim, central_lim_low + void,0.01)
    xfit2 = np.arange(intersect_low,up_lim,0.01)
    yfit1 = line_fit(optimizedParameters1, xfit1)
    yfit2 =  parabola_fit(optimizedParameters2, xfit2)
    yfit2_2 =  parabola_fit(optimizedParameters2_2, xfit2)

    plt.plot(xfit1, yfit1)
    #plt.plot(xfit2, yfit2)
    plt.plot(xfit2, yfit2_2)


    #-----------------
    #   INTERSECTION 1
    #-----------------
    div_int = 0.000001
    div_int_close = 0.0000001
    xfit_int = np.arange(intersect_low, intersect_up,div_int)

    yfit1_int = line_fit(optimizedParameters1, xfit_int)
    yfit2_int =  parabola_fit(optimizedParameters2, xfit_int)

    # Print the cross point
    idx = np.argwhere(np.isclose(yfit1_int, yfit2_int, atol=div_int_close)).reshape(-1)
    print ('\nIntersection_1:\t{}'.format(xfit_int[idx]))
    print ('Before_1:\t{}'.format(xfit_int[idx-1]))
    print ('After_1:\t{}'.format(xfit_int[idx+1]))
    Intersection_1 = xfit_int[idx][0]

    #ERROR PROPAGATION
    m=optimizedParameters1[0]
    q=optimizedParameters1[1]
    a = optimizedParameters2[0]
    b = optimizedParameters2[1]
    c = optimizedParameters2[2]

    #partial derivatives
    df_dm = 1/(2*a) * ( 1- ((b-m)**2-4*a*(c-q))**(-0.5)*(b-m) )
    df_dq = 1/np.sqrt((b - m)**2 - 4 *a* (c - q))
    df_da = -(-b + m + np.sqrt((b - m)**2 - 4* a* (c - q)))/(2* a**2) - (c - q)/(a *np.sqrt((b - m)**2 - 4* a* (c - q)))
    df_db = (-1 + (b - m)/np.sqrt((b - m)**2 - 4* a *(c - q)))/(2 *a)
    df_dc = -1/np.sqrt((b - m)**2 - 4 *a* (c - q))

    #sigma
    sigma_1 = 0
    sigma_1 =         pcov2ODR[0,0] * df_da * df_da + pcov2ODR[0,1] * df_da * df_db + pcov2ODR[0,2] * df_da * df_dc
    sigma_1 = sigma_1 + pcov2ODR[1,0] * df_db * df_da + pcov2ODR[1,1] * df_db * df_db + pcov2ODR[1,2] * df_db * df_dc
    sigma_1 = sigma_1 + pcov2ODR[2,0] * df_dc * df_da + pcov2ODR[2,1] * df_dc * df_db + pcov2ODR[2,2] * df_dc * df_dc

    sigma_1 = sigma_1 + pcov1ODR[0,0] * df_dm * df_dm + pcov1ODR[0,1] * df_dm * df_dq
    sigma_1 = sigma_1 + pcov1ODR[1,0] * df_dq * df_dm + pcov1ODR[1,1] * df_dq * df_dq
    
    sigma_1 = sigma_1**0.5

    print( '\nsigma_1 = '+ str(sigma_1))

    #-----------------
    #   INTERSECTION 2
    #-----------------
    yfit2_2_int =  parabola_fit(optimizedParameters2_2, xfit_int)


    # Print the cross point
    idx = np.argwhere(np.isclose(yfit1_int, yfit2_2_int, atol=div_int_close)).reshape(-1)
    Vbd_logI_V = xfit_int[idx]
    print ('\nIntersection_2:\t{}'.format(xfit_int[idx]))
    print ('Before_2:\t{}'.format(xfit_int[idx-1]))
    print ('After_2:\t{}'.format(xfit_int[idx+1]))
    Intersection_2 = xfit_int[idx][0]

    #ERROR PROPAGATION
    m=optimizedParameters1[0]
    q=optimizedParameters1[1]
    a = optimizedParameters2_2[0]
    b = optimizedParameters2_2[1]
    c = optimizedParameters2_2[2]

    #partial derivatives
    df_dm = 1/(2*a) * ( 1- ((b-m)**2-4*a*(c-q))**(-0.5)*(b-m) )
    df_dq = 1/np.sqrt((b - m)**2 - 4 *a* (c - q))
    df_da = -(-b + m + np.sqrt((b - m)**2 - 4* a* (c - q)))/(2* a**2) - (c - q)/(a *np.sqrt((b - m)**2 - 4* a* (c - q)))
    df_db = (-1 + (b - m)/np.sqrt((b - m)**2 - 4* a *(c - q)))/(2 *a)
    df_dc = -1/np.sqrt((b - m)**2 - 4 *a* (c - q))

    #sigma
    sigma_2 = 0
    sigma_2 =         pcov2ODR_2[0,0] * df_da * df_da + pcov2ODR_2[0,1] * df_da * df_db + pcov2ODR_2[0,2] * df_da * df_dc
    sigma_2 = sigma_2 + pcov2ODR_2[1,0] * df_db * df_da + pcov2ODR_2[1,1] * df_db * df_db + pcov2ODR_2[1,2] * df_db * df_dc
    sigma_2 = sigma_2 + pcov2ODR_2[2,0] * df_dc * df_da + pcov2ODR_2[2,1] * df_dc * df_db + pcov2ODR_2[2,2] * df_dc * df_dc

    sigma_2 = sigma_2 + pcov1ODR[0,0] * df_dm * df_dm + pcov1ODR[0,1] * df_dm * df_dq
    sigma_2 = sigma_2 + pcov1ODR[1,0] * df_dq * df_dm + pcov1ODR[1,1] * df_dq * df_dq
    
    sigma_2 = sigma_2**0.5
    
    errVbd_logI_V = sigma_2
    print ('\nsigma_2 = '+str(sigma_2))



    # MEAN
    nSigma = 1
    IntersectionDown = Intersection_1 - nSigma*sigma_1
    IntersectionUp = Intersection_2 + nSigma*sigma_2
    Intersection = (IntersectionUp + IntersectionDown) / 2
    downErr = Intersection - IntersectionDown
    upErr = IntersectionUp - Intersection
    if round(upErr,2) != round(downErr,2):
        print ('ERROR, check the code')
    
    print('\n\nMEAN WAY')
    print( 'Intersection = '+str( round(Intersection,2))+ '+-'+ str(round(upErr,2)))
    
    
    print('\n\n\n')
    print('------------------------------')
    print('-------[   RESULTS   ]--------')
    print('------------------------------')
    print('\n')
    print('log(I) - V')
    print('Vdb = ('+str(Vbd_logI_V)+' +- '+str(errVbd_logI_V)+') V')

    plt.show()

if __name__ == '__main__':
    args = PARSER.parse_args()
    main(**args.__dict__)

