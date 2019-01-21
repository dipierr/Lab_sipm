# I-V_Analysis.py


'''
---------------------------------------------------------------------------------
|                                                                               |
|   I-V_Analysis.py                                                             |
|                                                                               |
|   File to be analyzed: tab separated, written as:                             |
|   V (V)   I (A)                                                  |
|   with NO header                                                              |
|                                                                               |
|   Before Vdb: linear fit                                                      |
|   After Vbd: parabolic fit                                                    |
|                                                                               |
|   Reference:                                                                  |
|   Nagy Ferenc & al. - A model based DC analysis of SiPM breakdown voltages    |
|                                                                               |
|   Then compared to Sqrt(I) intersection                                       |
|                                                                               |
|   Davide Depaoli 2018 - 2019                                                  |
|                                                                               |
|                                                                               |
---------------------------------------------------------------------------------
'''



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



__description__ = 'ANALYSIS I-V CURVE for SiPM'

formatter = argparse.ArgumentDefaultsHelpFormatter
PARSER = argparse.ArgumentParser(description=__description__, formatter_class=formatter)
PARSER.add_argument('-f', '--input_file', type=str, required=False, default='file.txt', help='File to be analyzed.')
PARSER.add_argument('-NoFit', '--no_fit_gain', action='store_true', required=False, default=False, help='Disables histogram fitting')




def line_fit(p,x):
    m, q = p
    return m*x + q

def parabola_fit(p,x):
    a, b, c = p
    return a*(x)**2+b*x+c


def main(**kwargs):

    PlotSettings.PlotSettings()

#-------------------------------------------------------------------------------
#--------------------------------[   LOG(I) - V   ]-----------------------------
#-------------------------------------------------------------------------------

    low_lim = 24.2
    central_lim_low=26
    central_lim_high=26.6 #26.6
    #central_lim_high=27.5 # non va il test z con il confronto con GAIN-V
    up_lim = 30

    intersect_low = 26.2
    intersect_up = 40 #27;

    void=1.5

    errV = 0.01

    approx_Vbd_sqrt = 26.6 #26.6

    up_lim_ALL = 30



    file = open(kwargs['input_file'], "r")
    lines = file.read().split("\n")

    print(kwargs['input_file']+ '\n')


    #separate vectors before and after breakdown
    j=0
    before_Vbd=0 # lenght of vector V < Vbd
    after_Vbd=0 # lenght of vector V > Vbd, 1st case
    after_Vbd=0 # lenght of vector V > Vbd, 2nd case
    before_Vbd_sqrt=0 # sqrt case
    ALL = 0

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

                if(temp < up_lim_ALL):
                    ALL = ALL +1

                if low_lim < temp < central_lim_low:
                    before_Vbd=before_Vbd+1

                if central_lim_high < temp < up_lim:
                    after_Vbd=after_Vbd+1
                if temp < approx_Vbd_sqrt:
                   before_Vbd_sqrt = before_Vbd_sqrt+1
                j=j+1


    # define length of columns
    I = ["" for x in range(ALL)]
    V = ["" for x in range(ALL)]
    errI = ["" for x in range(ALL)]


    V1 = ["" for x in range(before_Vbd)]
    I1 = ["" for x in range(before_Vbd)]
    errI1 = ["" for x in range(before_Vbd)]

    V2 = ["" for x in range(after_Vbd)]
    I2 = ["" for x in range(after_Vbd)]
    errI2 = ["" for x in range(after_Vbd)]

    V2 = ["" for x in range(after_Vbd)]
    I2 = ["" for x in range(after_Vbd)]
    errI2 = ["" for x in range(after_Vbd)]

    ind_all = 0
    j=0
    k=0
    h=0

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
                    V[ind_all]=float(values[0])
                    I[ind_all]=float(values[1])
                    errI[ind_all]=float(0)
                    # errI[ind_all]=float(values[2])
                    ind_all=ind_all+1
                if low_lim < temp < central_lim_low:
                    V1[j]=float(values[0])
                    I1[j]=float(values[1])
                    errI1[j]=float(0)
                    # errI1[j]=float(values[2])
                    j=j+1
                if central_lim_high < temp < up_lim:
                    V2[h]=float(values[0])
                    I2[h]=float(values[1])
                    errI2[h]=float(0)
                    # errI2[h]=float(values[2])
                    h=h+1

    #plot I-V
    plt.figure(figsize=(10, 6))
    plt.plot(V, I, color='blue', marker='o', linestyle='None', markersize=1)
    plt.xlabel(r'$V (V)$')
    plt.ylabel(r'$I (A)$')
    plt.grid(True)
    plt.errorbar(V, I, xerr=errV, yerr=errI, linestyle='None')
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    # plt.ylim(10**(-3), 10**(1))
    # plt.yscale('log')


    if not kwargs["no_fit_gain"]:
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


        #FIT AFTER Vbd: PARABOLA 2
        print( '\n\nAFTER')
        fit2_model = Model(parabola_fit)
        data2 = RealData(V2, I2, sx=errV, sy=errI2)
        odr2 = ODR(data2, fit2_model, beta0=[1.,1.,1.])
        out = odr2.run()
        out.pprint()
        pcov2ODR = odr2.output.cov_beta # Covariance Matrix
        optimizedParameters2 = odr2.output.beta # fit parameters

        #PLOT
        xfit1 = np.arange(low_lim, central_lim_low + void,0.01)
        xfit2 = np.arange(intersect_low,up_lim,0.01)
        yfit1 = line_fit(optimizedParameters1, xfit1)
        yfit2 =  parabola_fit(optimizedParameters2, xfit2)

        plt.plot(xfit1, yfit1)
        #plt.plot(xfit2, yfit2)
        plt.plot(xfit2, yfit2)



        #------------------
        #   INTERSECTION
        #------------------
        div_int = 0.000002
        div_int_close = 0.0000001
        xfit_int = np.arange(intersect_low, intersect_up,div_int)
        yfit1_int = line_fit(optimizedParameters1, xfit_int)
        yfit2_int =  parabola_fit(optimizedParameters2, xfit_int)


        # Print the cross point
        idx = np.argwhere(np.isclose(yfit1_int, yfit2_int, atol=div_int_close)).reshape(-1)
        Vbd_logI_V = xfit_int[idx]
        print ('\nIntersection:\t{}'.format(xfit_int[idx]))
        print ('Before:\t{}'.format(xfit_int[idx-1]))
        print ('After:\t{}'.format(xfit_int[idx+1]))
        Intersection = xfit_int[idx][0]

        plt.plot(xfit_int[idx], yfit1_int[idx], color='red', marker='o', linestyle='None', markersize=5)

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
        sigma = 0
        sigma =         pcov2ODR[0,0] * df_da * df_da + pcov2ODR[0,1] * df_da * df_db + pcov2ODR[0,2] * df_da * df_dc
        sigma = sigma + pcov2ODR[1,0] * df_db * df_da + pcov2ODR[1,1] * df_db * df_db + pcov2ODR[1,2] * df_db * df_dc
        sigma = sigma + pcov2ODR[2,0] * df_dc * df_da + pcov2ODR[2,1] * df_dc * df_db + pcov2ODR[2,2] * df_dc * df_dc

        sigma = sigma + pcov1ODR[0,0] * df_dm * df_dm + pcov1ODR[0,1] * df_dm * df_dq
        sigma = sigma + pcov1ODR[1,0] * df_dq * df_dm + pcov1ODR[1,1] * df_dq * df_dq

        sigma = sigma**0.5

        errVbd_logI_V = sigma
        print ('\nsigma = '+str(sigma))


    #-------------------------------------------------------------------------------
    #------------------------------[   SQRT(I) - V   ]------------------------------
    #-------------------------------------------------------------------------------

        # All the plot
        SqrtI = np.sqrt(I)
        errSqrtI = errI/(2*SqrtI)

        # Select after Vbd:
        up_ind = before_Vbd_sqrt + 20;
        V_after = V[before_Vbd_sqrt:]
        I_after = I[before_Vbd_sqrt:]
        errI_after = errI[before_Vbd_sqrt:]

        SqrtI_after = np.sqrt(I_after)
        errSqrtI_after = errI_after/(2*SqrtI_after)

        #plot SqrtI-V
        plt.figure(figsize=(10, 6))
        plt.plot(V, SqrtI, color='blue', marker='o', linestyle='None', markersize=1)
        # plt.plot(V_after, SqrtI_after, color='red', marker='o', linestyle='None', markersize=2)
        plt.xlabel(r'$V (V)$')
        plt.ylabel(r'$\sqrt{I} (\sqrt{A})$')
        plt.grid(True)
        plt.errorbar(V, SqrtI, xerr=errV, yerr=errSqrtI, linestyle='None')

        #FIT
        print( 'SQRT(I) - V')
        fit_model_sqrt = Model(line_fit)
        data_sqrt = RealData(V_after, SqrtI_after, sx=errV, sy=errSqrtI_after)
        odr_sqrt = ODR(data_sqrt, fit_model_sqrt, beta0=[1.,1.])
        out = odr_sqrt.run()
        out.pprint()
        pcov_sqrt = odr_sqrt.output.cov_beta # Covariance Matrix
        optimizedParameters_sqrt = odr_sqrt.output.beta # fit parameters

        xfit_sqrt = np.arange(25, up_lim_ALL,0.01)
        plt.plot(xfit_sqrt,line_fit(optimizedParameters_sqrt, xfit_sqrt) )



        m=optimizedParameters_sqrt[0]
        q=optimizedParameters_sqrt[1]

        Vbd_sqrtI_V = -q/m


        plt.plot(Vbd_sqrtI_V, 0, color='red', marker='o', linestyle='None', markersize=5)


        #ERROR PROPAGATION


        #partial derivatives
        df_dm = q/(m**2)
        df_dq = 1/m

        #sigma
        sigma_sqrt = 0
        sigma_sqrt = sigma_sqrt + pcov1ODR[0,0] * df_dm * df_dm + pcov1ODR[0,1] * df_dm * df_dq
        sigma_sqrt = sigma_sqrt + pcov1ODR[1,0] * df_dq * df_dm + pcov1ODR[1,1] * df_dq * df_dq

        sigma_sqrt = sigma_sqrt**0.5

        errVbd_sqrtI_V = sigma_sqrt



        print('\n\n\n')
        print('------------------------------')
        print('-------[   RESULTS   ]--------')
        print('------------------------------')
        print('\n')
        # I only consider the second intersection:
        print('log(I) - V')
        print('\tVdb = ('+str(Vbd_logI_V)+' +- '+str(errVbd_logI_V)+') V')
        print('sqrt(I) - V')
        print('\tVdb = ('+str(Vbd_sqrtI_V)+' +- '+str(errVbd_sqrtI_V)+') V')



    # plt.show(block=False)
    plt.show()
    input("Press Enter to continue... ")
    plt.close()

if __name__ == '__main__':
    args = PARSER.parse_args()
    main(**args.__dict__)
