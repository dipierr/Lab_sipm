#GAIN-V_SiPM


import numpy as np
import matplotlib.pyplot as plt
from scipy.odr import *
import argparse


__description__ = 'Producing the azimuth and zenith angle for sim_telarray'

formatter = argparse.ArgumentDefaultsHelpFormatter
PARSER = argparse.ArgumentParser(description=__description__, formatter_class=formatter)
PARSER.add_argument('-o', '--option', type=int, required=False, default=1, help='Option')


def line_fit(p,x):
    m, q = p
    return m*x + q


def main(**kwargs):
    
    #------
    #SiPM 1
    #------
    plt.figure()
    HV_SiPM1 = np.array([34.00, 35.00, 36.00])
    err_HV_SiPM1 = np.array([0.01, 0.01, 0.01])
    GAIN_SiPM1 = np.array([0.0151408368519, 0.0171787459721, 0.0191618624342]) #from Waweform with 2000 points (at max sample rate) (20180209_DARK_34_03.txt, 20180209_DARK_35_02.txt, 20180209_DARK_36_01.txt)
    err_GAIN_SiPM1 = np.array([3.07424563391e-07, 1.13123828024e-06, 2.82262560196e-07]) #from Waweform with 2000 points (at max sample rate) (20180209_DARK_34_03.txt, 20180209_DARK_35_02.txt, 20180209_DARK_36_01.txt)
    plt.errorbar(HV_SiPM1, GAIN_SiPM1, xerr=err_HV_SiPM1, yerr=err_GAIN_SiPM1, linestyle='None')
    
    fit1_model = Model(line_fit)
    data1 = RealData(HV_SiPM1, GAIN_SiPM1, sx=err_HV_SiPM1, sy=err_GAIN_SiPM1)
    odr1 = ODR(data1, fit1_model, beta0=[1.,1.])
    out = odr1.run()
    out.pprint() 
    pcov1ODR = odr1.output.cov_beta # Covariance Matrix
    optimizedParameters1 = odr1.output.beta # fit parameters
    xfit1 = np.arange(25, 37,0.01)
    yfit1 = line_fit(optimizedParameters1, xfit1)
    plt.plot(xfit1, yfit1)
    plt.grid(True)
    
    print('Vbd_SiPM1\t'+str(round(-optimizedParameters1[1]/optimizedParameters1[0],2)))
        
    plt.show()


if __name__ == '__main__':
    args = PARSER.parse_args()
    main(**args.__dict__)
