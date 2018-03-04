#GAIN-V_SiPM.py

#-----------
#GAIN-V_SiPM
#-----------


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
    
def fit(HV_SiPM, GAIN_SiPM, err_HV_SiPM, err_GAIN_SiPM):
    plt.figure()
    plt.xlabel('HV (V)', fontsize = 18)
    plt.ylabel('GAIN (mV)', fontsize = 18)
    plt.grid(True)
    plt.errorbar(HV_SiPM, GAIN_SiPM, xerr=err_HV_SiPM, yerr=err_GAIN_SiPM, linestyle='None')
    fit_model = Model(line_fit)
    data = RealData(HV_SiPM, GAIN_SiPM, sx=err_HV_SiPM, sy=err_GAIN_SiPM)
    odr = ODR(data, fit_model, beta0=[1.,1.])
    out = odr.run()
    out.pprint() 
    pcovODR = odr.output.cov_beta # Covariance Matrix
    optimizedParameters = odr.output.beta # fit parameters
    xfit = np.arange(25, 37,0.01)
    yfit = line_fit(optimizedParameters, xfit)
    plt.plot(xfit, yfit)
    plt.grid(True)
    return optimizedParameters, pcovODR

def error_propagation_lin_fit(optimizedParameters, pcovODR, nSiPM):
    #linear fit: 
    #y = mx +q => @y=0: mx=-q => x=-q/m
    m=optimizedParameters[0]
    q=optimizedParameters[1]
    
    Vbd = -q/m
    
    #partial derivatives
    df_dm = q/(m**2)
    df_dq = -1/m
    
    #sigma
    sigma = 0
    sigma = sigma + pcovODR[0,0]*df_dm*df_dm + pcovODR[0,1]*df_dm*df_dq
    sigma = sigma + pcovODR[1,0]*df_dq*df_dm + pcovODR[1,1]*df_dq*df_dq
    
    sigma = sigma**0.5
    
    print('---------------------------------------------------------')
    print('Vbd_SiPM'+str(nSiPM)+'\t'+str(Vbd)+' +- '+str(sigma))
    print('---------------------------------------------------------\n\n')
    

def main(**kwargs):
    
    #----------------------------------------------------
    #---------------[   SiPM 1 (HD3_2)   ]---------------
    #----------------------------------------------------
    nSiPM = 1
    HV_SiPM1 = np.array([34.00, 35.00, 36.00])
    err_HV_SiPM1 = np.array([0.01, 0.01, 0.01])
    
    #from Waweform with 2000 points (at max sample rate) (20180212_1_DARK_34_02.txt, 20180212_1_DARK_35_01.txt, 20180212_1_DARK_36_01.txt):
    #GAIN_SiPM1 = np.array([0.012197229689097585, 0.014442865747665316, 0.01599250638039572]) 
    #err_GAIN_SiPM1 = np.array([0.00019212343238446072, 0.00024665495118529894, 7.956049599216814e-05])
    
    #from files 20180221_HD3-2_1_DARK_34_AS_2_01.txt, 20180221_HD3-2_1_DARK_35_AS_2_01.txt, 20180221_HD3-2_1_DARK_36_AS_2_01.txt considering 15000 windows of 1 musec each:
    GAIN_SiPM1 = np.array([17.1506, 19.6551, 21.9663]) 
    err_GAIN_SiPM1 = np.array([0.0209611, 0.0217898, 0.0217723])
    
    optimizedParameters1, pcovODR1 = fit(HV_SiPM1, GAIN_SiPM1, err_HV_SiPM1, err_GAIN_SiPM1)
    error_propagation_lin_fit(optimizedParameters1, pcovODR1, nSiPM)
    
    
    
    
    #----------------------------------------------------
    #---------------[   SiPM 2 (HD3_2)   ]---------------
    #----------------------------------------------------
    nSiPM = 2
    HV_SiPM2 = np.array([34.00, 35.00, 36.00])
    err_HV_SiPM2 = np.array([0.01, 0.01, 0.01])
    
    #from Waweform with 2000 points (at max sample rate) (20180212_2_DARK_34_01.txt, 20180212_2_DARK_35_01.txt, 20180212_2_DARK_36_01.txt):
    #GAIN_SiPM2 = np.array([0.011410721379920264,0.0135010146806218,0.0154419724055663]) 
    #err_GAIN_SiPM2 = np.array([0.0003764632506002512,0.00037003304446492515,0.00014876147139717083])
    
    #from files 20180221_HD3-2_2_DARK_34_AS_2_02.txt, 20180221_HD3-2_2_DARK_35_AS_2_02.txt, 20180221_HD3-2_2_DARK_36_AS_2_02.txt considering 15000 windows of 1 musec each:
    GAIN_SiPM2 = np.array([16.1435 ,18.5611,20.9333]) 
    err_GAIN_SiPM2 = np.array([0.0226831,0.0241206,0.0272634])
    
    optimizedParameters2, pcovODR2 = fit(HV_SiPM2, GAIN_SiPM2, err_HV_SiPM2, err_GAIN_SiPM2)
    error_propagation_lin_fit(optimizedParameters2, pcovODR2, nSiPM)
    
    
    
    #----------------------------------------------------
    #---------------[   SiPM 3 (HD3_2)   ]---------------
    #----------------------------------------------------
    nSiPM = 3
    HV_SiPM3 = np.array([34.00, 35.00, 36.00])
    err_HV_SiPM3 = np.array([0.01, 0.01, 0.01])
    
    #from Waweform with 2000 points (at max sample rate) (20180209_DARK_34_03.txt, 20180209_DARK_35_02.txt, 20180209_DARK_36_01.txt):
    #GAIN_SiPM3 = np.array([0.015223774484941218, 0.01757637806904217, 0.019537323488766054]) 
    #err_GAIN_SiPM3 = np.array([0.0002636386386728474, 0.00031677656441436237, 0.00019640674762372849])
    
    #from files 20180221_HD3-2_1_DARK_34_AS_2_03.txt, 20180221_HD3-2_1_DARK_35_AS_2_03.txt, 20180221_HD3-2_1_DARK_36_AS_2_03.txt considering 15000 windows of 1 musec each:
    GAIN_SiPM3 = np.array([17.5343,20.0972,22.6212]) 
    err_GAIN_SiPM3 = np.array([ 0.0201935,0.022037,0.0237572])
    
    optimizedParameters3, pcovODR3 = fit(HV_SiPM3, GAIN_SiPM3, err_HV_SiPM3, err_GAIN_SiPM3)
    error_propagation_lin_fit(optimizedParameters3, pcovODR3, nSiPM)
    plt.show()


if __name__ == '__main__':
    args = PARSER.parse_args()
    main(**args.__dict__)
