/******************************************************************************\
 * DCR_CrossTalk_FBK_HD3_2_LASER_data_2018_07
 *
 * Setup LED:
 *    > HV:         AGILENT E3641A
 *    > SUPPLY:     Kenwood Regulated DC Power Supply PW18-1.8Q
 *    > AMPLIFIER:  ADVANSID OUT 2
 *    > DIGITIZER:  DRS4 Evaluation Board
 *    > LASER:      Advanced Laser Diode System, Picosecond Laser System,
 *                  Controller EIG2000DX and a 406 nm laser head
 *
 * > File obtained using Ana_LED(...) function in Ana_Traces_SiPM.cxx using
 *   fit_hist_peaks_gaus_sum_012(...):
 *      [0]*TMath::Exp( - (x-[6])*(x-[6])/( 2*[4]*[4] ) ) +
 *      + [1]*TMath::Exp(-(x-[6]-[3])*(x-[6]-[3])/( 2*([4]*[4] + [5]*[5] )) ) +
 *      + [2]*TMath::Exp(-(x-[6]-2*[3])*(x-[6]-2*[3])/(2*([4]*[4] + 4*[5]*[5])))
 *
 * > Ana_Traces_SiPM.cxx github: 11/08/2018 # 1
 *
 * > Fit range and the initial values of the parameters are set manually for
 *   each file
 *
 * File used:
 *    20180614_HD3-2_1_LASER_PLS_81_PAPER_AGILENT_29_AS_2_50000_01.dat
 *    20180614_HD3-2_1_LASER_PLS_81_PAPER_AGILENT_30_AS_2_50000_01.dat
 *    20180614_HD3-2_1_LASER_PLS_81_PAPER_AGILENT_31_AS_2_50000_01.dat
 *    20180614_HD3-2_1_LASER_PLS_81_PAPER_AGILENT_32_AS_2_50000_01.dat
 *    20180614_HD3-2_1_LASER_PLS_81_PAPER_AGILENT_33_AS_2_50000_01.dat
 *    20180614_HD3-2_1_LASER_PLS_81_PAPER_AGILENT_34_AS_2_50000_01.dat
 *    20180614_HD3-2_1_LASER_PLS_81_PAPER_AGILENT_35_AS_2_50000_01.dat
 *    20180614_HD3-2_1_LASER_PLS_81_PAPER_AGILENT_36_AS_2_50000_01.dat
 *
\******************************************************************************/



#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <cstdio>
#include <fcntl.h>
#include <unistd.h>
#include <math.h>
#include <algorithm>
#include <functional>
#include <numeric>
#include "TMath.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TAxis.h"
#include "TPad.h"
#include "TLine.h"
#include "TSpectrum.h"
#include "TVirtualFitter.h"
#include "TFile.h"
#include "TTree.h"
#include "TSystem.h"
#include "TTask.h"
#include "TVirtualGraphPainter.h"
#include "TGaxis.h"
#include "TROOT.h"
#include "TString.h"

#define n_DCR_1 7
#define n_DCR_2 7
#define n_DCR_3 7
#define n_GAIN 6
#define n_MEAN 8

#define n_SiPM_tot 3



#define h 600
#define w 1000

TGraphErrors *gV_CT_1;
TGraphErrors *gV_CT_2;
TGraphErrors *gV_CT_3;
TGraphErrors *gV_CT_Del_1;
TGraphErrors *gV_CT_Del_2;
TGraphErrors *gV_CT_Del_3;
TGraphErrors *gV_MS;


int n_SiPM;

char title_DCR[80];
char title_DCR_mg[80];

void DCR_CrossTalk_FBK_HD3_2_LASER_data_2018_07();
// void DCR_CrossTalk_FBK_HD3_2_from_cnt_data_2018_07();
int find_index(double v[],int N, double value);
// void OperationPoint_LASER_PLS_FitGausSum_01();

void DCR_CrossTalk_FBK_HD3_2_LASER_data_2018_07(){

    // ERRORS
    bool fix_error_bool = true;

    // HV
    double HV_LASER[] =  {29.00,30.00,    31.00,        32.00,        33.00,        34.00,               35.00,               36.00};
    double errHV_LASER[]={0.01, 0.01,     0.01,         0.01,         0.01,         0.01,                0.01,                0.01};
    // double HV_LASER[] =  {29.00,30.00,    31.00,        32.00,        33.00,        34.00,        34.50,        35.00,        35.50,        36.00};
    // double errHV_LASER[]={0.01, 0.01,     0.01,         0.01,         0.01,         0.01,         0.01,         0.01,         0.01,         0.01};

    int n_SiPM_tot_LASER = 2; // MUST BE < n_SiPM_tot

    if(n_SiPM_tot_LASER>n_SiPM_tot){
        n_SiPM_tot_LASER = n_SiPM_tot;
    }

    // PEAK 0
    // H_peak_0[n_SiPM]
    double H_peak_0[n_SiPM_tot][n_GAIN] = {0.};
    double errH_peak_0[n_SiPM_tot][n_GAIN] = {0.};

    // Sigma_peak_0[n_SiPM_tot]
    double Sigma_peak_0[n_SiPM_tot][n_GAIN] = {0.};
    double errSigma_peak_0[n_SiPM_tot][n_GAIN] = {0.};

    // PEAK 1
    // H_peak_1[n_SiPM_tot]
    double H_peak_1[n_SiPM_tot][n_GAIN] = {0.};
    double errH_peak_1[n_SiPM_tot][n_GAIN] = {0.};

    // Sigma_peak_1[n_SiPM_tot]
    double Sigma_peak_1[n_SiPM_tot][n_GAIN] = {0.};
    double errSigma_peak_1[n_SiPM_tot][n_GAIN] = {0.};

    // Mean_peak_1[n_SiPM_tot]
    double Mean_peak_1[n_SiPM_tot][n_GAIN] = {0.};
    double errMean_peak_1[n_SiPM_tot][n_GAIN] = {0.};

    // PEAK 2
    // Mean_peak_2[n_SiPM_tot]
    double Mean_peak_2[n_SiPM_tot][n_GAIN] = {0.};
    double errMean_peak_2[n_SiPM_tot][n_GAIN] = {0.};

    // GLOBAL HIST
    // Mean hist global
    double Mean_hg[n_SiPM_tot][n_MEAN] = {0.};
    double errMean_hg[n_SiPM_tot][n_MEAN] = {0.};

    // Standard_dev_hist_global
    double Std_hg[n_SiPM_tot][n_MEAN] = {0.};
    double errStd_hg[n_SiPM_tot][n_MEAN] = {0.};

     // Integral of the average
    double Integral[n_SiPM_tot][n_GAIN] = {0.};
    double errIntegral[n_SiPM_tot][n_GAIN] = {0.};

    // Entries
    double Entries[n_SiPM_tot][n_GAIN] = {0.};

    // GAIN
    double GAIN[n_SiPM_tot][n_GAIN] = {0.};
    double errGAIN[n_SiPM_tot][n_GAIN] = {0.};

    // Other variables
    double GainWeighted[n_SiPM_tot][n_GAIN] = {0.};
    double errGainWeighted[n_SiPM_tot][n_GAIN] = {0.};
    double Mean_hist_St_dev[n_SiPM_tot][n_MEAN] = {0.};
    double errMean_hist_St_dev[n_SiPM_tot][n_MEAN] = {0.};
    double Integral_GAIN[n_SiPM_tot][n_GAIN] = {0.};
    double errIntegral_GAIN[n_SiPM_tot][n_GAIN] = {0.};
    double HV_LASER_GAIN[n_GAIN] = {0.};
    double errHV_LASER_GAIN[n_GAIN] = {0.};
    double Area0[n_SiPM_tot][n_GAIN] = {0.};
    double Area1[n_GAIN] = {0.};
    double errArea0[n_SiPM_tot][n_GAIN] = {0.};
    double errArea1[n_SiPM_tot][n_GAIN] = {0.};
    double Prob_0pe[n_SiPM_tot][n_GAIN] = {0.};
    double errProb_0pe[n_SiPM_tot][n_GAIN] = {0.};
    double Prob_1pe[n_SiPM_tot][n_GAIN] = {0.};
    double errProb_1pe[n_SiPM_tot][n_GAIN] = {0.};
    double Prob_1peS[n_SiPM_tot][n_GAIN] = {0.};
    double errProb_1peS[n_SiPM_tot][n_GAIN] = {0.};
    double Mu[n_SiPM_tot][n_GAIN] = {0.};
    double errMu[n_SiPM_tot][n_GAIN] = {0.};
    double Prob_Cross_Talk[n_SiPM_tot][n_GAIN] = {0.};
    double errProb_Cross_Talk[n_SiPM_tot][n_GAIN] = {0.};
    double errEntries[n_SiPM_tot][n_GAIN] = {0.};

    double Sigma_add[n_SiPM_tot][n_GAIN] = {0.};
    double errSigma_add[n_SiPM_tot][n_GAIN] = {0.};


    ///////////////////////////////////////////////////////////////////////////
    //      SiPM1
    ///////////////////////////////////////////////////////////////////////////

    n_SiPM = 0;


    // HV = 31 V [0]
    // Window for LED peak: (168, 180) ns  (minLED_amp, maxLED_amp)
    // Other: dleddt = 9; smooth_trace_bool = 0;
    //        histLED_low = -50; histLED_high = 700; histLED_binw = 0.7;
    //        fit_low = -10; fit_high = 28.7;
    H_peak_0[n_SiPM][0]          = 756.727;
    errH_peak_0[n_SiPM][0]       = 12.9998;
    Sigma_peak_0[n_SiPM][0]      = 2.23961;
    errSigma_peak_0[n_SiPM][0]   = 0.0332716;
    H_peak_1[n_SiPM][0]          = 1214.15;
    errH_peak_1[n_SiPM][0]       = 19.6876;
    Sigma_peak_1[n_SiPM][0]      = 2.81391;
    errSigma_peak_1[n_SiPM][0]   = 0.0445573;
    Mean_peak_1[n_SiPM][0]       = 6.73066;
    errMean_peak_1[n_SiPM][0]    = 0.0595634;
    Mean_peak_2[n_SiPM][0]       = 8.97026;
    errMean_peak_2[n_SiPM][0]    = 0.0682261;
    GAIN[n_SiPM][0]              = 2.23961;
    errGAIN[n_SiPM][0]           = 0.0332716;
    Sigma_add[n_SiPM][0]          = 1.7036;
    errSigma_add[n_SiPM][0]       = 0.0591894;
    Mean_hg[n_SiPM][2]           = 24.6891;
    errMean_hg[n_SiPM][2]        = 0.0694716;
    Std_hg[n_SiPM][2]            = 15.5342;
    errStd_hg[n_SiPM][2]         = 0.0491239;
    Entries[n_SiPM][0]           = 49999;


    // HV = 32 V [1]
    // Window for LED peak: (168, 180) ns  (minLED_amp, maxLED_amp)
    // Other: dleddt = 9; smooth_trace_bool = 0;
    //        histLED_low = -50; histLED_high = 700; histLED_binw = 0.7;
    //        fit_low = -10; fit_high = 34;
    H_peak_0[n_SiPM][1]          = 564.802;
    errH_peak_0[n_SiPM][1]       = 10.1474;
    Sigma_peak_0[n_SiPM][1]      = 2.56972;
    errSigma_peak_0[n_SiPM][1]   = 0.0326385;
    H_peak_1[n_SiPM][1]          = 979.947;
    errH_peak_1[n_SiPM][1]       = 12.4316;
    Sigma_peak_1[n_SiPM][1]      = 3.11386;
    errSigma_peak_1[n_SiPM][1]   = 0.0364095;
    Mean_peak_1[n_SiPM][1]       = 6.72997;
    errMean_peak_1[n_SiPM][1]    = 0.0554915;
    Mean_peak_2[n_SiPM][1]       = 9.29969;
    errMean_peak_2[n_SiPM][1]    = 0.0643784;
    GAIN[n_SiPM][1]              = 2.56972;
    errGAIN[n_SiPM][1]           = 0.0326385;
    Sigma_add[n_SiPM][1]         = 1.7586;
    errSigma_add[n_SiPM][1]      = 0.0433776;
    Mean_hg[n_SiPM][3]           = 32.991;
    errMean_hg[n_SiPM][3]        = 0.096609;
    Std_hg[n_SiPM][3]            = 21.6022;
    errStd_hg[n_SiPM][3]         = 0.0683129;
    Entries[n_SiPM][1]           = 49999;



    // HV = 33 V [2]
    // Window for LED peak: (168, 180) ns  (minLED_amp, maxLED_amp)
    // Other: dleddt = 9; smooth_trace_bool = 0;
    //        histLED_low = -50; histLED_high = 700; histLED_binw = 0.7;
    //        fit_low = -10; fit_high = 38;
    H_peak_0[n_SiPM][2]          = 443.125;
    errH_peak_0[n_SiPM][2]       = 8.93044;
    Sigma_peak_0[n_SiPM][2]      = 2.88207;
    errSigma_peak_0[n_SiPM][2]   = 0.0403952;
    H_peak_1[n_SiPM][2]          = 796.487;
    errH_peak_1[n_SiPM][2]       = 9.99911;
    Sigma_peak_1[n_SiPM][2]      = 3.46964;
    errSigma_peak_1[n_SiPM][2]   = 0.0434955;
    Mean_peak_1[n_SiPM][2]       = 6.61652;
    errMean_peak_1[n_SiPM][2]    = 0.0624507;
    Mean_peak_2[n_SiPM][2]       = 9.49859;
    errMean_peak_2[n_SiPM][2]    = 0.0743765;
    GAIN[n_SiPM][2]              = 2.88207;
    errGAIN[n_SiPM][2]           = 0.0403952;
    Sigma_add[n_SiPM][2]         = 1.93186;
    errSigma_add[n_SiPM][2]      = 0.0497065;
    Mean_hg[n_SiPM][4]           = 42.534;
    errMean_hg[n_SiPM][4]        = 0.129352;
    Std_hg[n_SiPM][4]            = 28.9238;
    errStd_hg[n_SiPM][4]         = 0.091466;
    Entries[n_SiPM][2]           = 49999;



    // HV = 34 V [3]
    // Window for LED peak: (168, 180) ns  (minLED_amp, maxLED_amp)
    // Other: dleddt = 9; smooth_trace_bool = 0;
    //        histLED_low = -50; histLED_high = 700; histLED_binw = 0.7;
    //        fit_low = -10; fit_high = 43;
    H_peak_0[n_SiPM][3]          = 351.779;
    errH_peak_0[n_SiPM][3]       = 7.1352;
    Sigma_peak_0[n_SiPM][3]      = 3.18861;
    errSigma_peak_0[n_SiPM][3]   = 0.0401266;
    H_peak_1[n_SiPM][3]          = 641.251;
    errH_peak_1[n_SiPM][3]       = 8.29728;
    Sigma_peak_1[n_SiPM][3]      = 3.89441;
    errSigma_peak_1[n_SiPM][3]   = 0.0440761;
    Mean_peak_1[n_SiPM][3]       = 6.16776;
    errMean_peak_1[n_SiPM][3]    = 0.0662144;
    Mean_peak_2[n_SiPM][3]       = 9.35637;
    errMean_peak_2[n_SiPM][3]    = 0.0774241;
    GAIN[n_SiPM][3]              = 3.18861;
    errGAIN[n_SiPM][3]           = 0.0401266;
    Sigma_add[n_SiPM][3]         = 2.23588;
    errSigma_add[n_SiPM][3]      = 0.0511766;
    Mean_hg[n_SiPM][5]           = 53.4209;
    errMean_hg[n_SiPM][5]        = 0.165406;
    Std_hg[n_SiPM][5]            = 36.9854;
    errStd_hg[n_SiPM][5]         = 0.116959;
    Entries[n_SiPM][3]           = 49999;


    // HV = 35.00 V [4]
    // Window for LED peak: (168, 180) ns  (minLED_amp, maxLED_amp)
    // Other: dleddt = 9; smooth_trace_bool = 0;
    //        histLED_low = -50; histLED_high = 700; histLED_binw = 0.7;
    //        fit_low = -10; fit_high = 49;
    H_peak_0[n_SiPM][4]          = 296.954;
    errH_peak_0[n_SiPM][4]       = 6.67128;
    Sigma_peak_0[n_SiPM][4]      = 3.40184;
    errSigma_peak_0[n_SiPM][4]   = 0.0500745;
    H_peak_1[n_SiPM][4]          = 530.895;
    errH_peak_1[n_SiPM][4]       = 7.18599;
    Sigma_peak_1[n_SiPM][4]      = 4.27458;
    errSigma_peak_1[n_SiPM][4]   = 0.0515737;
    Mean_peak_1[n_SiPM][4]       = 5.4591;
    errMean_peak_1[n_SiPM][4]    = 0.0750009;
    Mean_peak_2[n_SiPM][4]       = 8.86093;
    errMean_peak_2[n_SiPM][4]    = 0.0901808;
    GAIN[n_SiPM][4]              = 3.40184;
    errGAIN[n_SiPM][4]           = 0.0500745;
    Sigma_add[n_SiPM][4]         = 2.58834;
    errSigma_add[n_SiPM][4]      = 0.0540655;
    Mean_hg[n_SiPM][6]           = 65.0381;
    errMean_hg[n_SiPM][6]        = 0.204377;
    Std_hg[n_SiPM][6]            = 45.6996;
    errStd_hg[n_SiPM][6]         = 0.144516;
    Entries[n_SiPM][4]           = 49999;


    // HV = 36 V ([5])
    // Window for LED peak: (168, 180) ns  (minLED_amp, maxLED_amp)
    // Other: dleddt = 9; smooth_trace_bool = 0;
    //        histLED_low = -50; histLED_high = 700; histLED_binw = 0.7;
    //        fit_low = -10; fit_high = 60;
    H_peak_0[n_SiPM][5]          = 234.088;
    errH_peak_0[n_SiPM][5]       = 5.82906;
    Sigma_peak_0[n_SiPM][5]      = 3.79986;
    errSigma_peak_0[n_SiPM][5]   = 0.0672815;
    H_peak_1[n_SiPM][5]          = 393.178;
    errH_peak_1[n_SiPM][5]       = 5.96371;
    Sigma_peak_1[n_SiPM][5]      = 4.98827;
    errSigma_peak_1[n_SiPM][5]   = 0.0629806;
    Mean_peak_1[n_SiPM][5]       = 4.8069;
    errMean_peak_1[n_SiPM][5]    = 0.0952783;
    Mean_peak_2[n_SiPM][5]       = 8.60677;
    errMean_peak_2[n_SiPM][5]    = 0.116639;
    GAIN[n_SiPM][5]              = 3.79986;
    errGAIN[n_SiPM][5]           = 0.0672815;
    Sigma_add[n_SiPM][5]         = 3.2317;
    errSigma_add[n_SiPM][5]      = 0.0564978;
    Mean_hg[n_SiPM][7]           = 78.4931;
    errMean_hg[n_SiPM][7]        = 0.249375;
    Std_hg[n_SiPM][7]            = 55.7614;
    errStd_hg[n_SiPM][7]         = 0.176335;
    Entries[n_SiPM][5]           = 49999;


    ///////////////////////////////////////////////////////////////////////////
    //      SiPM2
    ///////////////////////////////////////////////////////////////////////////

    n_SiPM = 1;

    // HV = 31 [0]
    // Window for LED peak: (174, 184) ns  (minLED_amp, maxLED_amp)
    // Other: dleddt = 9; smooth_trace_bool = 0;
    //        histLED_low = -50; histLED_high = 700; histLED_binw = 0.7;
    //        fit_low = -20; fit_high = 27;
    H_peak_0[n_SiPM][0]          = 710.271;
    errH_peak_0[n_SiPM][0]       = 12.8005;
    Sigma_peak_0[n_SiPM][0]      = 2.00111;
    errSigma_peak_0[n_SiPM][0]   = 0.0274682;
    H_peak_1[n_SiPM][0]          = 1118.88;
    errH_peak_1[n_SiPM][0]       = 17.4156;
    Sigma_peak_1[n_SiPM][0]      = 2.63774;
    errSigma_peak_1[n_SiPM][0]   = 0.0339583;
    Mean_peak_1[n_SiPM][0]       = 4.49423;
    errMean_peak_1[n_SiPM][0]    = 0.0481392;
    Mean_peak_2[n_SiPM][0]       = 6.49534;
    errMean_peak_2[n_SiPM][0]    = 0.0554246;
    GAIN[n_SiPM][0]              = 2.00111;
    errGAIN[n_SiPM][0]           = 0.0274682;
    Sigma_add[n_SiPM][0]         = 1.7185;
    errSigma_add[n_SiPM][0]      = 0.0411549;
    Mean_hg[n_SiPM][2]           = 24.3421;
    errMean_hg[n_SiPM][2]        = 0.0710515;
    Std_hg[n_SiPM][2]            = 15.8874;
    errStd_hg[n_SiPM][2]         = 0.050241;
    Entries[n_SiPM][0]           = 49999;








    // HV = 32 [1]
    // Window for LED peak: (174, 184) ns  (minLED_amp, maxLED_amp)
    // Other: dleddt = 9; smooth_trace_bool = 0;
    //        histLED_low = -50; histLED_high = 700; histLED_binw = 0.7;
    //        fit_low = -20; fit_high = 32;
    H_peak_0[n_SiPM][1]          = 526.248;
    errH_peak_0[n_SiPM][1]       = 10.793;
    Sigma_peak_0[n_SiPM][1]      = 2.16341;
    errSigma_peak_0[n_SiPM][1]   = 0.0302084;
    H_peak_1[n_SiPM][1]          = 875.738;
    errH_peak_1[n_SiPM][1]       = 12.4199;
    Sigma_peak_1[n_SiPM][1]      = 2.88812;
    errSigma_peak_1[n_SiPM][1]   = 0.0326515;
    Mean_peak_1[n_SiPM][1]       = 4.31658;
    errMean_peak_1[n_SiPM][1]    = 0.0498871;
    Mean_peak_2[n_SiPM][1]       = 6.47998;
    errMean_peak_2[n_SiPM][1]    = 0.0583204;
    GAIN[n_SiPM][1]              = 2.16341;
    errGAIN[n_SiPM][1]           = 0.0302084;
    Sigma_add[n_SiPM][1]         = 1.91335;
    errSigma_add[n_SiPM][1]      = 0.035531;
    Mean_hg[n_SiPM][3]           = 33.6592;
    errMean_hg[n_SiPM][3]        = 0.0998505;
    Std_hg[n_SiPM][3]            = 22.327;
    errStd_hg[n_SiPM][3]         = 0.070605;
    Entries[n_SiPM][1]           = 49999;







    // HV = 33 V [2]
    // Window for LED peak: (174, 184) ns  (minLED_amp, maxLED_amp)
    // Other: dleddt = 9; smooth_trace_bool = 0;
    //        histLED_low = -50; histLED_high = 700; histLED_binw = 0.7;
    //        fit_low = -20; fit_high = 35;
    H_peak_0[n_SiPM][2]          = 423.027;
    errH_peak_0[n_SiPM][2]       = 9.31801;
    Sigma_peak_0[n_SiPM][2]      = 2.2906;
    errSigma_peak_0[n_SiPM][2]   = 0.0326721;
    H_peak_1[n_SiPM][2]          = 704.377;
    errH_peak_1[n_SiPM][2]       = 9.93279;
    Sigma_peak_1[n_SiPM][2]      = 3.12335;
    errSigma_peak_1[n_SiPM][2]   = 0.035086;
    Mean_peak_1[n_SiPM][2]       = 3.92169;
    errMean_peak_1[n_SiPM][2]    = 0.0533051;
    Mean_peak_2[n_SiPM][2]       = 6.21229;
    errMean_peak_2[n_SiPM][2]    = 0.0625212;
    GAIN[n_SiPM][2]              = 2.2906;
    errGAIN[n_SiPM][2]           = 0.0326721;
    Sigma_add[n_SiPM][2]         = 2.12332;
    errSigma_add[n_SiPM][2]      = 0.0377012;
    Mean_hg[n_SiPM][4]           = 44.2844;
    errMean_hg[n_SiPM][4]        = 0.133826;
    Std_hg[n_SiPM][4]            = 29.924;
    errStd_hg[n_SiPM][4]         = 0.094629;
    Entries[n_SiPM][2]           = 49999;





    // HV = 34 V [3]
    // Window for LED peak: (174, 184) ns  (minLED_amp, maxLED_amp)
    // Other: dleddt = 9; smooth_trace_bool = 0;
    //        histLED_low = -50; histLED_high = 700; histLED_binw = 0.7;
    //        fit_low = -20; fit_high = 40;
    H_peak_0[n_SiPM][3]          = 318.584;
    errH_peak_0[n_SiPM][3]       = 7.83841;
    Sigma_peak_0[n_SiPM][3]      = 2.63653;
    errSigma_peak_0[n_SiPM][3]   = 0.0459871;
    H_peak_1[n_SiPM][3]          = 537.034;
    errH_peak_1[n_SiPM][3]       = 8.0054;
    Sigma_peak_1[n_SiPM][3]      = 3.5292;
    errSigma_peak_1[n_SiPM][3]   = 0.04528;
    Mean_peak_1[n_SiPM][3]       = 3.67286;
    errMean_peak_1[n_SiPM][3]    = 0.0658013;
    Mean_peak_2[n_SiPM][3]       = 6.30939;
    errMean_peak_2[n_SiPM][3]    = 0.0802785;
    GAIN[n_SiPM][3]              = 2.63653;
    errGAIN[n_SiPM][3]           = 0.0459871;
    Sigma_add[n_SiPM][3]         = 2.34605;
    errSigma_add[n_SiPM][3]      = 0.0443708;
    Mean_hg[n_SiPM][5]           = 56.4931;
    errMean_hg[n_SiPM][5]        = 0.170828;
    Std_hg[n_SiPM][5]            = 38.1978;
    errStd_hg[n_SiPM][5]         = 0.120793;
    Entries[n_SiPM][3]           = 49999;





    // HV = 35.00 V [4]
    // Window for LED peak: (174, 184) ns  (minLED_amp, maxLED_amp)
    // Other: dleddt = 9; smooth_trace_bool = 0;
    //        histLED_low = -50; histLED_high = 700; histLED_binw = 0.7;
    //        fit_low = -20; fit_high = 44.5;
    H_peak_0[n_SiPM][4]          = 126.663;
    errH_peak_0[n_SiPM][4]       = 3.86506;
    Sigma_peak_0[n_SiPM][4]      = 4.81968;
    errSigma_peak_0[n_SiPM][4]   = 0.107605;
    H_peak_1[n_SiPM][4]          = 297.653;
    errH_peak_1[n_SiPM][4]       = 5.75533;
    Sigma_peak_1[n_SiPM][4]      = 5.67565;
    errSigma_peak_1[n_SiPM][4]   = 0.111621;
    Mean_peak_1[n_SiPM][4]       = 1.22104;
    errMean_peak_1[n_SiPM][4]    = 0.180563;
    Mean_peak_2[n_SiPM][4]       = 6.04072;
    errMean_peak_2[n_SiPM][4]    = 0.210195;
    GAIN[n_SiPM][4]              = 4.81968;
    errGAIN[n_SiPM][4]           = 0.107605;
    Sigma_add[n_SiPM][4]         = 2.99727;
    errSigma_add[n_SiPM][4]      = 0.121389;
    Mean_hg[n_SiPM][6]           = 66.7518;
    errMean_hg[n_SiPM][6]        = 0.209945;
    Std_hg[n_SiPM][6]            = 46.9446;
    errStd_hg[n_SiPM][6]         = 0.148453;
    Entries[n_SiPM][4]           = 49999;





    // HV = 36 V ([5])
    // Window for LED peak: (174, 184) ns  (minLED_amp, maxLED_amp)
    // Other: dleddt = 9; smooth_trace_bool = 0;
    //        histLED_low = -50; histLED_high = 700; histLED_binw = 0.7;
    //        fit_low = -20; fit_high = 44.5;
    H_peak_0[n_SiPM][5]          = 202.892;
    errH_peak_0[n_SiPM][5]       = 6.00704;
    Sigma_peak_0[n_SiPM][5]      = 3.1787;
    errSigma_peak_0[n_SiPM][5]   = 0.0678594;
    H_peak_1[n_SiPM][5]          = 352.986;
    errH_peak_1[n_SiPM][5]       = 6.05908;
    Sigma_peak_1[n_SiPM][5]      = 4.34391;
    errSigma_peak_1[n_SiPM][5]   = 0.0677024;
    Mean_peak_1[n_SiPM][5]       = 2.80963;
    errMean_peak_1[n_SiPM][5]    = 0.0974139;
    Mean_peak_2[n_SiPM][5]       = 5.98833;
    errMean_peak_2[n_SiPM][5]    = 0.11872;
    GAIN[n_SiPM][5]              = 3.1787;
    errGAIN[n_SiPM][5]           = 0.0678594;
    Sigma_add[n_SiPM][5]         = 2.96064;
    errSigma_add[n_SiPM][5]      = 0.067521;
    Mean_hg[n_SiPM][7]           = 84.0498;
    errMean_hg[n_SiPM][7]        = 0.257809;
    Std_hg[n_SiPM][7]            = 57.6472;
    errStd_hg[n_SiPM][7]         = 0.182298;
    Entries[n_SiPM][5]           = 49999;







    //------------------------------

    // select only points for n_GAIN
    // for the mean (NOT used at the time being) I use more HVs
    for (int i = 0; i < n_GAIN; i++) {
        HV_LASER_GAIN[i] = HV_LASER[i+n_MEAN-n_GAIN];
        errHV_LASER_GAIN[i] = errHV_LASER[i+n_MEAN-n_GAIN];
    }

    ///////////////////////////////////////////////////////////////////////////
    //      CROSS TALK PROBABILITY AND RATIO GAIN
    ///////////////////////////////////////////////////////////////////////////


    double R_gain[n_SiPM_tot][n_GAIN] = {0.};
    double errR_gain[n_SiPM_tot][n_GAIN] = {0.};

    for(int i=0; i<n_GAIN;i++){
        for(n_SiPM=0; n_SiPM<n_SiPM_tot_LASER; n_SiPM++){
            //////////////////////
            ///   CROSS TALK   ///
            //////////////////////
            errEntries[n_SiPM][i]=TMath::Power(Entries[n_SiPM][i],0.5);
            Area0[n_SiPM][i]=H_peak_0[n_SiPM][i]*Sigma_peak_0[n_SiPM][i]*TMath::Power(2*TMath::Pi(),0.5)/1;
            Prob_0pe[n_SiPM][i]=Area0[n_SiPM][i]/Entries[n_SiPM][i];
            Mu[n_SiPM][i]=-TMath::Log(Prob_0pe[n_SiPM][i]);
            Prob_1pe[n_SiPM][i]=Mu[n_SiPM][i]*TMath::Exp(-Mu[n_SiPM][i]);
            Area1[i]=H_peak_1[n_SiPM][i]*Sigma_peak_1[n_SiPM][i]*TMath::Power(2*TMath::Pi(),0.5)/1;
            Prob_1peS[n_SiPM][i]=Area1[i]/Entries[n_SiPM][i];
            Prob_Cross_Talk[n_SiPM][i]=1-(Prob_1peS[n_SiPM][i]/Prob_1pe[n_SiPM][i]);

            errArea1[n_SiPM][i]=Area1[i]*TMath::Power(TMath::Power(errH_peak_1[n_SiPM][i]/H_peak_1[n_SiPM][i],2)+TMath::Power(errSigma_peak_1[n_SiPM][i]/Sigma_peak_1[n_SiPM][i],2),0.5);
            errProb_1peS[n_SiPM][i]=errArea1[n_SiPM][i]/Entries[n_SiPM][i];
            errMu[n_SiPM][i]=errArea0[n_SiPM][i]/(Prob_0pe[n_SiPM][i]*Entries[n_SiPM][i]);
            errProb_1pe[n_SiPM][i]=TMath::Exp(-Mu[n_SiPM][i])*TMath::Abs(1-Mu[n_SiPM][i])*errMu[n_SiPM][i];
            errArea0[n_SiPM][i]=Area0[n_SiPM][i]*TMath::Power(TMath::Power(errH_peak_0[n_SiPM][i]/H_peak_0[n_SiPM][i],2)+TMath::Power(errSigma_peak_0[n_SiPM][i]/Sigma_peak_0[n_SiPM][i],2),0.5);
            errProb_0pe[n_SiPM][i]=errArea0[n_SiPM][i]/Entries[n_SiPM][i];
            errProb_Cross_Talk[n_SiPM][i]=TMath::Power(TMath::Power(errProb_1peS[n_SiPM][i]/Prob_1pe[n_SiPM][i],2)+TMath::Power(Prob_1peS[n_SiPM][i]*errProb_1pe[n_SiPM][i]/Prob_1pe[n_SiPM][i],2),0.5);

            //////////////////////
            ///   RATIO GAIN   ///
            //////////////////////
            R_gain[n_SiPM][i]    = GAIN[n_SiPM][i] / Sigma_add[n_SiPM][i];
            errR_gain[n_SiPM][i] = (errGAIN[n_SiPM][i]*errGAIN[n_SiPM][i])/(Sigma_add[n_SiPM][i]*Sigma_add[n_SiPM][i]);
            errR_gain[n_SiPM][i]+= (errSigma_add[n_SiPM][i]*errSigma_add[n_SiPM][i])*(GAIN[n_SiPM][i]*GAIN[n_SiPM][i])/TMath::Power(Sigma_add[n_SiPM][i],4);
            errR_gain[n_SiPM][i] = TMath::Sqrt(errR_gain[n_SiPM][i]);
        }
    }


    ///////////////////////////////////////////////////////////////////////////
    //          FIXED ERROR
    ///////////////////////////////////////////////////////////////////////////

    if(fix_error_bool){
        // double err_fix_CT = 0.015;
        double err_fix_CT = 0.02;
        for(int i=0; i<n_GAIN; i++){
            for(n_SiPM=0; n_SiPM<n_SiPM_tot_LASER; n_SiPM++){
                errProb_Cross_Talk[n_SiPM][i] = err_fix_CT;
            }
        }
    } // end fix_error_bool

    ///////////////////////////////////////////////////////////////////////////


    TGraphErrors *gV_Rg[n_SiPM_tot];
    TGraphErrors *gV_PCT[n_SiPM_tot];
    TCanvas *cV_Rg[n_SiPM_tot];
    TCanvas *cV_PCT[n_SiPM_tot];
    string title_cV_Rg;
    string title_cV_PCT;

    for(n_SiPM=0; n_SiPM<n_SiPM_tot_LASER; n_SiPM++){
        //------------------------------
        // R_Gain
        //------------------------------
        gV_Rg[n_SiPM]  = new TGraphErrors(n_GAIN, HV_LASER_GAIN, R_gain[n_SiPM], errHV_LASER_GAIN, errR_gain[n_SiPM]);
        gV_Rg[n_SiPM]->SetMarkerStyle(21);
        gV_Rg[n_SiPM]->SetMarkerColor(kOrange+2);
        gV_Rg[n_SiPM]->SetTitle();
        gV_Rg[n_SiPM]->GetXaxis()->SetTitle("Bias Voltage (V)");
        gV_Rg[n_SiPM]->GetYaxis()->SetTitle("R_{gain}");

        title_cV_Rg="cV_Rg_"+to_string(n_SiPM);
        cV_Rg[n_SiPM] = new TCanvas(title_cV_Rg.c_str(), title_cV_Rg.c_str(),w,h);
        cV_Rg[n_SiPM]->SetGrid();
        gV_Rg[n_SiPM]->Draw("AP");


        //------------------------------
        // Prob_Cross_Talk
        //------------------------------
        gV_PCT[n_SiPM]  = new TGraphErrors(n_GAIN, HV_LASER_GAIN, Prob_Cross_Talk[n_SiPM], errHV_LASER_GAIN, errProb_Cross_Talk[n_SiPM]);
        gV_PCT[n_SiPM]->SetMarkerStyle(21);
        // gV_PCT[n_SiPM]->SetMarkerSize(2);
        gV_PCT[n_SiPM]->SetMarkerColor(kOrange+2);
        gV_PCT[n_SiPM]->SetTitle();
        gV_PCT[n_SiPM]->GetXaxis()->SetTitle("Bias Voltage (V)");
        gV_PCT[n_SiPM]->GetYaxis()->SetTitle("Cross Talk");

        title_cV_PCT="cV_PCT_"+to_string(n_SiPM);
        cV_PCT[n_SiPM] = new TCanvas(title_cV_PCT.c_str(), title_cV_PCT.c_str(),w,h);
        cV_PCT[n_SiPM]->SetGrid();
        gV_PCT[n_SiPM]->Draw("AP");
    }


    cout<<endl<<endl;
    cout<<"//----------------"<<endl;
    cout<<"// FBK HD3-2 LASER"<<endl;
    cout<<"//----------------"<<endl;
    cout<<"double HV_FBK_LASER[] = {"; for(int i=0; i<n_GAIN-1; i++) cout<<HV_LASER_GAIN[i]<<", "; cout<<HV_LASER_GAIN[n_GAIN-1]<<"};"<<endl;
    cout<<"double errHV_FBK_LASER[] = {"; for(int i=0; i<n_GAIN-1; i++) cout<<errHV_LASER_GAIN[i]<<", "; cout<<errHV_LASER_GAIN[n_GAIN-1]<<"};"<<endl;
    cout<<"double CT_1_FBK_LASER[] = {"; for(int i=0; i<n_GAIN-1; i++) cout<<Prob_Cross_Talk[0][i]<<", "; cout<<Prob_Cross_Talk[0][n_GAIN-1]<<"};"<<endl;
    cout<<"double errCT_1_FBK_LASER[] = {"; for(int i=0; i<n_GAIN-1; i++) cout<<errProb_Cross_Talk[0][i]<<", "; cout<<errProb_Cross_Talk[0][n_GAIN-1]<<"};"<<endl;



  // /////////////////////////////////////////////////////////////////////////////
  // //    MERGE
  // /////////////////////////////////////////////////////////////////////////////
  //
  //     TCanvas *cV_CT_DARK_LASER = new TCanvas("cV_CT_DARK_LASER", "cV_CT_DARK_LASER",w,h);
  //     cV_CT_DARK_LASER->SetGrid();
  //     TMultiGraph *mgCT_DARK_LASER = new TMultiGraph("mgCT_DARK_LASER", ";Bias Voltage (V);Cross Talk");
  //     mgCT_DARK_LASER->Add(gV_CT_1);
  //     // mgCT_DARK_LASER->Add(gV_CT_2);
  //     // mgCT_DARK_LASER->Add(gV_CT_3);
  //     mgCT_DARK_LASER->Add(gV_CT_Del_1);
  //     // mgCT_DARK_LASER->Add(gV_CT_Del_2);
  //     // mgCT_DARK_LASER->Add(gV_CT_Del_3);
  //     mgCT_DARK_LASER->Add(gV_PCT[0]);
  //     mgCT_DARK_LASER->Draw("AP");
  //
  //     auto legendCT_DARK_LASER = new TLegend(0.15,0.70,0.35,0.85);
  //     legendCT_DARK_LASER->AddEntry(gV_CT_1,    "HD3-2 (1), from DCR","p");
  //     legendCT_DARK_LASER->AddEntry(gV_CT_Del_1,"HD3-2 (1), from DELAYS","p");
  //     legendCT_DARK_LASER->AddEntry(gV_PCT[0],     "HD3-2 (1), from LASER","p");
  //     legendCT_DARK_LASER->Draw();
  //

      //------------------------------
      // For LaTeX
      //------------------------------
      n_SiPM=0;
      cout<<endl<<endl;
      cout<<"For LaTeX"<<endl<<endl;
      cout<<"\% From file DCR_CrossTalk_FBK_HD3_2_DARK_LASER_data_2018_07.C"<<endl;
      cout<<"\% From LASER"<<endl;
      cout<<"\% HV & CT & R_{gain}   \\\\"<<endl;
      for(int i=0; i<n_GAIN; i++){
          printf("$ %.2f \\pm %.2f $ & $ %.2f \\pm %.2f $ & $ %.2f \\pm %.2f $ \\\\ \n", HV_LASER_GAIN[i], errHV_LASER_GAIN[i], Prob_Cross_Talk[n_SiPM][i], errProb_Cross_Talk[n_SiPM][i], R_gain[n_SiPM][i], errR_gain[n_SiPM][i]);
      }


//*****************************************************************************
//*****************************************************************************
//*****************************************************************************
//*****************************************************************************
//*****************************************************************************
//*****************************************************************************
//*****************************************************************************
//*****************************************************************************



/////////////////////////////////////////////////////////////////////////////
//    TEST Z
/////////////////////////////////////////////////////////////////////////////

    ///////////////////
    // DARK
    ///////////////////

    // cout<<endl<<endl;
    // cout<<"/////////////////////////////////////////////////"<<endl;
    // cout<<"//     TEST Z DARK"<<endl;
    // cout<<"/////////////////////////////////////////////////"<<endl;
    //
    // // TEST Z
    // int N_Z = n_DCR_1;
    // double testZ_DCR[3][N_Z];
    // double testZ_CT[3][N_Z];
    // n_SiPM = 0;
    //
    //
    // n_SiPM = 0;
    // for(int i=0; i<N_Z; i++){
    //   testZ_DCR[n_SiPM][i] = TMath::Abs((DCR_1[i] - DCR_Del_1[i]) / (TMath::Sqrt( errDCR_1[i]*errDCR_1[i] + errDCR_Del_1[i]*errDCR_Del_1[i] )));
    // }
    //
    // n_SiPM = 1;
    // for(int i=0; i<N_Z; i++){
    //   testZ_DCR[n_SiPM][i] = TMath::Abs((DCR_2[i] - DCR_Del_2[i]) / (TMath::Sqrt( errDCR_2[i]*errDCR_2[i] + errDCR_Del_2[i]*errDCR_Del_2[i] )));
    // }
    //
    // n_SiPM = 2;
    // for(int i=0; i<N_Z; i++){
    //   testZ_DCR[n_SiPM][i] = TMath::Abs((DCR_2[i] - DCR_Del_2[i]) / (TMath::Sqrt( errDCR_2[i]*errDCR_2[i] + errDCR_Del_2[i]*errDCR_Del_2[i] )));
    // }
    //
    // n_SiPM = 0;
    // for(int i=0; i<N_Z; i++){
    //   testZ_CT[n_SiPM][i] = TMath::Abs((CT_1[i] - CT_Del_1[i]) / (TMath::Sqrt( errCT_1[i]*errCT_1[i] + errCT_Del_1[i]*errCT_Del_1[i] )));
    // }
    //
    // n_SiPM = 1;
    // for(int i=0; i<N_Z; i++){
    //   testZ_CT[n_SiPM][i] = TMath::Abs((CT_2[i] - CT_Del_2[i]) / (TMath::Sqrt( errCT_2[i]*errCT_2[i] + errCT_Del_2[i]*errCT_Del_2[i] )));
    // }
    //
    // n_SiPM = 2;
    // for(int i=0; i<N_Z; i++){
    //   testZ_CT[n_SiPM][i] = TMath::Abs((CT_2[i] - CT_Del_2[i]) / (TMath::Sqrt( errCT_2[i]*errCT_2[i] + errCT_Del_2[i]*errCT_Del_2[i] )));
    // }
    //
    //
    // for(int n=0; n<n_SiPM_tot; n++){
    //   cout<<"SiPM "<<n+1<<endl;
    //   for(int i=0; i<N_Z; i++){
    //       cout<<testZ_DCR[n][i];
    //       if(testZ_DCR[n][i] > 1.96) cout<<"\tTEST NOT PASSED"<<endl;
    //       else cout<<endl;
    //
    //       cout<<testZ_CT[n][i];
    //       if(testZ_CT[n][i] > 1.96) cout<<"\tTEST NOT PASSED "<<endl;
    //       else cout<<endl;
    //   }
    //   cout<<endl;
    // }
    //
    //
    // /////////////////////////
    // // DARK vs LASER
    // /////////////////////////
    // cout<<endl<<endl;
    // cout<<"/////////////////////////////////////////////////"<<endl;
    // cout<<"//     TEST Z DARK and LASER"<<endl;
    // cout<<"/////////////////////////////////////////////////"<<endl;
    // N_Z = n_GAIN;
    // double testZ_CT_DarkLaser = 0.;
    //
    // n_SiPM = 0;
    //
    // cout<<endl;
    // cout<<"CT DARK CNT vs CT DARK DEL SiPM "<<n_SiPM+1<<endl;
    // for(int i=0; i<N_Z; i++){
    //   testZ_CT_DarkLaser = TMath::Abs((CT_1[i] - CT_Del_1[i]) / (TMath::Sqrt( errCT_1[i]*errCT_1[i] + errCT_Del_1[i]*errCT_Del_1[i] )));
    //
    //   cout<<testZ_CT_DarkLaser;
    //   if(testZ_CT_DarkLaser > 1.96) cout<<"\tTEST NOT PASSED"<<endl;
    //   else cout<<endl;
    // }
    //
    // cout<<endl;
    // cout<<"CT DARK CNT vs CT LASER SiPM "<<n_SiPM+1<<endl;
    // for(int i=0; i<N_Z; i++){
    //   testZ_CT_DarkLaser = TMath::Abs((CT_1[i] - Prob_Cross_Talk[n_SiPM][i]) / (TMath::Sqrt( errCT_1[i]*errCT_1[i] + errProb_Cross_Talk[n_SiPM][i]*errProb_Cross_Talk[n_SiPM][i] )));
    //
    //   cout<<testZ_CT_DarkLaser;
    //   if(testZ_CT_DarkLaser > 1.96) cout<<"\tTEST NOT PASSED"<<endl;
    //   else cout<<endl;
    // }
    //
    // cout<<endl;
    // cout<<"CT DARK DEL vs CT LASER SiPM "<<n_SiPM+1<<endl;
    // for(int i=0; i<N_Z; i++){
    //   testZ_CT_DarkLaser = TMath::Abs((CT_Del_1[i] - Prob_Cross_Talk[n_SiPM][i]) / (TMath::Sqrt( errCT_Del_1[i]*errCT_Del_1[i] + errProb_Cross_Talk[n_SiPM][i]*errProb_Cross_Talk[n_SiPM][i] )));
    //
    //   cout<<testZ_CT_DarkLaser;
    //   if(testZ_CT_DarkLaser > 1.96) cout<<"\tTEST NOT PASSED"<<endl;
    //   else cout<<endl;
    // }




    /////////////////////////
    // OLD STUFF:
    ////////////////////////

    // for (int i=0; i<n_GAIN; i++){
    //
    //     // GainWeighted[n_SiPM]
    //     // GainWeighted[n_SiPM][i] = GAIN[n_SiPM][i]/TMath::Sqrt(Sigma_peak_1[n_SiPM][i]*Sigma_peak_1[n_SiPM][i]-Sigma_peak_0[n_SiPM][i]*Sigma_peak_0[n_SiPM][i]);
    //     // errGainWeighted[n_SiPM][i] = errGAIN[n_SiPM][i]*errGAIN[n_SiPM][i]/(Sigma_peak_1[n_SiPM][i]*Sigma_peak_1[n_SiPM][i]-Sigma_peak_0[n_SiPM][i]*Sigma_peak_0[n_SiPM][i]);
    //     // errGainWeighted[n_SiPM][i] += errSigma_peak_1[n_SiPM][i]*errSigma_peak_1[n_SiPM][i]*Sigma_peak_1[n_SiPM][i]*Sigma_peak_1[n_SiPM][i]*GAIN[n_SiPM][i]*GAIN[n_SiPM][i]/TMath::Power(Sigma_peak_1[n_SiPM][i]*Sigma_peak_1[n_SiPM][i]-Sigma_peak_0[n_SiPM][i]*Sigma_peak_0[n_SiPM][i],3);
    //     // errGainWeighted[n_SiPM][i] += errSigma_peak_0[n_SiPM][i]*errSigma_peak_0[n_SiPM][i]*Sigma_peak_0[n_SiPM][i]*Sigma_peak_0[n_SiPM][i]*GAIN[n_SiPM][i]*GAIN[n_SiPM][i]/TMath::Power(Sigma_peak_1[n_SiPM][i]*Sigma_peak_1[n_SiPM][i]-Sigma_peak_0[n_SiPM][i]*Sigma_peak_0[n_SiPM][i],3);
    //     // errGainWeighted[n_SiPM][i] = TMath::Sqrt(errGainWeighted[n_SiPM][i]);
    //
    //     // Integral_GAIN
    //     // Integral_GAIN[n_SiPM][i]=TMath::Abs(Integral[n_SiPM][i])/(GAIN[n_SiPM][i]*Cross_Talk[i]);
    //     // errIntegral_GAIN[n_SiPM][i]=Integral[n_SiPM][i]/(GAIN[n_SiPM][i]*Cross_Talk[i])*TMath::Power(errGAIN[n_SiPM][i]*errGAIN[n_SiPM][i]/(GAIN[n_SiPM][i]*GAIN[n_SiPM][i])+errCross_Talk[i]*errCross_Talk[i]/(Cross_Talk[i]*Cross_Talk[i]),0.5);
    //     //Integral_GAIN[n_SiPM][i]=TMath::Abs(Integral[n_SiPM][i])/(GAIN[n_SiPM][i]);
    //     //errIntegral_GAIN[n_SiPM][i]=Integral[n_SiPM][i]*errGAIN[n_SiPM][i]/(GAIN[n_SiPM][i]*GAIN[n_SiPM][i]);
    //
    //     //Probabilità CrossTalk LED
    //
    // }
    // double temp;
    // for (int i = 0; i < n_MEAN; i++) {
    //   // Mean_hg[n_SiPM]/St_dev
    //   Mean_hist_St_dev[n_SiPM][i]= Mean_hg[n_SiPM][i]/Std_hg[n_SiPM][i];
    //   errMean_hist_St_dev[n_SiPM][i] = TMath::Power(Std_hg[n_SiPM][i],-1)*TMath::Power(errMean_hg[n_SiPM][i]*errMean_hg[n_SiPM][i]+Mean_hg[n_SiPM][i]*Mean_hg[n_SiPM][i]*errStd_hg[n_SiPM][i]*errStd_hg[n_SiPM][i]/(Std_hg[n_SiPM][i]*Std_hg[n_SiPM][i]),0.5);
    //
    //   // temp = errMean_hg[n_SiPM][i]*errMean_hg[n_SiPM][i]+Mean_hg[n_SiPM][i]*Mean_hg[n_SiPM][i]*errStd_hg[n_SiPM][i]*errStd_hg[n_SiPM][i]/(Std_hg[n_SiPM][i]*Std_hg[n_SiPM][i]);
    //   // temp = TMath::Power(temp, 0.5);
    //   // temp *= 1/Std_hg[n_SiPM][i];
    //   //
    // }

    //------------------------------
    // Gain_Weighted
    //------------------------------

    // TGraphErrors *gV_GW  = new TGraphErrors(n_GAIN, HV_LASER_GAIN, GainWeighted[0], errHV_LASER_GAIN, errGainWeighted[0]);
    //
    //
    // //------------------------------
    //
    // gV_GW->SetMarkerStyle(21);
    // gV_GW->SetMarkerColor(kOrange+2);
    // gV_GW->SetTitle();
    // gV_GW->GetXaxis()->SetTitle("Bias Voltage (V)");
    // gV_GW->GetYaxis()->SetTitle("Weighted Gain");
    //
    // //------------------------------
    //
    // TCanvas *cV_GW = new TCanvas("cV_GW", "cV_GW",w,h);
    // cV_GW->SetGrid();
    // gV_GW->Draw("AP");



    // check
    // TCanvas *cV_check_GW_Rg = new TCanvas("cV_check_GW_Rg", "cV_check_GW_Rg",w,h);
    // TMultiGraph *mgcheck_GW_Rg = new TMultiGraph("mgcheck_GW_Rg", ";Bias Voltage (V);Cross Talk");
    // gV_Rg->SetLineColor(kOrange+2);
    // mgcheck_GW_Rg->Add(gV_GW);
    // mgcheck_GW_Rg->Add(gV_Rg);
    // mgcheck_GW_Rg->Draw("AP");

    //------------------------------
    // MS
    //------------------------------
    // TGraphErrors *gV_MS  = new TGraphErrors(n_MEAN, HV_LASER, Mean_hist_St_dev[n_SiPM], errHV_LASER, errMean_hist_St_dev[n_SiPM]);
    //
    // gV_MS->SetMarkerStyle(21);
    // gV_MS->SetMarkerColor(kOrange+2);
    // gV_MS->SetTitle();
    // gV_MS->GetXaxis()->SetTitle("Bias Voltage (V)");
    // gV_MS->GetYaxis()->SetTitle("Mean/St_Dev ()");


    // TCanvas *cV_MS = new TCanvas("cV_MS", "cV_MS",w,h);
    // cV_MS->SetGrid();
    // gV_MS->Draw("AP");


    //*****************************************************************************
    //*****************************************************************************
    //*****************************************************************************
    //*****************************************************************************
    //*****************************************************************************
    //*****************************************************************************
    //*****************************************************************************
    //*****************************************************************************



}



int find_index(double v[], int N, double value){
    int index =0;


    double epsilon = 0.1;

    for(int i=0; i<N; i++){
        if((v[i]>value-epsilon) && (v[i]<value+epsilon)){
            index = i;
            // break;
        }
    }

    return index;

}
