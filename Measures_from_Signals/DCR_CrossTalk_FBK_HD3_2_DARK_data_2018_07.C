/******************************************************************************\
 * DCR_CrossTalk_FBK_HD3_2_from_cnt_data_2018_07.C
 *
 * Values obtained by Ana_Traces_SiPM.cxx (version of 07/08/2018, 1)
 *
 * KEY POINTS:
 *  > DCR_CT_1SiPM_nHVs(...)
 *  > dleddt = 6
 *  > NO trace smoothing
 *  > thr at 0.5pe and 1.5 pe set manually
 *  > min_thr_to_find_peaks = 8;  //first thr value in the DCR vs thr plot (mV)
 *  > max_thr_to_find_peaks = 80; //last thr value in the DCR vs thr plot (mV)
 *  > Area = 36; // 6x6 mm^2
 *
 *
 *
 * FILES ANALIZED:
 *
 *      20180725_HD3-2_01_DARK_AgilentE3641A_29.00_AS_2_100000ev_01.dat
 *      20180725_HD3-2_01_DARK_AgilentE3641A_30.00_AS_2_100000ev_01.dat
 *      20180725_HD3-2_01_DARK_AgilentE3641A_31.00_AS_2_100000ev_01.dat
 *      20180725_HD3-2_01_DARK_AgilentE3641A_32.00_AS_2_100000ev_01.dat
 *      20180725_HD3-2_01_DARK_AgilentE3641A_33.00_AS_2_100000ev_01.dat
 *      20180725_HD3-2_01_DARK_AgilentE3641A_34.00_AS_2_100000ev_01.dat
 *      20180725_HD3-2_01_DARK_AgilentE3641A_35.00_AS_2_100000ev_01.dat
 *      20180725_HD3-2_01_DARK_AgilentE3641A_36.00_AS_2_100000ev_01.dat
 *      20180725_HD3-2_01_DARK_AgilentE3641A_37.00_AS_2_100000ev_01.dat
 *      20180725_HD3-2_02_DARK_AgilentE3641A_29.00_AS_2_100000ev_01.dat
 *
 *      20180725_HD3-2_02_DARK_AgilentE3641A_30.00_AS_2_100000ev_01.dat
 *      20180725_HD3-2_02_DARK_AgilentE3641A_31.00_AS_2_100000ev_01.dat
 *      20180725_HD3-2_02_DARK_AgilentE3641A_32.00_AS_2_100000ev_01.dat
 *      20180725_HD3-2_02_DARK_AgilentE3641A_33.00_AS_2_100000ev_01.dat
 *      20180725_HD3-2_02_DARK_AgilentE3641A_34.00_AS_2_100000ev_01.dat
 *      20180725_HD3-2_02_DARK_AgilentE3641A_35.00_AS_2_100000ev_01.dat
 *      20180725_HD3-2_02_DARK_AgilentE3641A_36.00_AS_2_100000ev_01.dat
 *      20180725_HD3-2_02_DARK_AgilentE3641A_37.00_AS_2_100000ev_01.dat
 *
 *      20180725_HD3-2_03_DARK_AgilentE3641A_29.00_AS_2_100000ev_01.dat
 *      20180725_HD3-2_03_DARK_AgilentE3641A_30.00_AS_2_100000ev_01.dat
 *      20180725_HD3-2_03_DARK_AgilentE3641A_31.00_AS_2_100000ev_01.dat
 *      20180725_HD3-2_03_DARK_AgilentE3641A_32.00_AS_2_100000ev_01.dat
 *      20180725_HD3-2_03_DARK_AgilentE3641A_33.00_AS_2_100000ev_01.dat
 *      20180725_HD3-2_03_DARK_AgilentE3641A_34.00_AS_2_100000ev_01.dat
 *      20180725_HD3-2_03_DARK_AgilentE3641A_35.00_AS_2_100000ev_01.dat
 *      20180725_HD3-2_03_DARK_AgilentE3641A_36.00_AS_2_100000ev_01.dat
 *      20180725_HD3-2_03_DARK_AgilentE3641A_37.00_AS_2_100000ev_01.dat
 *
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

#define n_DCR 7

#define h 600
#define w 1000

void DCR_CrossTalk_FBK_HD3_2_from_cnt_data_2018_07();
int find_index(double v[],int N, double value);
double GetMaxOrPercentage(double num, double LimUp, double LimDown, double min_percentage, double max_percentage);

char title_DCR[80];
char title_DCR_mg[80];

void DCR_CrossTalk_FBK_HD3_2_from_cnt_data_2018_07(){

    // SiPM1:
    double HV_1[n_DCR], errHV_1[n_DCR];
    double DCR_1[n_DCR], errDCR_1[n_DCR];
    double CT_1[n_DCR], errCT_1[n_DCR];
    double DCR_Del_1[n_DCR], errDCR_Del_1[n_DCR];
    double CT_Del_1[n_DCR], errCT_Del_1[n_DCR];
    // SiPM2:
    double HV_2[n_DCR], errHV_2[n_DCR];
    double DCR_2[n_DCR], errDCR_2[n_DCR];
    double CT_2[n_DCR], errCT_2[n_DCR];
    double DCR_Del_2[n_DCR], errDCR_Del_2[n_DCR];
    double CT_Del_2[n_DCR], errCT_Del_2[n_DCR];
    // SiPM3:
    double HV_3[n_DCR], errHV_3[n_DCR];
    double DCR_3[n_DCR], errDCR_3[n_DCR];
    double CT_3[n_DCR], errCT_3[n_DCR];
    double DCR_Del_3[n_DCR], errDCR_Del_3[n_DCR];
    double CT_Del_3[n_DCR], errCT_Del_3[n_DCR];

    double HV = 0.;

    // index:
    int index = 0;

    // ERRORS:
    bool percentage_error_bool = false;
    bool fix_error_bool = false;

    double min_percentage = 0.;
    double max_percentage = 0.1;

    // EVALUARE DCR / AREA
    bool dcr_area = true;

    // DRAW ALL
    bool draw_all_bool = false;

    // Area:
    double Area = 36;

    // Initialization
    for(int i=0; i<n_DCR; i++){
        HV_1[i] = errHV_1[i] = DCR_1[i] = errDCR_1[i] = CT_1[i] = errCT_1[i] = DCR_Del_1[i] = errDCR_Del_1[i] = CT_Del_1[i] = errCT_Del_1[i] = 0.;
    }
    for(int i=0; i<n_DCR; i++){
        HV_2[i] = errHV_2[i] = DCR_2[i] = errDCR_2[i] = CT_2[i] = errCT_2[i] = DCR_Del_2[i] = errDCR_Del_2[i] = CT_Del_2[i] = errCT_Del_2[i] = 0.;
    }
    for(int i=0; i<n_DCR; i++){
        HV_3[i] = errHV_3[i] = DCR_3[i] = errDCR_3[i] = CT_3[i] = errCT_3[i] = DCR_Del_3[i] = errDCR_Del_3[i] = CT_Del_3[i] = errCT_Del_3[i] = 0.;
    }

    // HV
    HV_1[0]    = 31.00;
    errHV_1[0] =  0.01;
    for(int i=1; i<n_DCR; i++){
        HV_1[i]    = HV_1[i-1]+1.;
        errHV_1[i] = errHV_1[0];
    }

    HV_2[0]    = 31.00;
    errHV_2[0] =  0.01;
    for(int i=1; i<n_DCR; i++){
        HV_2[i]    = HV_2[i-1]+1.;
        errHV_2[i] = errHV_2[0];
    }

    HV_3[0]    = 31.00;
    errHV_3[0] =  0.01;
    for(int i=1; i<n_DCR; i++){
        HV_3[i]    = HV_3[i-1]+1.;
        errHV_3[i] = errHV_3[0];
    }

    // SiPM1:
    double HV_1_Ea[n_DCR], errHV_1_Ea[n_DCR];
    double DCR_1_Ea[n_DCR], errDCR_1_Ea[n_DCR];
    double CT_1_Ea[n_DCR], errCT_1_Ea[n_DCR];
    double DCR_Del_1_Ea[n_DCR], errDCR_Del_1_Ea[n_DCR];
    double CT_Del_1_Ea[n_DCR], errCT_Del_1_Ea[n_DCR];
    // SiPM2:
    double HV_2_Ea[n_DCR], errHV_2_Ea[n_DCR];
    double DCR_2_Ea[n_DCR], errDCR_2_Ea[n_DCR];
    double CT_2_Ea[n_DCR], errCT_2_Ea[n_DCR];
    double DCR_Del_2_Ea[n_DCR], errDCR_Del_2_Ea[n_DCR];
    double CT_Del_2_Ea[n_DCR], errCT_Del_2_Ea[n_DCR];
    // SiPM3:
    double HV_3_Ea[n_DCR], errHV_3_Ea[n_DCR];
    double DCR_3_Ea[n_DCR], errDCR_3_Ea[n_DCR];
    double CT_3_Ea[n_DCR], errCT_3_Ea[n_DCR];
    double DCR_Del_3_Ea[n_DCR], errDCR_Del_3_Ea[n_DCR];
    double CT_Del_3_Ea[n_DCR], errCT_Del_3_Ea[n_DCR];


    // Initialization
    for(int i=0; i<n_DCR; i++){
        HV_1_Ea[i] = errHV_1_Ea[i] = DCR_1_Ea[i] = errDCR_1_Ea[i] = CT_1_Ea[i] = errCT_1_Ea[i] = DCR_Del_1_Ea[i] = errDCR_Del_1_Ea[i] = CT_Del_1_Ea[i] = errCT_Del_1_Ea[i] = 0.;
    }
    for(int i=0; i<n_DCR; i++){
        HV_2_Ea[i] = errHV_2_Ea[i] = DCR_2_Ea[i] = errDCR_2_Ea[i] = CT_2_Ea[i] = errCT_2_Ea[i] = DCR_Del_2_Ea[i] = errDCR_Del_2_Ea[i] = CT_Del_2_Ea[i] = errCT_Del_2_Ea[i] = 0.;
    }
    for(int i=0; i<n_DCR; i++){
        HV_3_Ea[i] = errHV_3_Ea[i] = DCR_3_Ea[i] = errDCR_3_Ea[i] = CT_3_Ea[i] = errCT_3_Ea[i] = DCR_Del_3_Ea[i] = errDCR_Del_3_Ea[i] = CT_Del_3_Ea[i] = errCT_Del_3_Ea[i] = 0.;
    }

    // HV
    HV_1_Ea[0]    = 31.00;
    errHV_1_Ea[0] =  0.01;
    for(int i=1; i<n_DCR; i++){
        HV_1_Ea[i]    = HV_1_Ea[i-1]+1.;
        errHV_1_Ea[i] = errHV_1_Ea[0];
    }

    HV_2_Ea[0]    = 31.00;
    errHV_2_Ea[0] =  0.01;
    for(int i=1; i<n_DCR; i++){
        HV_2_Ea[i]    = HV_2_Ea[i-1]+1.;
        errHV_2_Ea[i] = errHV_2_Ea[0];
    }

    HV_3_Ea[0]    = 31.00;
    errHV_3_Ea[0] =  0.01;
    for(int i=1; i<n_DCR; i++){
        HV_3_Ea[i]    = HV_3_Ea[i-1]+1.;
        errHV_3_Ea[i] = errHV_3_Ea[0];
    }

    // SiPM1:
    double HV_1_Eb[n_DCR], errHV_1_Eb[n_DCR];
    double DCR_1_Eb[n_DCR], errDCR_1_Eb[n_DCR];
    double CT_1_Eb[n_DCR], errCT_1_Eb[n_DCR];
    double DCR_Del_1_Eb[n_DCR], errDCR_Del_1_Eb[n_DCR];
    double CT_Del_1_Eb[n_DCR], errCT_Del_1_Eb[n_DCR];
    // SiPM2:
    double HV_2_Eb[n_DCR], errHV_2_Eb[n_DCR];
    double DCR_2_Eb[n_DCR], errDCR_2_Eb[n_DCR];
    double CT_2_Eb[n_DCR], errCT_2_Eb[n_DCR];
    double DCR_Del_2_Eb[n_DCR], errDCR_Del_2_Eb[n_DCR];
    double CT_Del_2_Eb[n_DCR], errCT_Del_2_Eb[n_DCR];
    // SiPM3:
    double HV_3_Eb[n_DCR], errHV_3_Eb[n_DCR];
    double DCR_3_Eb[n_DCR], errDCR_3_Eb[n_DCR];
    double CT_3_Eb[n_DCR], errCT_3_Eb[n_DCR];
    double DCR_Del_3_Eb[n_DCR], errDCR_Del_3_Eb[n_DCR];
    double CT_Del_3_Eb[n_DCR], errCT_Del_3_Eb[n_DCR];


    // Initialization
    for(int i=0; i<n_DCR; i++){
        HV_1_Eb[i] = errHV_1_Eb[i] = DCR_1_Eb[i] = errDCR_1_Eb[i] = CT_1_Eb[i] = errCT_1_Eb[i] = DCR_Del_1_Eb[i] = errDCR_Del_1_Eb[i] = CT_Del_1_Eb[i] = errCT_Del_1_Eb[i] = 0.;
    }
    for(int i=0; i<n_DCR; i++){
        HV_2_Eb[i] = errHV_2_Eb[i] = DCR_2_Eb[i] = errDCR_2_Eb[i] = CT_2_Eb[i] = errCT_2_Eb[i] = DCR_Del_2_Eb[i] = errDCR_Del_2_Eb[i] = CT_Del_2_Eb[i] = errCT_Del_2_Eb[i] = 0.;
    }
    for(int i=0; i<n_DCR; i++){
        HV_3_Eb[i] = errHV_3_Eb[i] = DCR_3_Eb[i] = errDCR_3_Eb[i] = CT_3_Eb[i] = errCT_3_Eb[i] = DCR_Del_3_Eb[i] = errDCR_Del_3_Eb[i] = CT_Del_3_Eb[i] = errCT_Del_3_Eb[i] = 0.;
    }

    // HV
    HV_1_Eb[0]    = 31.00;
    errHV_1_Eb[0] =  0.01;
    for(int i=1; i<n_DCR; i++){
        HV_1_Eb[i]    = HV_1_Eb[i-1]+1.;
        errHV_1_Eb[i] = errHV_1_Eb[0];
    }

    HV_2_Eb[0]    = 31.00;
    errHV_2_Eb[0] =  0.01;
    for(int i=1; i<n_DCR; i++){
        HV_2_Eb[i]    = HV_2_Eb[i-1]+1.;
        errHV_2_Eb[i] = errHV_2_Eb[0];
    }

    HV_3_Eb[0]    = 31.00;
    errHV_3_Eb[0] =  0.01;
    for(int i=1; i<n_DCR; i++){
        HV_3_Eb[i]    = HV_3_Eb[i-1]+1.;
        errHV_3_Eb[i] = errHV_3_Eb[0];
    }

    //--------------------------------------------------------------------------
    // For errors
    // double diff_temp_0 = 1;
    // double diff_temp_1 = 1*(1+3./(double)(nfiletot-1)*(double)i);
    // diff_temp_1 = (int)(diff_temp_1*10);
    // diff_temp_1/=10;
    //--------------------------------------------------------------------------



    ///////////////////////////////////////////////////////////////////////////
    //      SiPM1
    ///////////////////////////////////////////////////////////////////////////
    // pe_0_5_vect[7] = {8, 8, 8, 9, 9, 9, 10}
    // pe_1_5_vect[7] = {19, 23, 26, 30, 35, 37, 40}
    //minyhistDelays = 20;  maxyhistDelays = 150
    //expDelLow_max = 20;  expDelHigh_max = 150
    HV = 31.00;
    index = find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV);
    DCR_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 10.8468;
    errDCR_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.0110521;
    CT_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.194892;
    DCR_Del_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 10.6742;
    errDCR_Del_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.0332437;
    CT_Del_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.206139;

    HV = 32.00;
    index = find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV);
    DCR_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 14.7845;
    errDCR_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.0131704;
    CT_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.20516;
    DCR_Del_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 15.1647;
    errDCR_Del_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.0290423;
    CT_Del_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.226054;

    HV = 33.00;
    index = find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV);
    DCR_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 17.7625;
    errDCR_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.0146536;
    CT_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.240606;
    DCR_Del_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 18.458;
    errDCR_Del_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.0277113;
    CT_Del_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.263384;

    HV = 34.00;
    index = find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV);
    DCR_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 20.2909;
    errDCR_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.0158566;
    CT_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.269749;
    DCR_Del_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 21.1589;
    errDCR_Del_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.0271707;
    CT_Del_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.287259;

    HV = 35.00;
    index = find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV);
    DCR_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 22.7001;
    errDCR_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.0169655;
    CT_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.29342;
    DCR_Del_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 23.717;
    errDCR_Del_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.0269739;
    CT_Del_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.306948;

    HV = 36.00;
    index = find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV);
    DCR_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 24.9459;
    errDCR_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.0179725;
    CT_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.328475;
    DCR_Del_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 26.0245;
    errDCR_Del_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.0269753;
    CT_Del_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.340825;

    HV = 37.00;
    index = find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV);
    DCR_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 27.2241;
    errDCR_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.0189719;
    CT_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.358377;
    DCR_Del_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 28.4227;
    errDCR_Del_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.02718;
    CT_Del_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.366465;


    //////// Ea ////////
    HV = 31.00;
    index = find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV);
    DCR_1_Ea[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)] = DCR_1[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)];
    errDCR_1_Ea[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)] = errDCR_1[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)];
    CT_1_Ea[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)] = CT_1[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)];
    DCR_Del_1_Ea[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)] = DCR_Del_1[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)];
    errDCR_Del_1_Ea[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)] = errDCR_Del_1[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)];
    CT_Del_1_Ea[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)] = CT_Del_1[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)];

    HV = 32.00;
    index = find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV);
    DCR_1_Ea[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)] = DCR_1[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)];
    errDCR_1_Ea[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)] = errDCR_1[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)];
    CT_1_Ea[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)] = CT_1[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)];
    DCR_Del_1_Ea[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)] = DCR_Del_1[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)];
    errDCR_Del_1_Ea[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)] = errDCR_Del_1[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)];
    CT_Del_1_Ea[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)] = CT_Del_1[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)];

    HV = 33.00;
    index = find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV);
    DCR_1_Ea[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)] = DCR_1[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)];
    errDCR_1_Ea[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)] = errDCR_1[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)];
    CT_1_Ea[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)] = CT_1[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)];
    DCR_Del_1_Ea[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)] = DCR_Del_1[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)];
    errDCR_Del_1_Ea[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)] = errDCR_Del_1[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)];
    CT_Del_1_Ea[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)] = CT_Del_1[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)];

    HV = 34.00;
    index = find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV);
    DCR_1_Ea[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)] = 20.3256;
    errDCR_1_Ea[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)] = 0.0158728;
    CT_1_Ea[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)] = 0.281489;
    DCR_Del_1_Ea[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)] = 21.186;
    errDCR_Del_1_Ea[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)] = 0.027144;
    CT_Del_1_Ea[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)] = 0.297041;

    HV = 35.00;
    index = find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV);
    DCR_1_Ea[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)] = 22.7297;
    errDCR_1_Ea[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)] = 0.0169789;
    CT_1_Ea[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)] = 0.303196;
    DCR_Del_1_Ea[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)] = 23.741;
    errDCR_Del_1_Ea[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)] = 0.0269559;
    CT_Del_1_Ea[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)] = 0.321314;

    HV = 36.00;
    index = find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV);
    DCR_1_Ea[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)] = 24.9651;
    errDCR_1_Ea[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)] = 0.017981;
    CT_1_Ea[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)] = 0.33856;
    DCR_Del_1_Ea[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)] = 26.0412;
    errDCR_Del_1_Ea[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)] = 0.0269684;
    CT_Del_1_Ea[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)] = 0.352023;

    HV = 37.00;
    index = find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV);
    DCR_1_Ea[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)] = 27.2669;
    errDCR_1_Ea[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)] = 0.0189905;
    CT_1_Ea[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)] = 0.368358;
    DCR_Del_1_Ea[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)] = 28.4689;
    errDCR_Del_1_Ea[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)] = 0.0271892;
    CT_Del_1_Ea[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)] = 0.378462;


    //////// Eb ////////
    HV = 31.00;
    index = find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV);
    DCR_1_Eb[find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV)] = 10.7205;
    errDCR_1_Eb[find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV)] = 0.0109802;
    CT_1_Eb[find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV)] = 0.178874;
    DCR_Del_1_Eb[find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV)] = 10.5428;
    errDCR_Del_1_Eb[find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV)] = 0.033473;
    CT_Del_1_Eb[find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV)] = 0.189372;

    HV = 32.00;
    index = find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV);
    DCR_1_Eb[find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV)] = 14.7069;
    errDCR_1_Eb[find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV)] = 0.0131305;
    CT_1_Eb[find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV)] = 0.192051;
    DCR_Del_1_Eb[find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV)] = 15.0977;
    errDCR_Del_1_Eb[find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV)] = 0.0291175;
    CT_Del_1_Eb[find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV)] = 0.20332;

    HV = 33.00;
    index = find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV);
    DCR_1_Eb[find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV)] = 17.7188;
    errDCR_1_Eb[find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV)] = 0.0146324;
    CT_1_Eb[find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV)] = 0.229271;
    DCR_Del_1_Eb[find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV)] = 18.4259;
    errDCR_Del_1_Eb[find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV)] = 0.0277465;
    CT_Del_1_Eb[find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV)] = 0.248407;

    HV = 34.00;
    index = find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV);
    DCR_1_Eb[find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV)] = 20.2242;
    errDCR_1_Eb[find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV)] = 0.0158254;
    CT_1_Eb[find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV)] = 0.259681;
    DCR_Del_1_Eb[find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV)] = 21.1079;
    errDCR_Del_1_Eb[find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV)] = 0.0272029;
    CT_Del_1_Eb[find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV)] = 0.268703;

    HV = 35.00;
    index = find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV);
    DCR_1_Eb[find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV)] = 22.6376;
    errDCR_1_Eb[find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV)] = 0.0169372;
    CT_1_Eb[find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV)] = 0.279383;
    DCR_Del_1_Eb[find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV)] = 23.6671;
    errDCR_Del_1_Eb[find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV)] = 0.0269916;
    CT_Del_1_Eb[find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV)] = 0.280708;

    HV = 36.00;
    index = find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV);
    DCR_1_Eb[find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV)] = 24.8999;
    errDCR_1_Eb[find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV)] = 0.0179521;
    CT_1_Eb[find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV)] = 0.316264;
    DCR_Del_1_Eb[find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV)] = 25.9781;
    errDCR_Del_1_Eb[find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV)] = 0.0269812;
    CT_Del_1_Eb[find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV)] = 0.322032;

    HV = 37.00;
    index = find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV);
    DCR_1_Eb[find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV)] = 27.1703;
    errDCR_1_Eb[find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV)] = 0.0189485;
    CT_1_Eb[find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV)] = 0.345828;
    DCR_Del_1_Eb[find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV)] = 28.3609;
    errDCR_Del_1_Eb[find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV)] = 0.0271687;
    CT_Del_1_Eb[find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV)] = 0.34792;






    ///////////////////////////////////////////////////////////////////////////
    //      SiPM2
    ///////////////////////////////////////////////////////////////////////////
    // pe_0_5_vect[7] = {8, 9, 9, 10, 10, 10, 11}
    // pe_1_5_vect[7] = {19, 22, 25, 29, 33, 34, 38}
    //minyhistDelays = 20;  maxyhistDelays = 150
    //expDelLow_max = 20;  expDelHigh_max = 150
    HV = 31.00;
    index = find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV);
    DCR_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 10.9316;
    errDCR_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.0111002;
    CT_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.176788;
    DCR_Del_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 10.7255;
    errDCR_Del_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.0331164;
    CT_Del_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.19483;

    HV = 32.00;
    index = find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV);
    DCR_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 14.9826;
    errDCR_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.0132717;
    CT_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.206578;
    DCR_Del_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 15.2545;
    errDCR_Del_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.0286105;
    CT_Del_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.227325;

    HV = 33.00;
    index = find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV);
    DCR_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 18.7346;
    errDCR_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.0151214;
    CT_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.236805;
    DCR_Del_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 19.3657;
    errDCR_Del_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.0270007;
    CT_Del_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.257828;

    HV = 34.00;
    index = find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV);
    DCR_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 21.4701;
    errDCR_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.0164034;
    CT_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.263922;
    DCR_Del_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 22.3912;
    errDCR_Del_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.0266961;
    CT_Del_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.277089;

    HV = 35.00;
    index = find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV);
    DCR_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 24.1236;
    errDCR_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.0176065;
    CT_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.289775;
    DCR_Del_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 25.1994;
    errDCR_Del_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.0266617;
    CT_Del_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.300092;

    HV = 36.00;
    index = find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV);
    DCR_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 26.6274;
    errDCR_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.0187121;
    CT_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.328639;
    DCR_Del_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 27.837;
    errDCR_Del_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.0268445;
    CT_Del_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.337381;

    HV = 37.00;
    index = find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV);
    DCR_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 28.9859;
    errDCR_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.0197316;
    CT_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.354963;
    DCR_Del_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 30.3086;
    errDCR_Del_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.0271966;
    CT_Del_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.35954;


    //////// Ea ////////
    HV = 31.00;
    index = find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV);
    DCR_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = DCR_2[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)];
    errDCR_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = errDCR_2[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)];
    CT_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = CT_2[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)];
    DCR_Del_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = DCR_Del_2[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)];
    errDCR_Del_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = errDCR_Del_2[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)];
    CT_Del_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = CT_Del_2[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)];

    HV = 32.00;
    index = find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV);
    DCR_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 15.1275;
    errDCR_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 0.0133455;
    CT_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 0.2284;
    DCR_Del_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 15.4069;
    errDCR_Del_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 0.0284775;
    CT_Del_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 0.247424;

    HV = 33.00;
    index = find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV);
    DCR_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 18.8724;
    errDCR_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 0.0151872;
    CT_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 0.259608;
    DCR_Del_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 19.4985;
    errDCR_Del_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 0.0269203;
    CT_Del_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 0.27471;

    HV = 34.00;
    index = find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV);
    DCR_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 21.6926;
    errDCR_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 0.0165057;
    CT_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 0.279133;
    DCR_Del_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 22.601;
    errDCR_Del_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 0.0266066;
    CT_Del_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 0.292766;

    HV = 35.00;
    index = find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV);
    DCR_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 24.3031;
    errDCR_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 0.0176867;
    CT_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 0.302267;
    DCR_Del_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 25.3718;
    errDCR_Del_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 0.0265973;
    CT_Del_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 0.314508;

    HV = 36.00;
    index = find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV);
    DCR_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 26.7797;
    errDCR_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 0.0187785;
    CT_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 0.347915;
    DCR_Del_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 27.9802;
    errDCR_Del_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 0.026794;
    CT_Del_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 0.350269;

    HV = 37.00;
    index = find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV);
    DCR_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 29.0874;
    errDCR_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 0.0197751;
    CT_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 0.36895;
    DCR_Del_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 30.4092;
    errDCR_Del_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 0.0271947;
    CT_Del_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 0.373502;


    //////// Eb ////////
    HV = 31.00;
    index = find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV);
    DCR_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 10.6895;
    errDCR_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 0.0109625;
    CT_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 0.157942;
    DCR_Del_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 10.4833;
    errDCR_Del_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 0.0335342;
    CT_Del_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 0.176065;

    HV = 32.00;
    index = find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV);
    DCR_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 14.5568;
    errDCR_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 0.0130533;
    CT_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 0.189521;
    DCR_Del_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 14.8065;
    errDCR_Del_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 0.0289394;
    CT_Del_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 0.203812;

    HV = 33.00;
    index = find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV);
    DCR_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 18.4081;
    errDCR_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 0.014965;
    CT_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 0.221794;
    DCR_Del_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 19.0392;
    errDCR_Del_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 0.0271436;
    CT_Del_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 0.233831;

    HV = 34.00;
    index = find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV);
    DCR_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 21.2632;
    errDCR_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 0.0163081;
    CT_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 0.247879;
    DCR_Del_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 22.1801;
    errDCR_Del_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 0.0267407;
    CT_Del_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 0.251234;

    HV = 35.00;
    index = find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV);
    DCR_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 23.9642;
    errDCR_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 0.0175352;
    CT_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 0.271288;
    DCR_Del_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 25.0477;
    errDCR_Del_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 0.0266949;
    CT_Del_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 0.269809;

    HV = 36.00;
    index = find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV);
    DCR_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 26.4917;
    errDCR_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 0.0186528;
    CT_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 0.31498;
    DCR_Del_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 27.7111;
    errDCR_Del_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 0.0268662;
    CT_Del_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 0.318684;

    HV = 37.00;
    index = find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV);
    DCR_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 28.875;
    errDCR_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 0.0196841;
    CT_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 0.339219;
    DCR_Del_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 30.1893;
    errDCR_Del_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 0.0271704;
    CT_Del_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 0.337829;


    ////////////////////////////////////////////////////////////////////////////
    //      SiPM3
    ////////////////////////////////////////////////////////////////////////////
    // pe_0_5_vect[7] = {8, 8, 9, 10, 11, 12, 14}
    // pe_1_5_vect[7] = {19.5, 23.5, 26.5, 30, 33.5, 36, 40}
    //minyhistDelays = 20;  maxyhistDelays = 150
    //expDelLow_max = 20;  expDelHigh_max = 150
    HV = 31.00;
    index = find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV);
    DCR_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 13.2422;
    errDCR_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 0.012366;
    CT_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 0.192589;
    DCR_Del_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 13.1078;
    errDCR_Del_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 0.0298438;
    CT_Del_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 0.189957;

    HV = 32.00;
    index = find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV);
    DCR_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 17.6465;
    errDCR_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 0.0145972;
    CT_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 0.208325;
    DCR_Del_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 18.0967;
    errDCR_Del_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 0.0271804;
    CT_Del_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 0.222822;

    HV = 33.00;
    index = find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV);
    DCR_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 21.2145;
    errDCR_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 0.0162856;
    CT_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 0.24549;
    DCR_Del_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 22.0211;
    errDCR_Del_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 0.0266151;
    CT_Del_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 0.26071;

    HV = 34.00;
    index = find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV);
    DCR_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 24.272;
    errDCR_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 0.0176727;
    CT_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 0.276746;
    DCR_Del_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 25.3283;
    errDCR_Del_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 0.0266077;
    CT_Del_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 0.286266;

    HV = 35.00;
    index = find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV);
    DCR_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 27.1553;
    errDCR_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 0.018942;
    CT_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 0.306715;
    DCR_Del_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 28.3349;
    errDCR_Del_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 0.0267659;
    CT_Del_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 0.310631;

    HV = 36.00;
    index = find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV);
    DCR_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 29.8447;
    errDCR_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 0.0200983;
    CT_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 0.339257;
    DCR_Del_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 31.2002;
    errDCR_Del_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 0.0271952;
    CT_Del_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 0.341255;

    HV = 37.00;
    index = find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV);
    DCR_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 32.3161;
    errDCR_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 0.0211412;
    CT_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 0.369078;
    DCR_Del_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 33.8594;
    errDCR_Del_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 0.027641;
    CT_Del_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 0.366327;


    //////// Ea ////////
    HV = 31.00;
    index = find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV);
    DCR_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = DCR_3[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)];
    errDCR_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = errDCR_3[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)];
    CT_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = CT_3[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)];
    DCR_Del_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = DCR_Del_3[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)];
    errDCR_Del_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = errDCR_Del_3[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)];
    CT_Del_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = CT_Del_3[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)];

    HV = 32.00;
    index = find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV);
    DCR_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = DCR_3[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)];
    errDCR_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = errDCR_3[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)];
    CT_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = CT_3[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)];
    DCR_Del_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = DCR_Del_3[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)];
    errDCR_Del_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = errDCR_Del_3[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)];
    CT_Del_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = CT_Del_3[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)];

    HV = 33.00;
    index = find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV);
    DCR_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = 21.2565;
    errDCR_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = 0.016305;
    CT_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = 0.259796;
    DCR_Del_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = 22.0518;
    errDCR_Del_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = 0.0265947;
    CT_Del_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = 0.269417;

    HV = 34.00;
    index = find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV);
    DCR_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = 24.3413;
    errDCR_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = 0.0177037;
    CT_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = 0.288257;
    DCR_Del_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = 25.3939;
    errDCR_Del_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = 0.0266003;
    CT_Del_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = 0.297125;

    HV = 35.00;
    index = find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV);
    DCR_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = 27.2304;
    errDCR_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = 0.0189746;
    CT_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = 0.317317;
    DCR_Del_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = 28.4193;
    errDCR_Del_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = 0.0267863;
    CT_Del_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = 0.322109;

    HV = 36.00;
    index = find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV);
    DCR_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = 29.9533;
    errDCR_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = 0.0201445;
    CT_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = 0.350983;
    DCR_Del_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = 31.3187;
    errDCR_Del_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = 0.0272407;
    CT_Del_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = 0.351507;

    HV = 37.00;
    index = find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV);
    DCR_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = 32.4947;
    errDCR_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = 0.0212159;
    CT_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = 0.379357;
    DCR_Del_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = 34.0568;
    errDCR_Del_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = 0.027737;
    CT_Del_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = 0.378994;


    //////// Eb ////////
    HV = 31.00;
    index = find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV);
    DCR_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 13.1021;
    errDCR_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 0.0122915;
    CT_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 0.176183;
    DCR_Del_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 12.9489;
    errDCR_Del_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 0.0299836;
    CT_Del_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 0.171953;

    HV = 32.00;
    index = find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV);
    DCR_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 17.5785;
    errDCR_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 0.0145641;
    CT_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 0.194778;
    DCR_Del_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 18.0338;
    errDCR_Del_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 0.0272164;
    CT_Del_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 0.194892;

    HV = 33.00;
    index = find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV);
    DCR_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 21.1104;
    errDCR_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 0.0162375;
    CT_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 0.23415;
    DCR_Del_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 21.9188;
    errDCR_Del_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 0.0266291;
    CT_Del_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 0.240303;

    HV = 34.00;
    index = find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV);
    DCR_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 24.1866;
    errDCR_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 0.0176346;
    CT_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 0.26584;
    DCR_Del_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 25.2341;
    errDCR_Del_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 0.0265887;
    CT_Del_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 0.267889;

    HV = 35.00;
    index = find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV);
    DCR_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 27.0406;
    errDCR_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 0.0188921;
    CT_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 0.294349;
    DCR_Del_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 28.2044;
    errDCR_Del_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 0.0267136;
    CT_Del_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 0.291051;

    HV = 36.00;
    index = find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV);
    DCR_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 29.7127;
    errDCR_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 0.0200421;
    CT_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 0.327356;
    DCR_Del_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 31.0564;
    errDCR_Del_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 0.0271189;
    CT_Del_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 0.32386;

    HV = 37.00;
    index = find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV);
    DCR_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 32.0933;
    errDCR_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 0.0210479;
    CT_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 0.355082;
    DCR_Del_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 33.6209;
    errDCR_Del_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 0.0274934;
    CT_Del_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 0.345098;

    //------------------------------

    // FIXED ERROR
    if(fix_error_bool){ // fix_error_bool
        double err_fix_DCR = 1;
        double err_fix_CT = 0.02;
        for(int i=0; i<n_DCR; i++){
            errDCR_1[i] = err_fix_DCR;
        }
        for(int i=0; i<n_DCR; i++){
            errDCR_2[i] = err_fix_DCR;
        }
        for(int i=0; i<n_DCR; i++){
            errDCR_3[i] = err_fix_DCR;
        }
        for(int i=0; i<n_DCR; i++){
            errCT_1[i] = err_fix_CT;
        }
        for(int i=0; i<n_DCR; i++){
            errCT_2[i] = err_fix_CT;
        }
        for(int i=0; i<n_DCR; i++){
            errCT_3[i] = err_fix_CT;
        }
        for(int i=0; i<n_DCR; i++){
            errDCR_Del_1[i] = err_fix_DCR;
        }
        for(int i=0; i<n_DCR; i++){
            errDCR_Del_2[i] = err_fix_DCR;
        }
        for(int i=0; i<n_DCR; i++){
            errDCR_Del_3[i] = err_fix_DCR;
        }
        for(int i=0; i<n_DCR; i++){
            errCT_Del_1[i] = err_fix_CT;
        }
        for(int i=0; i<n_DCR; i++){
            errCT_Del_2[i] = err_fix_CT;
        }
        for(int i=0; i<n_DCR; i++){
            errCT_Del_3[i] = err_fix_CT;
        }
    } // end fix_error_bool
    // error from vect
    for(int i=0; i<n_DCR; i++){
        errDCR_1[i] = GetMaxOrPercentage(DCR_1[i], DCR_1_Ea[i], DCR_1_Eb[i], min_percentage, max_percentage);
        errDCR_2[i] = GetMaxOrPercentage(DCR_2[i], DCR_2_Ea[i], DCR_2_Eb[i], min_percentage, max_percentage);
        errDCR_3[i] = GetMaxOrPercentage(DCR_3[i], DCR_3_Ea[i], DCR_3_Eb[i], min_percentage, max_percentage);
        errDCR_Del_1[i] = GetMaxOrPercentage(DCR_Del_1[i], DCR_Del_1_Ea[i], DCR_Del_1_Eb[i], min_percentage, max_percentage);
        errDCR_Del_2[i] = GetMaxOrPercentage(DCR_Del_2[i], DCR_Del_2_Ea[i], DCR_Del_2_Eb[i], min_percentage, max_percentage);
        errDCR_Del_3[i] = GetMaxOrPercentage(DCR_Del_3[i], DCR_Del_3_Ea[i], DCR_Del_3_Eb[i], min_percentage, max_percentage);
        errCT_1[i] = GetMaxOrPercentage(CT_1[i], CT_1_Ea[i], CT_1_Eb[i], min_percentage, max_percentage);
        errCT_2[i] = GetMaxOrPercentage(CT_2[i], CT_2_Ea[i], CT_2_Eb[i], min_percentage, max_percentage);
        errCT_3[i] = GetMaxOrPercentage(CT_3[i], CT_3_Ea[i], CT_3_Eb[i], min_percentage, max_percentage);
        errCT_Del_1[i] = GetMaxOrPercentage(CT_Del_1[i], CT_Del_1_Ea[i], CT_Del_1_Eb[i], min_percentage, max_percentage);
        errCT_Del_2[i] = GetMaxOrPercentage(CT_Del_2[i], CT_Del_2_Ea[i], CT_Del_2_Eb[i], min_percentage, max_percentage);
        errCT_Del_3[i] = GetMaxOrPercentage(CT_Del_3[i], CT_Del_3_Ea[i], CT_Del_3_Eb[i], min_percentage, max_percentage);
    } // end error from vect


    //------------------------------

    if(dcr_area){
        for(int i=0; i< n_DCR; i++){
            DCR_1[i]  /= Area;
            errDCR_1[i] /= Area;
            DCR_1[i]  *= 1e3;
            errDCR_1[i] *= 1e3;
        }
        for(int i=0; i< n_DCR; i++){
            DCR_2[i]  /= Area;
            errDCR_2[i] /= Area;
            DCR_2[i]  *= 1e3;
            errDCR_2[i] *= 1e3;
        }
        for(int i=0; i< n_DCR; i++){
            DCR_3[i]  /= Area;
            errDCR_3[i] /= Area;
            DCR_3[i]  *= 1e3;
            errDCR_3[i] *= 1e3;
        }

        //------------------------------

        for(int i=0; i< n_DCR; i++){
            DCR_Del_1[i]  /= Area;
            errDCR_Del_1[i] /= Area;
            DCR_Del_1[i]  *= 1e3;
            errDCR_Del_1[i] *= 1e3;
        }
        for(int i=0; i< n_DCR; i++){
            DCR_Del_2[i]  /= Area;
            errDCR_Del_2[i] /= Area;
            DCR_Del_2[i]  *= 1e3;
            errDCR_Del_2[i] *= 1e3;
        }
        for(int i=0; i< n_DCR; i++){
            DCR_Del_3[i]  /= Area;
            errDCR_Del_3[i] /= Area;
            DCR_Del_3[i]  *= 1e3;
            errDCR_Del_3[i] *= 1e3;
        }

    }


    //------------------------------
    //------------------------------


    TGraphErrors *gV_DCR_1  = new TGraphErrors(n_DCR, HV_1, DCR_1, errHV_1, errDCR_1);
    TGraphErrors *gV_DCR_2  = new TGraphErrors(n_DCR, HV_2, DCR_2, errHV_2, errDCR_2);
    TGraphErrors *gV_DCR_3  = new TGraphErrors(n_DCR, HV_3, DCR_3, errHV_3, errDCR_3);

    TGraphErrors *gV_CT_1  = new TGraphErrors(n_DCR, HV_1, CT_1, errHV_1, errCT_1);
    TGraphErrors *gV_CT_2  = new TGraphErrors(n_DCR, HV_2, CT_2, errHV_2, errCT_2);
    TGraphErrors *gV_CT_3  = new TGraphErrors(n_DCR, HV_3, CT_3, errHV_3, errCT_3);

    TGraphErrors *gV_DCR_Del_1  = new TGraphErrors(n_DCR, HV_1, DCR_Del_1, errHV_1, errDCR_Del_1);
    TGraphErrors *gV_DCR_Del_2  = new TGraphErrors(n_DCR, HV_2, DCR_Del_2, errHV_2, errDCR_Del_2);
    TGraphErrors *gV_DCR_Del_3  = new TGraphErrors(n_DCR, HV_3, DCR_Del_3, errHV_3, errDCR_3);

    TGraphErrors *gV_CT_Del_1  = new TGraphErrors(n_DCR, HV_1, CT_Del_1, errHV_1, errCT_Del_1);
    TGraphErrors *gV_CT_Del_2  = new TGraphErrors(n_DCR, HV_2, CT_Del_2, errHV_2, errCT_Del_2);
    TGraphErrors *gV_CT_Del_3  = new TGraphErrors(n_DCR, HV_3, CT_Del_3, errHV_3, errCT_Del_3);


    //------------------------------

    if(dcr_area) strcpy(title_DCR, "\\frac{DCR}{Area} \\left(\\frac{kHz}{mm^2}\\right)");
    else         strcpy(title_DCR, "DCR (MHz)");

    strcpy(title_DCR_mg, ";Bias Voltage (V); ");
    strcat(title_DCR_mg, title_DCR);

    gV_DCR_1->SetMarkerStyle(20);
    gV_DCR_1->SetMarkerColor(kOrange+1);
    gV_DCR_1->SetTitle();
    gV_DCR_1->GetXaxis()->SetTitle("Bias Voltage (V)");
    gV_DCR_1->GetYaxis()->SetTitle(title_DCR);

    gV_DCR_2->SetMarkerStyle(20);
    gV_DCR_2->SetMarkerColor(kRed);
    gV_DCR_2->SetTitle();
    gV_DCR_2->GetXaxis()->SetTitle("Bias Voltage (V)");
    gV_DCR_2->GetYaxis()->SetTitle(title_DCR);

    gV_DCR_3->SetMarkerStyle(20);
    gV_DCR_3->SetMarkerColor(kMagenta);
    gV_DCR_3->SetTitle();
    gV_DCR_3->GetXaxis()->SetTitle("Bias Voltage (V)");
    gV_DCR_3->GetYaxis()->SetTitle(title_DCR);

    //------------------------------

    gV_CT_1->SetMarkerStyle(20);
    gV_CT_1->SetMarkerColor(kOrange+1);
    gV_CT_1->SetTitle();
    gV_CT_1->GetXaxis()->SetTitle("Bias Voltage (V)");
    gV_CT_1->GetYaxis()->SetTitle("P_{CT}");

    gV_CT_2->SetMarkerStyle(20);
    gV_CT_2->SetMarkerColor(kRed);
    gV_CT_2->SetTitle();
    gV_CT_2->GetXaxis()->SetTitle("Bias Voltage (V)");
    gV_CT_2->GetYaxis()->SetTitle("P_{CT}");

    gV_CT_3->SetMarkerStyle(20);
    gV_CT_3->SetMarkerColor(kMagenta);
    gV_CT_3->SetTitle();
    gV_CT_3->GetXaxis()->SetTitle("Bias Voltage (V)");
    gV_CT_3->GetYaxis()->SetTitle("P_{CT}");

    //------------------------------

    //------------------------------

    gV_DCR_Del_1->SetMarkerStyle(22);
    gV_DCR_Del_1->SetMarkerColor(kOrange+1);
    gV_DCR_Del_1->SetTitle();
    gV_DCR_Del_1->GetXaxis()->SetTitle("Bias Voltage (V)");
    gV_DCR_Del_1->GetYaxis()->SetTitle(title_DCR);

    gV_DCR_Del_2->SetMarkerStyle(22);
    gV_DCR_Del_2->SetMarkerColor(kRed);
    gV_DCR_Del_2->SetTitle();
    gV_DCR_Del_2->GetXaxis()->SetTitle("Bias Voltage (V)");
    gV_DCR_Del_2->GetYaxis()->SetTitle(title_DCR);

    gV_DCR_Del_3->SetMarkerStyle(22);
    gV_DCR_Del_3->SetMarkerColor(kMagenta);
    gV_DCR_Del_3->SetTitle();
    gV_DCR_Del_3->GetXaxis()->SetTitle("Bias Voltage (V)");
    gV_DCR_Del_3->GetYaxis()->SetTitle(title_DCR);

    //------------------------------

    gV_CT_Del_1->SetMarkerStyle(22);
    gV_CT_Del_1->SetMarkerColor(kOrange+1);
    gV_CT_Del_1->SetTitle();
    gV_CT_Del_1->GetXaxis()->SetTitle("Bias Voltage (V)");
    gV_CT_Del_1->GetYaxis()->SetTitle("P_{CT}");

    gV_CT_Del_2->SetMarkerStyle(22);
    gV_CT_Del_2->SetMarkerColor(kRed);
    gV_CT_Del_2->SetTitle();
    gV_CT_Del_2->GetXaxis()->SetTitle("Bias Voltage (V)");
    gV_CT_Del_2->GetYaxis()->SetTitle("P_{CT}");

    gV_CT_Del_3->SetMarkerStyle(22);
    gV_CT_Del_3->SetMarkerColor(kMagenta);
    gV_CT_Del_3->SetTitle();
    gV_CT_Del_3->GetXaxis()->SetTitle("Bias Voltage (V)");
    gV_CT_Del_3->GetYaxis()->SetTitle("P_{CT}");

    //------------------------------

    if(draw_all_bool){
        TCanvas *cV_DCR_1 = new TCanvas("cV_DCR_1", "cV_DCR_1",w,h);
        cV_DCR_1->SetGrid();
        gV_DCR_1->Draw("AP");

        TCanvas *cV_DCR_2 = new TCanvas("cV_DCR_2", "cV_DCR_2",w,h);
        cV_DCR_2->SetGrid();
        gV_DCR_2->Draw("AP");

        TCanvas *cV_DCR_3 = new TCanvas("cV_DCR_3", "cV_DCR_3",w,h);
        cV_DCR_3->SetGrid();
        gV_DCR_3->Draw("AP");
    }

    //------------------------------

    if(draw_all_bool){
        TCanvas *cV_CT_1 = new TCanvas("cV_CT_1", "cV_CT_1",w,h);
        cV_CT_1->SetGrid();
        gV_CT_1->Draw("AP");

        TCanvas *cV_CT_2 = new TCanvas("cV_CT_2", "cV_CT_2",w,h);
        cV_CT_2->SetGrid();
        gV_CT_2->Draw("AP");

        TCanvas *cV_CT_3 = new TCanvas("cV_CT_3", "cV_CT_3",w,h);
        cV_CT_3->SetGrid();
        gV_CT_3->Draw("AP");
    }


    //------------------------------

    auto legendDCR = new TLegend(0.15,0.70,0.35,0.85);
    legendDCR->AddEntry(gV_DCR_1,"HD3-2 (1)","p");
    legendDCR->AddEntry(gV_DCR_2,"HD3-2 (2)","p");
    legendDCR->AddEntry(gV_DCR_3,"HD3-2 (3)","p");

    auto legendCT = new TLegend(0.15,0.70,0.35,0.85);
    legendCT->AddEntry(gV_CT_1,"HD3-2 (1)","p");
    legendCT->AddEntry(gV_CT_2,"HD3-2 (2)","p");
    legendCT->AddEntry(gV_CT_3,"HD3-2 (3)","p");



    //------------------------------


    TCanvas *cDCR = new TCanvas("cDCR", "cDCR",w,h);
    cDCR->SetGrid();
    TMultiGraph *mgDCR = new TMultiGraph("mgDCR", title_DCR_mg);
    mgDCR->Add(gV_DCR_1);
    mgDCR->Add(gV_DCR_2);
    mgDCR->Add(gV_DCR_3);
    mgDCR->Draw("AP");
    legendDCR->Draw();

    //------------------------------

    TCanvas *cCT = new TCanvas("cCT", "cCT",w,h);
    cCT->SetGrid();
    TMultiGraph *mgCT = new TMultiGraph("mgCT", ";Bias Voltage (V); P_{CT}");
    mgCT->Add(gV_CT_1);
    mgCT->Add(gV_CT_2);
    mgCT->Add(gV_CT_3);
    mgCT->Draw("AP");
    legendCT->Draw();

    ///////////////////////////////////////////////////////////////////////////

    //------------------------------

    if(draw_all_bool){
        TCanvas *cV_DCR_Del_1 = new TCanvas("cV_DCR_Del_1", "cV_DCR_Del_1",w,h);
        cV_DCR_Del_1->SetGrid();
        gV_DCR_Del_1->Draw("AP");

        TCanvas *cV_DCR_Del_2 = new TCanvas("cV_DCR_Del_2", "cV_DCR_Del_2",w,h);
        cV_DCR_Del_2->SetGrid();
        gV_DCR_Del_2->Draw("AP");

        TCanvas *cV_DCR_Del_3 = new TCanvas("cV_DCR_Del_3", "cV_DCR_Del_3",w,h);
        cV_DCR_Del_3->SetGrid();
        gV_DCR_Del_3->Draw("AP");
    }


    //------------------------------

    if(draw_all_bool){
        TCanvas *cV_CT_Del_1 = new TCanvas("cV_CT_Del_1", "cV_CT_Del_1",w,h);
        cV_CT_Del_1->SetGrid();
        gV_CT_Del_1->Draw("AP");

        TCanvas *cV_CT_Del_2 = new TCanvas("cV_CT_Del_2", "cV_CT_Del_2",w,h);
        cV_CT_Del_2->SetGrid();
        gV_CT_Del_2->Draw("AP");

        TCanvas *cV_CT_Del_3 = new TCanvas("cV_CT_Del_3", "cV_CT_Del_3",w,h);
        cV_CT_Del_3->SetGrid();
        gV_CT_Del_3->Draw("AP");
    }


    //------------------------------


    auto legendDCR_Del = new TLegend(0.15,0.70,0.35,0.85);
    legendDCR_Del->AddEntry(gV_DCR_Del_1,"HD3-2 (1)","p");
    legendDCR_Del->AddEntry(gV_DCR_Del_2,"HD3-2 (2)","p");
    legendDCR_Del->AddEntry(gV_DCR_Del_3,"HD3-2 (3)","p");

    auto legendCT_Del = new TLegend(0.15,0.70,0.35,0.85);
    legendCT_Del->AddEntry(gV_CT_Del_1,"HD3-2 (1)","p");
    legendCT_Del->AddEntry(gV_CT_Del_2,"HD3-2 (2)","p");
    legendCT_Del->AddEntry(gV_CT_Del_3,"HD3-2 (3)","p");

    //------------------------------

    auto legendDCR_CNT_Del = new TLegend(0.15,0.70,0.35,0.85);
    legendDCR_CNT_Del->AddEntry(gV_DCR_1,"","p");
    legendDCR_CNT_Del->AddEntry(gV_DCR_Del_1,"HD3-2 (1)","p");
    legendDCR_CNT_Del->AddEntry(gV_DCR_2,"","p");
    legendDCR_CNT_Del->AddEntry(gV_DCR_Del_2,"HD3-2 (2)","p");
    legendDCR_CNT_Del->AddEntry(gV_DCR_3,"","p");
    legendDCR_CNT_Del->AddEntry(gV_DCR_Del_3,"HD3-2 (3)","p");
    legendDCR_CNT_Del->SetNColumns(2);

    auto legendCT_CNT_Del = new TLegend(0.15,0.70,0.35,0.85);
    legendCT_CNT_Del->AddEntry(gV_CT_1,"","p");
    legendCT_CNT_Del->AddEntry(gV_CT_Del_1,"HD3-2 (1)","p");
    legendCT_CNT_Del->AddEntry(gV_CT_2,"","p");
    legendCT_CNT_Del->AddEntry(gV_CT_Del_2,"HD3-2 (2)","p");
    legendCT_CNT_Del->AddEntry(gV_CT_3,"","p");
    legendCT_CNT_Del->AddEntry(gV_CT_Del_3,"HD3-2 (3)","p");
    legendCT_CNT_Del->SetNColumns(2);


    //------------------------------

    TCanvas *cDCR_Del = new TCanvas("cDCR_Del", "cDCR_Del",w,h);
    cDCR_Del->SetGrid();
    TMultiGraph *mgDCR_Del = new TMultiGraph("mgDCR_Del", title_DCR_mg);
    mgDCR_Del->Add(gV_DCR_Del_1);
    mgDCR_Del->Add(gV_DCR_Del_2);
    mgDCR_Del->Add(gV_DCR_Del_3);
    mgDCR_Del->Draw("AP");
    legendDCR_Del->Draw();


    //------------------------------

    TCanvas *cCT_Del = new TCanvas("cCT_Del", "cCT_Del",w,h);
    cCT_Del->SetGrid();
    TMultiGraph *mgCT_Del = new TMultiGraph("mgCT_Del", ";Bias Voltage (V); P_{CT}");
    mgCT_Del->Add(gV_CT_Del_1);
    mgCT_Del->Add(gV_CT_Del_2);
    mgCT_Del->Add(gV_CT_Del_3);
    mgCT_Del->Draw("AP");
    legendCT_Del->Draw();



    ///////////////////////////////////////////////////////////////////////////

    if(draw_all_bool){
        //------------------------------

        TCanvas *cDCR_CNT_Del_1 = new TCanvas("cDCR_CNT_Del_1", "cDCR_CNT_Del_1",w,h);
        cDCR_CNT_Del_1->SetGrid();
        TMultiGraph *mgDCR_CNT_Del_1 = new TMultiGraph("mgDCR_CNT_Del_1", title_DCR_mg);
        mgDCR_CNT_Del_1->Add(gV_DCR_1);
        mgDCR_CNT_Del_1->Add(gV_DCR_Del_1);
        mgDCR_CNT_Del_1->Draw("AP");

        //------------------------------

        TCanvas *cCT_CNT_Del_1 = new TCanvas("cCT_CNT_Del_1", "cCT_CNT_Del_1",w,h);
        cCT_CNT_Del_1->SetGrid();
        TMultiGraph *mgCT_CNT_Del_1 = new TMultiGraph("mgCT_CNT_Del_1", ";Bias Voltage (V); P_{CT}");
        mgCT_CNT_Del_1->Add(gV_CT_1);
        mgCT_CNT_Del_1->Add(gV_CT_Del_1);
        mgCT_CNT_Del_1->Draw("AP");


        //------------------------------
        //------------------------------

        TCanvas *cDCR_CNT_Del_2 = new TCanvas("cDCR_CNT_Del_2", "cDCR_CNT_Del_2",w,h);
        cDCR_CNT_Del_2->SetGrid();
        TMultiGraph *mgDCR_CNT_Del_2 = new TMultiGraph("mgDCR_CNT_Del_2", title_DCR_mg);
        mgDCR_CNT_Del_2->Add(gV_DCR_2);
        mgDCR_CNT_Del_2->Add(gV_DCR_Del_2);
        mgDCR_CNT_Del_2->Draw("AP");

        //------------------------------

        TCanvas *cCT_CNT_Del_2 = new TCanvas("cCT_CNT_Del_2", "cCT_CNT_Del_2",w,h);
        cCT_CNT_Del_2->SetGrid();
        TMultiGraph *mgCT_CNT_Del_2 = new TMultiGraph("mgCT_CNT_Del_2", ";Bias Voltage (V); P_{CT}");
        mgCT_CNT_Del_2->Add(gV_CT_2);
        mgCT_CNT_Del_2->Add(gV_CT_Del_2);
        mgCT_CNT_Del_2->Draw("AP");


        //------------------------------
        //------------------------------

        TCanvas *cDCR_CNT_Del_3 = new TCanvas("cDCR_CNT_Del_3", "cDCR_CNT_Del_3",w,h);
        cDCR_CNT_Del_3->SetGrid();
        TMultiGraph *mgDCR_CNT_Del_3 = new TMultiGraph("mgDCR_CNT_Del_3", title_DCR_mg);
        mgDCR_CNT_Del_3->Add(gV_DCR_3);
        mgDCR_CNT_Del_3->Add(gV_DCR_Del_3);
        mgDCR_CNT_Del_3->Draw("AP");

        //------------------------------

        TCanvas *cCT_CNT_Del_3 = new TCanvas("cCT_CNT_Del_3", "cCT_CNT_Del_3",w,h);
        cCT_CNT_Del_3->SetGrid();
        TMultiGraph *mgCT_CNT_Del_3 = new TMultiGraph("mgCT_CNT_Del_3", ";Bias Voltage (V); P_{CT}");
        mgCT_CNT_Del_3->Add(gV_CT_3);
        mgCT_CNT_Del_3->Add(gV_CT_Del_3);
        mgCT_CNT_Del_3->Draw("AP");


    }


    //------------------------------
    //------------------------------

    TCanvas *cDCR_CNT_Del = new TCanvas("cDCR_CNT_Del", "cDCR_CNT_Del",w,h);
    cDCR_CNT_Del->SetGrid();
    TMultiGraph *mgDCR_CNT_Del = new TMultiGraph("mgDCR_CNT_Del", title_DCR_mg);
    mgDCR_CNT_Del->Add(gV_DCR_1);
    mgDCR_CNT_Del->Add(gV_DCR_2);
    mgDCR_CNT_Del->Add(gV_DCR_3);
    mgDCR_CNT_Del->Add(gV_DCR_Del_1);
    mgDCR_CNT_Del->Add(gV_DCR_Del_2);
    mgDCR_CNT_Del->Add(gV_DCR_Del_3);
    mgDCR_CNT_Del->Draw("AP");
    legendDCR_CNT_Del->Draw();


    //------------------------------

    TCanvas *cCT_CNT_Del = new TCanvas("cCT_CNT_Del", "cCT_CNT_Del",w,h);
    cCT_CNT_Del->SetGrid();
    TMultiGraph *mgCT_CNT_Del = new TMultiGraph("mgCT_CNT_Del", ";Bias Voltage (V); P_{CT}");
    mgCT_CNT_Del->Add(gV_CT_1);
    mgCT_CNT_Del->Add(gV_CT_2);
    mgCT_CNT_Del->Add(gV_CT_3);
    mgCT_CNT_Del->Add(gV_CT_Del_1);
    mgCT_CNT_Del->Add(gV_CT_Del_2);
    mgCT_CNT_Del->Add(gV_CT_Del_3);
    mgCT_CNT_Del->Draw("AP");
    legendCT_CNT_Del->Draw();

    cout<<"###############################################################################"<<endl;
    cout<<" WARNING "<<endl;
    if(dcr_area) cout<<"DCR / AREA, please check axis title"<<endl;
    else         cout<<" DCR global, not / area"<<endl;
    cout<<"###############################################################################"<<endl;




    //------------------------------

    if(dcr_area){ // /= Area
        for(int i=0; i< n_DCR; i++){
            DCR_1_Ea[i]  /= Area;
            errDCR_1_Ea[i] /= Area;
            DCR_1_Ea[i]  *= 1e3;
            errDCR_1_Ea[i] *= 1e3;
        }
        for(int i=0; i< n_DCR; i++){
            DCR_2_Ea[i]  /= Area;
            errDCR_2_Ea[i] /= Area;
            DCR_2_Ea[i]  *= 1e3;
            errDCR_2_Ea[i] *= 1e3;
        }
        for(int i=0; i< n_DCR; i++){
            DCR_3_Ea[i]  /= Area;
            errDCR_3_Ea[i] /= Area;
            DCR_3_Ea[i]  *= 1e3;
            errDCR_3_Ea[i] *= 1e3;
        }

        //------------------------------

        for(int i=0; i< n_DCR; i++){
            DCR_Del_1_Ea[i]  /= Area;
            errDCR_Del_1_Ea[i] /= Area;
            DCR_Del_1_Ea[i]  *= 1e3;
            errDCR_Del_1_Ea[i] *= 1e3;
        }
        for(int i=0; i< n_DCR; i++){
            DCR_Del_2_Ea[i]  /= Area;
            errDCR_Del_2_Ea[i] /= Area;
            DCR_Del_2_Ea[i]  *= 1e3;
            errDCR_Del_2_Ea[i] *= 1e3;
        }
        for(int i=0; i< n_DCR; i++){
            DCR_Del_3_Ea[i]  /= Area;
            errDCR_Del_3_Ea[i] /= Area;
            DCR_Del_3_Ea[i]  *= 1e3;
            errDCR_Del_3_Ea[i] *= 1e3;
        }


        for(int i=0; i< n_DCR; i++){
            DCR_1_Eb[i]  /= Area;
            errDCR_1_Eb[i] /= Area;
            DCR_1_Eb[i]  *= 1e3;
            errDCR_1_Eb[i] *= 1e3;
        }
        for(int i=0; i< n_DCR; i++){
            DCR_2_Eb[i]  /= Area;
            errDCR_2_Eb[i] /= Area;
            DCR_2_Eb[i]  *= 1e3;
            errDCR_2_Eb[i] *= 1e3;
        }
        for(int i=0; i< n_DCR; i++){
            DCR_3_Eb[i]  /= Area;
            errDCR_3_Eb[i] /= Area;
            DCR_3_Eb[i]  *= 1e3;
            errDCR_3_Eb[i] *= 1e3;
        }

        //------------------------------

        for(int i=0; i< n_DCR; i++){
            DCR_Del_1_Eb[i]  /= Area;
            errDCR_Del_1_Eb[i] /= Area;
            DCR_Del_1_Eb[i]  *= 1e3;
            errDCR_Del_1_Eb[i] *= 1e3;
        }
        for(int i=0; i< n_DCR; i++){
            DCR_Del_2_Eb[i]  /= Area;
            errDCR_Del_2_Eb[i] /= Area;
            DCR_Del_2_Eb[i]  *= 1e3;
            errDCR_Del_2_Eb[i] *= 1e3;
        }
        for(int i=0; i< n_DCR; i++){
            DCR_Del_3_Eb[i]  /= Area;
            errDCR_Del_3_Eb[i] /= Area;
            DCR_Del_3_Eb[i]  *= 1e3;
            errDCR_Del_3_Eb[i] *= 1e3;
        }

    } // end /= Area


    TGraphErrors *gV_DCR_1_Ea  = new TGraphErrors(n_DCR, HV_1_Ea, DCR_1_Ea, errHV_1_Ea, errDCR_1_Ea);
    TGraphErrors *gV_DCR_2_Ea  = new TGraphErrors(n_DCR, HV_2_Ea, DCR_2_Ea, errHV_2_Ea, errDCR_2_Ea);
    TGraphErrors *gV_DCR_3_Ea  = new TGraphErrors(n_DCR, HV_3_Ea, DCR_3_Ea, errHV_3_Ea, errDCR_3_Ea);
    TGraphErrors *gV_CT_1_Ea  = new TGraphErrors(n_DCR, HV_1_Ea, CT_1_Ea, errHV_1_Ea, errCT_1_Ea);
    TGraphErrors *gV_CT_2_Ea  = new TGraphErrors(n_DCR, HV_2_Ea, CT_2_Ea, errHV_2_Ea, errCT_2_Ea);
    TGraphErrors *gV_CT_3_Ea  = new TGraphErrors(n_DCR, HV_3_Ea, CT_3_Ea, errHV_3_Ea, errCT_3_Ea);
    TGraphErrors *gV_DCR_Del_1_Ea  = new TGraphErrors(n_DCR, HV_1_Ea, DCR_Del_1_Ea, errHV_1_Ea, errDCR_Del_1_Ea);
    TGraphErrors *gV_DCR_Del_2_Ea  = new TGraphErrors(n_DCR, HV_2_Ea, DCR_Del_2_Ea, errHV_2_Ea, errDCR_Del_2_Ea);
    TGraphErrors *gV_DCR_Del_3_Ea  = new TGraphErrors(n_DCR, HV_3_Ea, DCR_Del_3_Ea, errHV_3_Ea, errDCR_3_Ea);
    TGraphErrors *gV_CT_Del_1_Ea  = new TGraphErrors(n_DCR, HV_1_Ea, CT_Del_1_Ea, errHV_1_Ea, errCT_Del_1_Ea);
    TGraphErrors *gV_CT_Del_2_Ea  = new TGraphErrors(n_DCR, HV_2_Ea, CT_Del_2_Ea, errHV_2_Ea, errCT_Del_2_Ea);
    TGraphErrors *gV_CT_Del_3_Ea  = new TGraphErrors(n_DCR, HV_3_Ea, CT_Del_3_Ea, errHV_3_Ea, errCT_Del_3_Ea);

    TGraphErrors *gV_DCR_1_Eb  = new TGraphErrors(n_DCR, HV_1_Eb, DCR_1_Eb, errHV_1_Eb, errDCR_1_Eb);
    TGraphErrors *gV_DCR_2_Eb  = new TGraphErrors(n_DCR, HV_2_Eb, DCR_2_Eb, errHV_2_Eb, errDCR_2_Eb);
    TGraphErrors *gV_DCR_3_Eb  = new TGraphErrors(n_DCR, HV_3_Eb, DCR_3_Eb, errHV_3_Eb, errDCR_3_Eb);
    TGraphErrors *gV_CT_1_Eb  = new TGraphErrors(n_DCR, HV_1_Eb, CT_1_Eb, errHV_1_Eb, errCT_1_Eb);
    TGraphErrors *gV_CT_2_Eb  = new TGraphErrors(n_DCR, HV_2_Eb, CT_2_Eb, errHV_2_Eb, errCT_2_Eb);
    TGraphErrors *gV_CT_3_Eb  = new TGraphErrors(n_DCR, HV_3_Eb, CT_3_Eb, errHV_3_Eb, errCT_3_Eb);
    TGraphErrors *gV_DCR_Del_1_Eb  = new TGraphErrors(n_DCR, HV_1_Eb, DCR_Del_1_Eb, errHV_1_Eb, errDCR_Del_1_Eb);
    TGraphErrors *gV_DCR_Del_2_Eb  = new TGraphErrors(n_DCR, HV_2_Eb, DCR_Del_2_Eb, errHV_2_Eb, errDCR_Del_2_Eb);
    TGraphErrors *gV_DCR_Del_3_Eb  = new TGraphErrors(n_DCR, HV_3_Eb, DCR_Del_3_Eb, errHV_3_Eb, errDCR_3_Eb);
    TGraphErrors *gV_CT_Del_1_Eb  = new TGraphErrors(n_DCR, HV_1_Eb, CT_Del_1_Eb, errHV_1_Eb, errCT_Del_1_Eb);
    TGraphErrors *gV_CT_Del_2_Eb  = new TGraphErrors(n_DCR, HV_2_Eb, CT_Del_2_Eb, errHV_2_Eb, errCT_Del_2_Eb);
    TGraphErrors *gV_CT_Del_3_Eb  = new TGraphErrors(n_DCR, HV_3_Eb, CT_Del_3_Eb, errHV_3_Eb, errCT_Del_3_Eb);

    //----------

    gV_DCR_1_Ea->SetMarkerStyle(4);
    gV_DCR_2_Ea->SetMarkerStyle(4);
    gV_DCR_3_Ea->SetMarkerStyle(4);
    gV_CT_1_Ea->SetMarkerStyle(4);
    gV_CT_2_Ea->SetMarkerStyle(4);
    gV_CT_3_Ea->SetMarkerStyle(4);
    gV_DCR_Del_1_Ea->SetMarkerStyle(4);
    gV_DCR_Del_2_Ea->SetMarkerStyle(4);
    gV_DCR_Del_3_Ea->SetMarkerStyle(4);
    gV_CT_Del_1_Ea->SetMarkerStyle(4);
    gV_CT_Del_2_Ea->SetMarkerStyle(4);
    gV_CT_Del_3_Ea->SetMarkerStyle(4);
    gV_DCR_1_Eb->SetMarkerStyle(4);
    gV_DCR_2_Eb->SetMarkerStyle(4);
    gV_DCR_3_Eb->SetMarkerStyle(4);
    gV_CT_1_Eb->SetMarkerStyle(4);
    gV_CT_2_Eb->SetMarkerStyle(4);
    gV_CT_3_Eb->SetMarkerStyle(4);
    gV_DCR_Del_1_Eb->SetMarkerStyle(4);
    gV_DCR_Del_2_Eb->SetMarkerStyle(4);
    gV_DCR_Del_3_Eb->SetMarkerStyle(4);
    gV_CT_Del_1_Eb->SetMarkerStyle(4);
    gV_CT_Del_2_Eb->SetMarkerStyle(4);
    gV_CT_Del_3_Eb->SetMarkerStyle(4);


    bool draw_check_bool = false;

    if(draw_check_bool){
        TCanvas *cV_DCR_1_E = new TCanvas("cV_DCR_1_E", "cV_DCR_1_E",w,h);
        cV_DCR_1_E->SetGrid();
        TMultiGraph *mgV_DCR_1_E = new TMultiGraph("mgV_DCR_1_E", ";Bias Voltage (V); DCR");
        mgV_DCR_1_E->Add(gV_DCR_1);
        mgV_DCR_1_E->Add(gV_DCR_1_Ea);
        mgV_DCR_1_E->Add(gV_DCR_1_Eb);
        mgV_DCR_1_E->Draw("AP");

        TCanvas *cV_DCR_2_E = new TCanvas("cV_DCR_2_E", "cV_DCR_2_E",w,h);
        cV_DCR_2_E->SetGrid();
        TMultiGraph *mgV_DCR_2_E = new TMultiGraph("mgV_DCR_2_E", ";Bias Voltage (V); DCR");
        mgV_DCR_2_E->Add(gV_DCR_2);
        mgV_DCR_2_E->Add(gV_DCR_2_Ea);
        mgV_DCR_2_E->Add(gV_DCR_2_Eb);
        mgV_DCR_2_E->Draw("AP");

        TCanvas *cV_DCR_3_E = new TCanvas("cV_DCR_3_E", "cV_DCR_3_E",w,h);
        cV_DCR_3_E->SetGrid();
        TMultiGraph *mgV_DCR_3_E = new TMultiGraph("mgV_DCR_3_E", ";Bias Voltage (V); DCR");
        mgV_DCR_3_E->Add(gV_DCR_3);
        mgV_DCR_3_E->Add(gV_DCR_3_Ea);
        mgV_DCR_3_E->Add(gV_DCR_3_Eb);
        mgV_DCR_3_E->Draw("AP");

        //----------

        TCanvas *cV_CT_1_E = new TCanvas("cV_CT_1_E", "cV_CT_1_E",w,h);
        cV_CT_1_E->SetGrid();
        TMultiGraph *mgV_CT_1_E = new TMultiGraph("mgV_CT_1_E", ";Bias Voltage (V); CT");
        mgV_CT_1_E->Add(gV_CT_1);
        mgV_CT_1_E->Add(gV_CT_1_Ea);
        mgV_CT_1_E->Add(gV_CT_1_Eb);
        mgV_CT_1_E->Draw("AP");

        TCanvas *cV_CT_2_E = new TCanvas("cV_CT_2_E", "cV_CT_2_E",w,h);
        cV_CT_2_E->SetGrid();
        TMultiGraph *mgV_CT_2_E = new TMultiGraph("mgV_CT_2_E", ";Bias Voltage (V); CT");
        mgV_CT_2_E->Add(gV_CT_2);
        mgV_CT_2_E->Add(gV_CT_2_Ea);
        mgV_CT_2_E->Add(gV_CT_2_Eb);
        mgV_CT_2_E->Draw("AP");

        TCanvas *cV_CT_3_E = new TCanvas("cV_CT_3_E", "cV_CT_3_E",w,h);
        cV_CT_3_E->SetGrid();
        TMultiGraph *mgV_CT_3_E = new TMultiGraph("mgV_CT_3_E", ";Bias Voltage (V); CT");
        mgV_CT_3_E->Add(gV_CT_3);
        mgV_CT_3_E->Add(gV_CT_3_Ea);
        mgV_CT_3_E->Add(gV_CT_3_Eb);
        mgV_CT_3_E->Draw("AP");

        //----------

        TCanvas *cV_DCR_Del_1_E = new TCanvas("cV_DCR_Del_1_E", "cV_DCR_Del_1_E",w,h);
        cV_DCR_Del_1_E->SetGrid();
        TMultiGraph *mgV_DCR_Del_1_E = new TMultiGraph("mgV_DCR_Del_1_E", ";Bias Voltage (V); DCR_Del");
        mgV_DCR_Del_1_E->Add(gV_DCR_Del_1);
        mgV_DCR_Del_1_E->Add(gV_DCR_Del_1_Ea);
        mgV_DCR_Del_1_E->Add(gV_DCR_Del_1_Eb);
        mgV_DCR_Del_1_E->Draw("AP");

        TCanvas *cV_DCR_Del_2_E = new TCanvas("cV_DCR_Del_2_E", "cV_DCR_Del_2_E",w,h);
        cV_DCR_Del_2_E->SetGrid();
        TMultiGraph *mgV_DCR_Del_2_E = new TMultiGraph("mgV_DCR_Del_2_E", ";Bias Voltage (V); DCR_Del");
        mgV_DCR_Del_2_E->Add(gV_DCR_Del_2);
        mgV_DCR_Del_2_E->Add(gV_DCR_Del_2_Ea);
        mgV_DCR_Del_2_E->Add(gV_DCR_Del_2_Eb);
        mgV_DCR_Del_2_E->Draw("AP");

        TCanvas *cV_DCR_Del_3_E = new TCanvas("cV_DCR_Del_3_E", "cV_DCR_Del_3_E",w,h);
        cV_DCR_Del_3_E->SetGrid();
        TMultiGraph *mgV_DCR_Del_3_E = new TMultiGraph("mgV_DCR_Del_3_E", ";Bias Voltage (V); DCR_Del");
        mgV_DCR_Del_3_E->Add(gV_DCR_Del_3);
        mgV_DCR_Del_3_E->Add(gV_DCR_Del_3_Ea);
        mgV_DCR_Del_3_E->Add(gV_DCR_Del_3_Eb);
        mgV_DCR_Del_3_E->Draw("AP");

        //----------

        TCanvas *cV_CT_Del_1_E = new TCanvas("cV_CT_Del_1_E", "cV_CT_Del_1_E",w,h);
        cV_CT_Del_1_E->SetGrid();
        TMultiGraph *mgV_CT_Del_1_E = new TMultiGraph("mgV_CT_Del_1_E", ";Bias Voltage (V); CT_Del");
        mgV_CT_Del_1_E->Add(gV_CT_Del_1);
        mgV_CT_Del_1_E->Add(gV_CT_Del_1_Ea);
        mgV_CT_Del_1_E->Add(gV_CT_Del_1_Eb);
        mgV_CT_Del_1_E->Draw("AP");

        TCanvas *cV_CT_Del_2_E = new TCanvas("cV_CT_Del_2_E", "cV_CT_Del_2_E",w,h);
        cV_CT_Del_2_E->SetGrid();
        TMultiGraph *mgV_CT_Del_2_E = new TMultiGraph("mgV_CT_Del_2_E", ";Bias Voltage (V); CT_Del");
        mgV_CT_Del_2_E->Add(gV_CT_Del_2);
        mgV_CT_Del_2_E->Add(gV_CT_Del_2_Ea);
        mgV_CT_Del_2_E->Add(gV_CT_Del_2_Eb);
        mgV_CT_Del_2_E->Draw("AP");

        TCanvas *cV_CT_Del_3_E = new TCanvas("cV_CT_Del_3_E", "cV_CT_Del_3_E",w,h);
        cV_CT_Del_3_E->SetGrid();
        TMultiGraph *mgV_CT_Del_3_E = new TMultiGraph("mgV_CT_Del_3_E", ";Bias Voltage (V); CT_Del");
        mgV_CT_Del_3_E->Add(gV_CT_Del_3);
        mgV_CT_Del_3_E->Add(gV_CT_Del_3_Ea);
        mgV_CT_Del_3_E->Add(gV_CT_Del_3_Eb);
        mgV_CT_Del_3_E->Draw("AP");

    }

    //------------------------------
    // For LaTeX
    //------------------------------
    // cout<<endl<<endl;
    // cout<<"For LaTeX"<<endl<<endl;
    // cout<<"\% From file DCR_CrossTalk_FBK_HD3_2_from_cnt_data_2018_07.C"<<endl;
    // cout<<"\% From cnt"<<endl;
    // cout<<"\% HV & DCR_1 (kHz/mm^2) & DCR_2 (kHz/mm^2) & DCR_3 (kHz/mm^2) & CT_1  & CT_2  & CT_3  \\\\"<<endl;
    // for(int i=0; i<n_DCR; i++){
    //     printf("$ %.2f \\pm %.2f $ & $ %.2f \\pm %.2f $ & $ %.2f \\pm %.2f $ & $ %.2f \\pm %.2f $ & $ %.2f \\pm %.2f $ & $ %.2f \\pm %.2f $ & $ %.2f \\pm %.2f $ \\\\ \n", HV_1[i], errHV_1[i], DCR_1[i], errDCR_1[i],DCR_2[i], errDCR_2[i],DCR_3[i], errDCR_3[i], CT_1[i], errCT_1[i],CT_2[i], errCT_2[i],CT_3[i], errCT_3[i]);
    // }
    cout<<endl<<endl;
    cout<<"-------------------------------------------------"<<endl;
    cout<<endl<<endl;
    cout<<"For LaTeX"<<endl<<endl;
    cout<<"\% From file DCR_CrossTalk_FBK_HD3_2_from_cnt_data_2018_07.C"<<endl;
    cout<<"\% From cnt"<<endl;
    cout<<"\% HV & DCR_1 (kHz/mm^2) & DCR_2 (kHz/mm^2) & DCR_3 (kHz/mm^2) \\\\"<<endl;
    for(int i=0; i<n_DCR; i++){
        printf("$ %.2f \\pm %.2f $ & $ %.2f \\pm %.2f $ & $ %.2f \\pm %.2f $ & $ %.2f \\pm %.2f $  \\\\ \n", HV_1[i], errHV_1[i], DCR_1[i], errDCR_1[i],DCR_2[i], errDCR_2[i],DCR_3[i], errDCR_3[i]);
    }
    cout<<endl<<endl;
    cout<<"\% From file DCR_CrossTalk_FBK_HD3_2_from_cnt_data_2018_07.C"<<endl;
    cout<<"\% From cnt"<<endl;
    cout<<"\% HV & CT_1 & CT_2  & CT_3  \\\\"<<endl;
    for(int i=0; i<n_DCR; i++){
        printf("$ %.2f \\pm %.2f $ & $ %.2f \\pm %.2f $ & $ %.2f \\pm %.2f $ & $ %.2f \\pm %.2f $  \\\\ \n", HV_1[i], errHV_1[i], CT_1[i], errCT_1[i],CT_2[i], errCT_2[i],CT_3[i], errCT_3[i]);
    }

    cout<<endl<<endl;
    cout<<"-------------------------------------------------"<<endl;
    cout<<endl<<endl;
    cout<<"\% From file DCR_CrossTalk_FBK_HD3_2_from_cnt_data_2018_07.C"<<endl;
    cout<<"\% From Del"<<endl;
    cout<<"\% HV & DCR_Del_1 (kHz/mm^2) & DCR_Del_2 (kHz/mm^2) & DCR_Del_3 (kHz/mm^2) \\\\"<<endl;
    for(int i=0; i<n_DCR; i++){
        printf("$ %.2f \\pm %.2f $ & $ %.2f \\pm %.2f $ & $ %.2f \\pm %.2f $ & $ %.2f \\pm %.2f $  \\\\ \n", HV_1[i], errHV_1[i], DCR_Del_1[i], errDCR_Del_1[i],DCR_Del_2[i], errDCR_Del_2[i],DCR_Del_3[i], errDCR_Del_3[i]);
    }
    cout<<endl<<endl;
    cout<<"\% From file DCR_CrossTalk_FBK_HD3_2_from_cnt_data_2018_07.C"<<endl;
    cout<<"\% From Del"<<endl;
    cout<<"\% HV & CT_Del_1 & CT_Del_2  & CT_Del_3  \\\\"<<endl;
    for(int i=0; i<n_DCR; i++){
        printf("$ %.2f \\pm %.2f $ & $ %.2f \\pm %.2f $ & $ %.2f \\pm %.2f $ & $ %.2f \\pm %.2f $  \\\\ \n", HV_1[i], errHV_1[i], CT_Del_1[i], errCT_Del_1[i],CT_Del_2[i], errCT_Del_2[i],CT_Del_3[i], errCT_Del_3[i]);
    }





    ////////////////////////////////////////////////////////////////////////////
    //  PRINT
    ////////////////////////////////////////////////////////////////////////////
    cout<<endl<<endl;
    cout<<"//---------------"<<endl;
    cout<<"// FBK HD3-2 DARK"<<endl;
    cout<<"//---------------"<<endl;
    cout<<"double HV_FBK_DARK[] = {"; for(int i=0; i<n_DCR-1; i++) cout<<HV_1[i]<<", "; cout<<HV_1[n_DCR-1]<<"};"<<endl;
    cout<<"double errHV_FBK_DARK[] = {"; for(int i=0; i<n_DCR-1; i++) cout<<errHV_1[i]<<", "; cout<<errHV_1[n_DCR-1]<<"};"<<endl;
    cout<<endl;
    cout<<"// DCR"<<endl;
    cout<<"double DCR_1_FBK_DARK[] = {"; for(int i=0; i<n_DCR-1; i++) cout<<DCR_1[i]<<", "; cout<<DCR_1[n_DCR-1]<<"};"<<endl;
    cout<<"double errDCR_1_FBK_DARK[] = {"; for(int i=0; i<n_DCR-1; i++) cout<<errDCR_1[i]<<", "; cout<<errDCR_1[n_DCR-1]<<"};"<<endl;
    cout<<"double DCR_Del_1_FBK_DARK[] = {"; for(int i=0; i<n_DCR-1; i++) cout<<DCR_Del_1[i]<<", "; cout<<DCR_Del_1[n_DCR-1]<<"};"<<endl;
    cout<<"double errDCR_Del_1_FBK_DARK[] = {"; for(int i=0; i<n_DCR-1; i++) cout<<errDCR_Del_1[i]<<", "; cout<<errDCR_Del_1[n_DCR-1]<<"};"<<endl;

    cout<<"double DCR_2_FBK_DARK[] = {"; for(int i=0; i<n_DCR-1; i++) cout<<DCR_2[i]<<", "; cout<<DCR_2[n_DCR-1]<<"};"<<endl;
    cout<<"double errDCR_2_FBK_DARK[] = {"; for(int i=0; i<n_DCR-1; i++) cout<<errDCR_2[i]<<", "; cout<<errDCR_2[n_DCR-1]<<"};"<<endl;
    cout<<"double DCR_Del_2_FBK_DARK[] = {"; for(int i=0; i<n_DCR-1; i++) cout<<DCR_Del_2[i]<<", "; cout<<DCR_Del_2[n_DCR-1]<<"};"<<endl;
    cout<<"double errDCR_Del_2_FBK_DARK[] = {"; for(int i=0; i<n_DCR-1; i++) cout<<errDCR_Del_2[i]<<", "; cout<<errDCR_Del_2[n_DCR-1]<<"};"<<endl;

    cout<<"double DCR_3_FBK_DARK[] = {"; for(int i=0; i<n_DCR-1; i++) cout<<DCR_3[i]<<", "; cout<<DCR_3[n_DCR-1]<<"};"<<endl;
    cout<<"double errDCR_3_FBK_DARK[] = {"; for(int i=0; i<n_DCR-1; i++) cout<<errDCR_3[i]<<", "; cout<<errDCR_3[n_DCR-1]<<"};"<<endl;
    cout<<"double DCR_Del_3_FBK_DARK[] = {"; for(int i=0; i<n_DCR-1; i++) cout<<DCR_Del_3[i]<<", "; cout<<DCR_Del_3[n_DCR-1]<<"};"<<endl;
    cout<<"double errDCR_Del_3_FBK_DARK[] = {"; for(int i=0; i<n_DCR-1; i++) cout<<errDCR_Del_3[i]<<", "; cout<<errDCR_Del_3[n_DCR-1]<<"};"<<endl;

    cout<<endl;
    cout<<"// CT"<<endl;
    cout<<"double CT_1_FBK_DARK[] = {"; for(int i=0; i<n_DCR-1; i++) cout<<CT_1[i]<<", "; cout<<CT_1[n_DCR-1]<<"};"<<endl;
    cout<<"double errCT_1_FBK_DARK[] = {"; for(int i=0; i<n_DCR-1; i++) cout<<errCT_1[i]<<", "; cout<<errCT_1[n_DCR-1]<<"};"<<endl;
    cout<<"double CT_Del_1_FBK_DARK[] = {"; for(int i=0; i<n_DCR-1; i++) cout<<CT_Del_1[i]<<", "; cout<<CT_Del_1[n_DCR-1]<<"};"<<endl;
    cout<<"double errCT_Del_1_FBK_DARK[] = {"; for(int i=0; i<n_DCR-1; i++) cout<<errCT_Del_1[i]<<", "; cout<<errCT_Del_1[n_DCR-1]<<"};"<<endl;

    cout<<"double CT_2_FBK_DARK[] = {"; for(int i=0; i<n_DCR-1; i++) cout<<CT_2[i]<<", "; cout<<CT_2[n_DCR-1]<<"};"<<endl;
    cout<<"double errCT_2_FBK_DARK[] = {"; for(int i=0; i<n_DCR-1; i++) cout<<errCT_2[i]<<", "; cout<<errCT_2[n_DCR-1]<<"};"<<endl;
    cout<<"double CT_Del_2_FBK_DARK[] = {"; for(int i=0; i<n_DCR-1; i++) cout<<CT_Del_2[i]<<", "; cout<<CT_Del_2[n_DCR-1]<<"};"<<endl;
    cout<<"double errCT_Del_2_FBK_DARK[] = {"; for(int i=0; i<n_DCR-1; i++) cout<<errCT_Del_2[i]<<", "; cout<<errCT_Del_2[n_DCR-1]<<"};"<<endl;

    cout<<"double CT_3_FBK_DARK[] = {"; for(int i=0; i<n_DCR-1; i++) cout<<CT_3[i]<<", "; cout<<CT_3[n_DCR-1]<<"};"<<endl;
    cout<<"double errCT_3_FBK_DARK[] = {"; for(int i=0; i<n_DCR-1; i++) cout<<errCT_3[i]<<", "; cout<<errCT_3[n_DCR-1]<<"};"<<endl;
    cout<<"double CT_Del_3_FBK_DARK[] = {"; for(int i=0; i<n_DCR-1; i++) cout<<CT_Del_3[i]<<", "; cout<<CT_Del_3[n_DCR-1]<<"};"<<endl;
    cout<<"double errCT_Del_3_FBK_DARK[] = {"; for(int i=0; i<n_DCR-1; i++) cout<<errCT_Del_3[i]<<", "; cout<<errCT_Del_3[n_DCR-1]<<"};"<<endl;



    ////////////////////////////////////////////////////////////////////////////
    //          OLD
    ////////////////////////////////////////////////////////////////////////////

    // // FIXED ERROR
    // if(fix_error_bool){ // fix_error_bool
    //     double err_fix_DCR = 1;
    //     double err_fix_CT = 0.02;
    //     for(int i=0; i<n_DCR; i++){
    //         errDCR_1[i] = err_fix_DCR;
    //     }
    //     for(int i=0; i<n_DCR; i++){
    //         errDCR_2[i] = err_fix_DCR;
    //     }
    //     for(int i=0; i<n_DCR; i++){
    //         errDCR_3[i] = err_fix_DCR;
    //     }
    //     for(int i=0; i<n_DCR; i++){
    //         errCT_1[i] = err_fix_CT;
    //     }
    //     for(int i=0; i<n_DCR; i++){
    //         errCT_2[i] = err_fix_CT;
    //     }
    //     for(int i=0; i<n_DCR; i++){
    //         errCT_3[i] = err_fix_CT;
    //     }
    //     for(int i=0; i<n_DCR; i++){
    //         errDCR_Del_1[i] = err_fix_DCR;
    //     }
    //     for(int i=0; i<n_DCR; i++){
    //         errDCR_Del_2[i] = err_fix_DCR;
    //     }
    //     for(int i=0; i<n_DCR; i++){
    //         errDCR_Del_3[i] = err_fix_DCR;
    //     }
    //     for(int i=0; i<n_DCR; i++){
    //         errCT_Del_1[i] = err_fix_CT;
    //     }
    //     for(int i=0; i<n_DCR; i++){
    //         errCT_Del_2[i] = err_fix_CT;
    //     }
    //     for(int i=0; i<n_DCR; i++){
    //         errCT_Del_3[i] = err_fix_CT;
    //     }
    // }

    // // PERCENTAGE ERROR
    // if(percentage_error_bool){
    //     double err_rel = 0.05;
    //     for(int i=0; i<n_DCR; i++){
    //         errDCR_1[i] = err_rel * DCR_1[i];
    //     }
    //     for(int i=0; i<n_DCR; i++){
    //         errDCR_2[i] = err_rel * DCR_2[i];
    //     }
    //     for(int i=0; i<n_DCR; i++){
    //         errDCR_3[i] = err_rel * DCR_3[i];
    //     }
    //     for(int i=0; i<n_DCR; i++){
    //         errCT_1[i] = err_rel * CT_1[i];
    //     }
    //     for(int i=0; i<n_DCR; i++){
    //         errCT_2[i] = err_rel * CT_2[i];
    //     }
    //     for(int i=0; i<n_DCR; i++){
    //         errCT_3[i] = err_rel * CT_3[i];
    //     }
    //     for(int i=0; i<n_DCR; i++){
    //         errDCR_Del_1[i] = err_rel * DCR_Del_1[i];
    //     }
    //     for(int i=0; i<n_DCR; i++){
    //         errDCR_Del_2[i] = err_rel * DCR_Del_2[i];
    //     }
    //     for(int i=0; i<n_DCR; i++){
    //         errDCR_Del_3[i] = err_rel * DCR_Del_3[i];
    //     }
    //     for(int i=0; i<n_DCR; i++){
    //         errCT_Del_1[i] = err_rel * CT_Del_1[i];
    //     }
    //     for(int i=0; i<n_DCR; i++){
    //         errCT_Del_2[i] = err_rel * CT_Del_2[i];
    //     }
    //     for(int i=0; i<n_DCR; i++){
    //         errCT_Del_3[i] = err_rel * CT_Del_3[i];
    //     }
    // }

    // TEST Z
    // int N_Z = n_DCR;
    // double testZ_DCR[3][N_Z];
    // double testZ_CT[3][N_Z];
    // int n_SiPM = 0;
    // int n_SiPM_tot = 3;
    //
    // n_SiPM = 0;
    // for(int i=0; i<N_Z; i++){
    //     testZ_DCR[n_SiPM][i] = TMath::Abs((DCR_1[i] - DCR_Del_1[i]) / (TMath::Sqrt( errDCR_1[i]*errDCR_1[i] + errDCR_Del_1[i]*errDCR_Del_1[i] )));
    // }
    //
    // n_SiPM = 1;
    // for(int i=0; i<N_Z; i++){
    //     testZ_DCR[n_SiPM][i] = TMath::Abs((DCR_2[i] - DCR_Del_2[i]) / (TMath::Sqrt( errDCR_2[i]*errDCR_2[i] + errDCR_Del_2[i]*errDCR_Del_2[i] )));
    // }
    //
    // n_SiPM = 2;
    // for(int i=0; i<N_Z; i++){
    //     testZ_DCR[n_SiPM][i] = TMath::Abs((DCR_2[i] - DCR_Del_2[i]) / (TMath::Sqrt( errDCR_2[i]*errDCR_2[i] + errDCR_Del_2[i]*errDCR_Del_2[i] )));
    // }
    //
    // n_SiPM = 0;
    // for(int i=0; i<N_Z; i++){
    //     testZ_CT[n_SiPM][i] = TMath::Abs((CT_1[i] - CT_Del_1[i]) / (TMath::Sqrt( errCT_1[i]*errCT_1[i] + errCT_Del_1[i]*errCT_Del_1[i] )));
    // }
    //
    // n_SiPM = 1;
    // for(int i=0; i<N_Z; i++){
    //     testZ_CT[n_SiPM][i] = TMath::Abs((CT_2[i] - CT_Del_2[i]) / (TMath::Sqrt( errCT_2[i]*errCT_2[i] + errCT_Del_2[i]*errCT_Del_2[i] )));
    // }
    //
    // n_SiPM = 2;
    // for(int i=0; i<N_Z; i++){
    //     testZ_CT[n_SiPM][i] = TMath::Abs((CT_2[i] - CT_Del_2[i]) / (TMath::Sqrt( errCT_2[i]*errCT_2[i] + errCT_Del_2[i]*errCT_Del_2[i] )));
    // }

    // cout<<endl;
    // cout<<"TEST Z"<<endl;
    // for(int n=0; n<n_SiPM_tot; n++){
    //     cout<<"SiPM "<<n+1<<endl;
    //     for(int i=0; i<N_Z; i++){
    //         cout<<testZ_DCR[n][i];
    //         if(testZ_DCR[n][i] > 1.96) cout<<"\tTEST NOT PASSED"<<endl;
    //         else cout<<endl;
    //
    //         cout<<testZ_CT[n][i];
    //         if(testZ_CT[n][i] > 1.96) cout<<"\tTEST NOT PASSED "<<endl;
    //         else cout<<endl;
    //     }
    //     cout<<endl;
    // }




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


double GetMaxOrPercentage(double num, double LimUp, double LimDown, double min_percentage, double max_percentage){
    double err, err_rel;

    err = max(TMath::Abs(num-LimUp), TMath::Abs(num-LimDown));
    err_rel = err/num;

    if(err_rel > max_percentage){
        // err = num*max_percentage;
        cout<<"ERROR > max_percentage; err_rel = "<<err_rel<<endl;
    }else{
        if(err_rel < min_percentage){
            // err = num*min_percentage;
            cout<<"ERROR < min_percentage; err_rel = "<<err_rel<<endl;
        }
    }

    return err;

}
