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

    double min_percentage = 0.02;
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





    ///////////////////////////////////////////////////////////////////////////
    //      SiPM1
    ///////////////////////////////////////////////////////////////////////////
    // pe_0_5_vect[7] = {8, 8, 8, 9, 9, 9, 10}
    // pe_1_5_vect[7] = {19, 23, 27, 30, 35, 37, 40}
    //minyhistDelays = 20;  maxyhistDelays = 150
    //expDelLow_max = 20;  expDelHigh_max = 150
    HV = 31.00;
    index = find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV);
    DCR_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 10.8468;
    errDCR_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.0110521;
    CT_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.194892;
    DCR_Del_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 10.6742;
    errDCR_Del_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.0332437;
    CT_Del_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.206141;

    HV = 32.00;
    index = find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV);
    DCR_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 14.7845;
    errDCR_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.0131704;
    CT_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.20516;
    DCR_Del_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 15.1647;
    errDCR_Del_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.0290423;
    CT_Del_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.226057;

    HV = 33.00;
    index = find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV);
    DCR_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 17.7625;
    errDCR_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.0146536;
    CT_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.234807;
    DCR_Del_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 18.458;
    errDCR_Del_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.0277113;
    CT_Del_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.256393;

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
    CT_Del_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.366466;



    ///////////////////////////////////////////////////////////////////////////
    //      SiPM2
    ///////////////////////////////////////////////////////////////////////////
    // pe_0_5_vect[7] = {8, 8, 8, 9, 9, 9, 10}
    // pe_1_5_vect[7] = {19, 23, 26, 30, 34, 35, 40}
    //minyhistDelays = 20;  maxyhistDelays = 150
    //expDelLow_max = 20;  expDelHigh_max = 150
    HV = 31.00;
    index = find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV);
    DCR_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 10.9316;
    errDCR_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.0111002;
    CT_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.176788;
    DCR_Del_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 10.7255;
    errDCR_Del_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.0331164;
    CT_Del_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.194828;

    HV = 32.00;
    index = find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV);
    DCR_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 15.1275;
    errDCR_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.0133455;
    CT_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.189882;
    DCR_Del_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 15.407;
    errDCR_Del_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.028496;
    CT_Del_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.204558;

    HV = 33.00;
    index = find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV);
    DCR_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 18.8724;
    errDCR_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.0151872;
    CT_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.226319;
    DCR_Del_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 19.4985;
    errDCR_Del_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.0269202;
    CT_Del_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.244278;

    HV = 34.00;
    index = find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV);
    DCR_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 21.6926;
    errDCR_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.0165057;
    CT_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.254493;
    DCR_Del_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 22.601;
    errDCR_Del_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.0266066;
    CT_Del_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.265904;

    HV = 35.00;
    index = find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV);
    DCR_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 24.3031;
    errDCR_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.0176867;
    CT_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.281917;
    DCR_Del_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 25.3718;
    errDCR_Del_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.0265973;
    CT_Del_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.289128;

    HV = 36.00;
    index = find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV);
    DCR_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 26.7797;
    errDCR_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.0187785;
    CT_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.322969;
    DCR_Del_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 27.9802;
    errDCR_Del_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.026794;
    CT_Del_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.331032;

    HV = 37.00;
    index = find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV);
    DCR_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 29.0874;
    errDCR_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.0197751;
    CT_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.346014;
    DCR_Del_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 30.4092;
    errDCR_Del_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.0271947;
    CT_Del_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.348417;




    ///////////////////////////////////////////////////////////////////////////
    //      SiPM3
    ///////////////////////////////////////////////////////////////////////////

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
    errDCR_Del_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 0.0271953;
    CT_Del_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 0.341255;

    HV = 37.00;
    index = find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV);
    DCR_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 32.3161;
    errDCR_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 0.0211412;
    CT_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 0.369078;
    DCR_Del_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 33.8594;
    errDCR_Del_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 0.027641;
    CT_Del_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 0.366327;


    //##########################################################################
    //          FIND ERRORS
    //##########################################################################


    ///////////////////////////////////////////////////////////////////////////
    //      SiPM1 Ea
    ///////////////////////////////////////////////////////////////////////////
    // pe_0_5_vect[7] = {8, 8, 8, 8, 8, 8, 10}
    // pe_1_5_vect[7] = {21, 25, 28, 33, 38, 43, 47}
    //minyhistDelays = 20;  maxyhistDelays = 150
    //expDelLow_max = 20;  expDelHigh_max = 150
    HV = 31.00;
    index = find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV);
    DCR_1_Ea[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)] = 10.8468;
    errDCR_1_Ea[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)] = 0.0110521;
    CT_1_Ea[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)] = 0.159048;
    DCR_Del_1_Ea[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)] = 10.6742;
    errDCR_Del_1_Ea[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)] = 0.0332437;
    CT_Del_1_Ea[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)] = 0.15804;

    HV = 32.00;
    index = find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV);
    DCR_1_Ea[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)] = 14.7845;
    errDCR_1_Ea[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)] = 0.0131704;
    CT_1_Ea[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)] = 0.185495;
    DCR_Del_1_Ea[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)] = 15.1647;
    errDCR_Del_1_Ea[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)] = 0.0290423;
    CT_Del_1_Ea[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)] = 0.189544;

    HV = 33.00;
    index = find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV);
    DCR_1_Ea[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)] = 17.7625;
    errDCR_1_Ea[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)] = 0.0146536;
    CT_1_Ea[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)] = 0.228707;
    DCR_Del_1_Ea[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)] = 18.458;
    errDCR_Del_1_Ea[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)] = 0.0277113;
    CT_Del_1_Ea[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)] = 0.247977;

    HV = 34.00;
    index = find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV);
    DCR_1_Ea[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)] = 20.3256;
    errDCR_1_Ea[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)] = 0.0158728;
    CT_1_Ea[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)] = 0.255342;
    DCR_Del_1_Ea[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)] = 21.186;
    errDCR_Del_1_Ea[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)] = 0.027144;
    CT_Del_1_Ea[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)] = 0.261282;

    HV = 35.00;
    index = find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV);
    DCR_1_Ea[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)] = 22.7297;
    errDCR_1_Ea[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)] = 0.0169789;
    CT_1_Ea[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)] = 0.278252;
    DCR_Del_1_Ea[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)] = 23.741;
    errDCR_Del_1_Ea[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)] = 0.0269559;
    CT_Del_1_Ea[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)] = 0.279835;

    HV = 36.00;
    index = find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV);
    DCR_1_Ea[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)] = 24.9651;
    errDCR_1_Ea[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)] = 0.017981;
    CT_1_Ea[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)] = 0.301545;
    DCR_Del_1_Ea[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)] = 26.0412;
    errDCR_Del_1_Ea[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)] = 0.0269684;
    CT_Del_1_Ea[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)] = 0.295666;

    HV = 37.00;
    index = find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV);
    DCR_1_Ea[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)] = 27.2241;
    errDCR_1_Ea[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)] = 0.0189719;
    CT_1_Ea[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)] = 0.329982;
    DCR_Del_1_Ea[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)] = 28.4227;
    errDCR_Del_1_Ea[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)] = 0.02718;
    CT_Del_1_Ea[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)] = 0.320734;


    ///////////////////////////////////////////////////////////////////////////
    //      SiPM2 Ea
    ///////////////////////////////////////////////////////////////////////////
    // pe_0_5_vect[7] = {9, 10, 11, 12, 13, 15, 18}
    // pe_1_5_vect[7] = {18, 20, 25, 27, 30, 30, 32}
    //minyhistDelays = 20;  maxyhistDelays = 150
    //expDelLow_max = 20;  expDelHigh_max = 150
    HV = 31.00;
    index = find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV);
    DCR_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 10.6895;
    errDCR_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 0.0109625;
    CT_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 0.210911;
    DCR_Del_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 10.4833;
    errDCR_Del_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 0.0335341;
    CT_Del_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 0.232251;

    HV = 32.00;
    index = find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV);
    DCR_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 14.5568;
    errDCR_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 0.0130533;
    CT_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 0.249134;
    DCR_Del_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 14.8065;
    errDCR_Del_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 0.0289394;
    CT_Del_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 0.265598;

    HV = 33.00;
    index = find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV);
    DCR_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 18.0494;
    errDCR_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 0.0147924;
    CT_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 0.245794;
    DCR_Del_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 18.6533;
    errDCR_Del_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 0.0272456;
    CT_Del_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 0.267676;

    HV = 34.00;
    index = find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV);
    DCR_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 21.015;
    errDCR_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 0.0161934;
    CT_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 0.283325;
    DCR_Del_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 21.9092;
    errDCR_Del_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 0.0267194;
    CT_Del_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 0.298341;

    HV = 35.00;
    index = find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV);
    DCR_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 23.6327;
    errDCR_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 0.0173866;
    CT_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 0.310841;
    DCR_Del_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 24.6869;
    errDCR_Del_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 0.0266172;
    CT_Del_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 0.323234;

    HV = 36.00;
    index = find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV);
    DCR_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 25.8183;
    errDCR_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 0.0183576;
    CT_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 0.365755;
    DCR_Del_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 26.977;
    errDCR_Del_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 0.0265593;
    CT_Del_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 0.365745;

    HV = 37.00;
    index = find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV);
    DCR_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 27.6054;
    errDCR_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 0.0191372;
    CT_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 0.404789;
    DCR_Del_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 28.7942;
    errDCR_Del_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 0.0263734;
    CT_Del_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 0.402225;



    ///////////////////////////////////////////////////////////////////////////
    //      SiPM3 Ea
    ///////////////////////////////////////////////////////////////////////////
    // pe_0_5_vect[7] = {9, 10, 11, 12, 14, 15, 18}
    // pe_1_5_vect[7] = {18, 22, 25, 27, 30, 32, 35}
    //minyhistDelays = 20;  maxyhistDelays = 150
    //expDelLow_max = 20;  expDelHigh_max = 150
    HV = 31.00;
    index = find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV);
    DCR_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = 13.1021;
    errDCR_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = 0.0122915;
    CT_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = 0.232156;
    DCR_Del_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = 12.9489;
    errDCR_Del_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = 0.0299836;
    CT_Del_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = 0.216705;

    HV = 32.00;
    index = find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV);
    DCR_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = 17.3636;
    errDCR_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = 0.0144594;
    CT_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = 0.228256;
    DCR_Del_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = 17.7919;
    errDCR_Del_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = 0.0272692;
    CT_Del_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = 0.244039;

    HV = 33.00;
    index = find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV);
    DCR_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = 20.9634;
    errDCR_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = 0.0161695;
    CT_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = 0.258495;
    DCR_Del_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = 21.7458;
    errDCR_Del_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = 0.0265984;
    CT_Del_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = 0.270949;

    HV = 34.00;
    index = find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV);
    DCR_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = 24.0479;
    errDCR_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = 0.0175726;
    CT_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = 0.295889;
    DCR_Del_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = 25.0715;
    errDCR_Del_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = 0.0265189;
    CT_Del_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = 0.302424;

    HV = 35.00;
    index = find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV);
    DCR_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = 26.6753;
    errDCR_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = 0.018733;
    CT_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = 0.326226;
    DCR_Del_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = 27.7887;
    errDCR_Del_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = 0.0265087;
    CT_Del_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = 0.330481;

    HV = 36.00;
    index = find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV);
    DCR_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = 29.279;
    errDCR_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = 0.019857;
    CT_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = 0.36154;
    DCR_Del_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = 30.5843;
    errDCR_Del_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = 0.0268489;
    CT_Del_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = 0.360884;

    HV = 37.00;
    index = find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV);
    DCR_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = 31.0595;
    errDCR_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = 0.0206131;
    CT_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = 0.40121;
    DCR_Del_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = 32.4955;
    errDCR_Del_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = 0.0267661;
    CT_Del_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = 0.399899;


    //-------------------------------------------------------------------------

    ///////////////////////////////////////////////////////////////////////////
    //      SiPM1 Eb
    ///////////////////////////////////////////////////////////////////////////
    // pe_0_5_vect[7] = {9, 10, 11, 12, 13, 15, 18}
    // pe_1_5_vect[7] = {18, 22, 25, 27, 30, 33, 35}
    //minyhistDelays = 20;  maxyhistDelays = 150
    //expDelLow_max = 20;  expDelHigh_max = 150
    HV = 31.00;
    index = find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV);
    DCR_1_Eb[find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV)] = 10.7205;
    errDCR_1_Eb[find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV)] = 0.0109802;
    CT_1_Eb[find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV)] = 0.225367;
    DCR_Del_1_Eb[find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV)] = 10.5428;
    errDCR_Del_1_Eb[find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV)] = 0.033473;
    CT_Del_1_Eb[find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV)] = 0.220219;

    HV = 32.00;
    index = find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV);
    DCR_1_Eb[find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV)] = 14.4917;
    errDCR_1_Eb[find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV)] = 0.0130198;
    CT_1_Eb[find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV)] = 0.220962;
    DCR_Del_1_Eb[find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV)] = 14.8765;
    errDCR_Del_1_Eb[find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV)] = 0.0292807;
    CT_Del_1_Eb[find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV)] = 0.239717;

    HV = 33.00;
    index = find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV);
    DCR_1_Eb[find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV)] = 17.5137;
    errDCR_1_Eb[find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV)] = 0.0145326;
    CT_1_Eb[find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV)] = 0.250471;
    DCR_Del_1_Eb[find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV)] = 18.215;
    errDCR_Del_1_Eb[find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV)] = 0.0277938;
    CT_Del_1_Eb[find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV)] = 0.270254;

    HV = 34.00;
    index = find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV);
    DCR_1_Eb[find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV)] = 20.0712;
    errDCR_1_Eb[find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV)] = 0.0157538;
    CT_1_Eb[find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV)] = 0.289202;
    DCR_Del_1_Eb[find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV)] = 20.9357;
    errDCR_Del_1_Eb[find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV)] = 0.0271786;
    CT_Del_1_Eb[find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV)] = 0.302458;

    HV = 35.00;
    index = find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV);
    DCR_1_Eb[find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV)] = 22.4171;
    errDCR_1_Eb[find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV)] = 0.0168369;
    CT_1_Eb[find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV)] = 0.316738;
    DCR_Del_1_Eb[find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV)] = 23.4197;
    errDCR_Del_1_Eb[find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV)] = 0.026935;
    CT_Del_1_Eb[find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV)] = 0.331369;

    HV = 36.00;
    index = find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV);
    DCR_1_Eb[find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV)] = 24.4549;
    errDCR_1_Eb[find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV)] = 0.0177543;
    CT_1_Eb[find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV)] = 0.347784;
    DCR_Del_1_Eb[find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV)] = 25.4725;
    errDCR_Del_1_Eb[find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV)] = 0.0268055;
    CT_Del_1_Eb[find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV)] = 0.361245;

    HV = 37.00;
    index = find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV);
    DCR_1_Eb[find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV)] = 26.2797;
    errDCR_1_Eb[find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV)] = 0.0185601;
    CT_1_Eb[find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV)] = 0.386231;
    DCR_Del_1_Eb[find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV)] = 27.3814;
    errDCR_Del_1_Eb[find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV)] = 0.0266365;
    CT_Del_1_Eb[find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV)] = 0.395822;


    ///////////////////////////////////////////////////////////////////////////
    //      SiPM2 Eb
    ///////////////////////////////////////////////////////////////////////////
    // pe_0_5_vect[7] = {8, 8, 8, 8, 8, 8, 10}
    // pe_1_5_vect[7] = {21, 25, 28, 32, 37, 40, 45}
    //minyhistDelays = 20;  maxyhistDelays = 150
    //expDelLow_max = 20;  expDelHigh_max = 150
    HV = 31.00;
    index = find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV);
    DCR_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 10.9316;
    errDCR_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 0.0111002;
    CT_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 0.131453;
    DCR_Del_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 10.7255;
    errDCR_Del_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 0.0331164;
    CT_Del_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 0.136572;

    HV = 32.00;
    index = find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV);
    DCR_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 15.1275;
    errDCR_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 0.0133455;
    CT_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 0.158579;
    DCR_Del_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 15.407;
    errDCR_Del_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 0.028496;
    CT_Del_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 0.155831;

    HV = 33.00;
    index = find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV);
    DCR_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 18.8724;
    errDCR_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 0.0151872;
    CT_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 0.20467;
    DCR_Del_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 19.4985;
    errDCR_Del_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 0.0269202;
    CT_Del_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 0.208428;

    HV = 34.00;
    index = find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV);
    DCR_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 21.7974;
    errDCR_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 0.0165538;
    CT_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 0.236793;
    DCR_Del_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 22.6931;
    errDCR_Del_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 0.0261976;
    CT_Del_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 0.237478;

    HV = 35.00;
    index = find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV);
    DCR_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 24.3871;
    errDCR_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 0.0177241;
    CT_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 0.257462;
    DCR_Del_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 25.4487;
    errDCR_Del_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 0.0265585;
    CT_Del_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 0.250932;

    HV = 36.00;
    index = find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV);
    DCR_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 26.8525;
    errDCR_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 0.0188102;
    CT_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 0.294755;
    DCR_Del_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 28.0447;
    errDCR_Del_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 0.0267643;
    CT_Del_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 0.289827;

    HV = 37.00;
    index = find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV);
    DCR_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 29.0874;
    errDCR_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 0.0197751;
    CT_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 0.316284;
    DCR_Del_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 30.4092;
    errDCR_Del_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 0.0271947;
    CT_Del_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 0.303547;



    ///////////////////////////////////////////////////////////////////////////
    //      SiPM3 Eb
    ///////////////////////////////////////////////////////////////////////////
    // pe_0_5_vect[7] = {8, 8, 8, 8, 8, 8, 10}
    // pe_1_5_vect[7] = {21, 25, 28, 32, 37, 40, 45}
    //minyhistDelays = 20;  maxyhistDelays = 150
    //expDelLow_max = 20;  expDelHigh_max = 150
    HV = 31.00;
    index = find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV);
    DCR_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 13.2422;
    errDCR_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 0.012366;
    CT_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 0.165705;
    DCR_Del_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 13.1078;
    errDCR_Del_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 0.0298438;
    CT_Del_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 0.156278;

    HV = 32.00;
    index = find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV);
    DCR_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 17.6465;
    errDCR_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 0.0145972;
    CT_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 0.194027;
    DCR_Del_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 18.0967;
    errDCR_Del_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 0.0271804;
    CT_Del_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 0.194214;

    HV = 33.00;
    index = find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV);
    DCR_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 21.2565;
    errDCR_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 0.016305;
    CT_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 0.235449;
    DCR_Del_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 22.0518;
    errDCR_Del_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 0.0265947;
    CT_Del_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 0.244405;

    HV = 34.00;
    index = find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV);
    DCR_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 24.3689;
    errDCR_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 0.017716;
    CT_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 0.266216;
    DCR_Del_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 25.4166;
    errDCR_Del_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 0.026592;
    CT_Del_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 0.270176;

    HV = 35.00;
    index = find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV);
    DCR_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 27.3095;
    errDCR_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 0.019009;
    CT_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 0.289045;
    DCR_Del_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 28.5021;
    errDCR_Del_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 0.0267917;
    CT_Del_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 0.283831;

    HV = 36.00;
    index = find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV);
    DCR_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 30.1011;
    errDCR_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 0.0202073;
    CT_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 0.320699;
    DCR_Del_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 31.4716;
    errDCR_Del_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 0.0272803;
    CT_Del_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 0.315756;

    HV = 37.00;
    index = find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV);
    DCR_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 32.8109;
    errDCR_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 0.021348;
    CT_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 0.341771;
    DCR_Del_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 34.4032;
    errDCR_Del_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 0.0279002;
    CT_Del_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 0.328943;



    //------------------------------




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


    // TEST Z
    int N_Z = n_DCR;
    double testZ_DCR[3][N_Z];
    double testZ_CT[3][N_Z];
    int n_SiPM = 0;
    int n_SiPM_tot = 3;

    n_SiPM = 0;
    for(int i=0; i<N_Z; i++){
        testZ_DCR[n_SiPM][i] = TMath::Abs((DCR_1[i] - DCR_Del_1[i]) / (TMath::Sqrt( errDCR_1[i]*errDCR_1[i] + errDCR_Del_1[i]*errDCR_Del_1[i] )));
    }

    n_SiPM = 1;
    for(int i=0; i<N_Z; i++){
        testZ_DCR[n_SiPM][i] = TMath::Abs((DCR_2[i] - DCR_Del_2[i]) / (TMath::Sqrt( errDCR_2[i]*errDCR_2[i] + errDCR_Del_2[i]*errDCR_Del_2[i] )));
    }

    n_SiPM = 2;
    for(int i=0; i<N_Z; i++){
        testZ_DCR[n_SiPM][i] = TMath::Abs((DCR_2[i] - DCR_Del_2[i]) / (TMath::Sqrt( errDCR_2[i]*errDCR_2[i] + errDCR_Del_2[i]*errDCR_Del_2[i] )));
    }

    n_SiPM = 0;
    for(int i=0; i<N_Z; i++){
        testZ_CT[n_SiPM][i] = TMath::Abs((CT_1[i] - CT_Del_1[i]) / (TMath::Sqrt( errCT_1[i]*errCT_1[i] + errCT_Del_1[i]*errCT_Del_1[i] )));
    }

    n_SiPM = 1;
    for(int i=0; i<N_Z; i++){
        testZ_CT[n_SiPM][i] = TMath::Abs((CT_2[i] - CT_Del_2[i]) / (TMath::Sqrt( errCT_2[i]*errCT_2[i] + errCT_Del_2[i]*errCT_Del_2[i] )));
    }

    n_SiPM = 2;
    for(int i=0; i<N_Z; i++){
        testZ_CT[n_SiPM][i] = TMath::Abs((CT_2[i] - CT_Del_2[i]) / (TMath::Sqrt( errCT_2[i]*errCT_2[i] + errCT_Del_2[i]*errCT_Del_2[i] )));
    }

    cout<<endl;
    cout<<"TEST Z"<<endl;
    for(int n=0; n<n_SiPM_tot; n++){
        cout<<"SiPM "<<n+1<<endl;
        for(int i=0; i<N_Z; i++){
            cout<<testZ_DCR[n][i];
            if(testZ_DCR[n][i] > 1.96) cout<<"\tTEST NOT PASSED"<<endl;
            else cout<<endl;

            cout<<testZ_CT[n][i];
            if(testZ_CT[n][i] > 1.96) cout<<"\tTEST NOT PASSED "<<endl;
            else cout<<endl;
        }
        cout<<endl;
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


    bool draw_check_bool = true;

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



    ////////////////////////////////////////////////////////////////////////////
    //          OLD
    ////////////////////////////////////////////////////////////////////////////

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
    double err;

    err = max(TMath::Abs(num-LimUp), TMath::Abs(num-LimDown));

    if(err > num*max_percentage){
        err = num*max_percentage;
    }else{
        if(err < num*min_percentage){
            err = num*min_percentage;
        }
    }

    return err;

}
