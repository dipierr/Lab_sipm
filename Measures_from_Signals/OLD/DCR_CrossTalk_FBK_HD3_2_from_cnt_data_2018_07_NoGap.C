/******************************************************************************\
 * DCR_CrossTalk_FBK_HD3_2_from_cnt_data_2018_07_NoGap.C
 *
 * Values obtained by Ana_Traces_SiPM.cxx (version of 01/09/2018)
 *
 * KEY POINTS:
 *  > DCR_CT_1SiPM_nHVs(...)
 *  > dleddt = 6
 *  > NO trace smoothing
 *  > thr at 0.5pe and 1.5 pe set manually
 *  > min_thr_to_find_peaks = 8;  //first thr value in the DCR vs thr plot (mV)
 *  > max_thr_to_find_peaks = 80; //last thr value in the DCR vs thr plot (mV)
 *  > Area = 36; // 6x6 mm^2
 *  > No gap between peaks (when I found a peak I don't jump of a certain value
 *    but I only do i++, i.e. gap_between_peaks=1)
 *  > trace_time = SUM(trace_DLED[0][trace_DLED_length-1] - trace_DLED[0][0])/
 *                 / n_ev_tot
 *
 *  > minyhistDelays = 50;  maxyhistDelays = 250; bins_Delays = 100;
 *    expDelLow_max  = minyhistDelays; expDelHigh_max = maxyhistDelays;
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

#define n_DCR_1 7
#define n_DCR_2 7
#define n_DCR_3 7


#define h 600
#define w 1000

void DCR_CrossTalk_FBK_HD3_2_from_cnt_data_2018_07_NoGap();
int find_index(double v[],int N, double value);

char title_DCR[80];
char title_DCR_mg[80];

void DCR_CrossTalk_FBK_HD3_2_from_cnt_data_2018_07_NoGap(){

    // SiPM1:
    double HV_1[n_DCR_1], errHV_1[n_DCR_1];
    double DCR_1[n_DCR_1], errDCR_1[n_DCR_1];
    double CT_1[n_DCR_1], errCT_1[n_DCR_1];
    double DCR_Del_1[n_DCR_1], errDCR_Del_1[n_DCR_1];
    double CT_Del_1[n_DCR_1], errCT_Del_1[n_DCR_1];
    // SiPM2:
    double HV_2[n_DCR_2], errHV_2[n_DCR_2];
    double DCR_2[n_DCR_2], errDCR_2[n_DCR_2];
    double CT_2[n_DCR_2], errCT_2[n_DCR_2];
    double DCR_Del_2[n_DCR_2], errDCR_Del_2[n_DCR_2];
    double CT_Del_2[n_DCR_2], errCT_Del_2[n_DCR_2];
    // SiPM3:
    double HV_3[n_DCR_3], errHV_3[n_DCR_3];
    double DCR_3[n_DCR_3], errDCR_3[n_DCR_3];
    double CT_3[n_DCR_3], errCT_3[n_DCR_3];
    double DCR_Del_3[n_DCR_3], errDCR_Del_3[n_DCR_3];
    double CT_Del_3[n_DCR_3], errCT_Del_3[n_DCR_3];

    double HV = 0.;

    // index:
    int index = 0;

    // ERRORS:
    bool percentage_error_bool = false;
    bool fix_error_bool = true;

    // EVALUARE DCR / AREA
    bool dcr_area = true;

    // DRAW ALL
    bool draw_all_bool = true;

    // Area:
    double Area = 36;

    // Initialization
    for(int i=0; i<n_DCR_1; i++){
        HV_1[i] = errHV_1[i] = DCR_1[i] = errDCR_1[i] = CT_1[i] = errCT_1[i] = DCR_Del_1[i] = errDCR_Del_1[i] = CT_Del_1[i] = errCT_Del_1[i] = 0.;
    }
    for(int i=0; i<n_DCR_2; i++){
        HV_2[i] = errHV_2[i] = DCR_2[i] = errDCR_2[i] = CT_2[i] = errCT_2[i] = DCR_Del_2[i] = errDCR_Del_2[i] = CT_Del_2[i] = errCT_Del_2[i] = 0.;
    }
    for(int i=0; i<n_DCR_3; i++){
        HV_3[i] = errHV_3[i] = DCR_3[i] = errDCR_3[i] = CT_3[i] = errCT_3[i] = DCR_Del_3[i] = errDCR_Del_3[i] = CT_Del_3[i] = errCT_Del_3[i] = 0.;
    }

    // HV
    HV_1[0]    = 31.00;
    errHV_1[0] =  0.01;
    for(int i=1; i<n_DCR_1; i++){
        HV_1[i]    = HV_1[i-1]+1.;
        errHV_1[i] = errHV_1[0];
    }

    HV_2[0]    = 31.00;
    errHV_2[0] =  0.01;
    for(int i=1; i<n_DCR_2; i++){
        HV_2[i]    = HV_2[i-1]+1.;
        errHV_2[i] = errHV_2[0];
    }

    HV_3[0]    = 31.00;
    errHV_3[0] =  0.01;
    for(int i=1; i<n_DCR_3; i++){
        HV_3[i]    = HV_3[i-1]+1.;
        errHV_3[i] = errHV_3[0];
    }



    ///////////////////////////////////////////////////////////////////////////
    //      SiPM1
    ///////////////////////////////////////////////////////////////////////////
    // pe_0_5_vect[7] = {8, 8, 8, 10, 11, 12, 13}
    // pe_1_5_vect[7] = {20, 23, 26, 31, 35, 40, 41}
    HV = 31.00;
    index = find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV);
    DCR_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 9.95933;
    errDCR_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.00993211;
    CT_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.195298;
    DCR_Del_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 12.4792;
    errDCR_Del_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.0274071;
    CT_Del_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.247793;

    HV = 32.00;
    index = find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV);
    DCR_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 13.2638;
    errDCR_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.0114621;
    CT_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.23259;
    DCR_Del_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 16.682;
    errDCR_Del_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.0276985;
    CT_Del_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.260092;

    HV = 33.00;
    index = find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV);
    DCR_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 15.6956;
    errDCR_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.0124687;
    CT_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.276528;
    DCR_Del_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 19.7206;
    errDCR_Del_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.0293058;
    CT_Del_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.292369;

    HV = 34.00;
    index = find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV);
    DCR_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 17.6614;
    errDCR_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.0132266;
    CT_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.308771;
    DCR_Del_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 22.2145;
    errDCR_Del_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.0311016;
    CT_Del_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.313847;

    HV = 35.00;
    index = find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV);
    DCR_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 19.4671;
    errDCR_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.0138863;
    CT_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.343641;
    DCR_Del_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 24.5922;
    errDCR_Del_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.0332963;
    CT_Del_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.340381;

    HV = 36.00;
    index = find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV);
    DCR_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 21.0983;
    errDCR_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.0144565;
    CT_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.374096;
    DCR_Del_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 26.8789;
    errDCR_Del_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.0356678;
    CT_Del_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.363921;

    HV = 37.00;
    index = find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV);
    DCR_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 22.7339;
    errDCR_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.0150064;
    CT_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.421334;
    DCR_Del_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 29.2745;
    errDCR_Del_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.0386523;
    CT_Del_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.401166;


    ///////////////////////////////////////////////////////////////////////////
    //      SiPM2
    ///////////////////////////////////////////////////////////////////////////
    // pe_0_5_vect[7] = {8, 8, 8, 10, 11, 12, 13}
    // pe_1_5_vect[7] = {20, 23, 26, 31, 35, 40, 41}
    HV = 31.00;
    index = find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV);
    DCR_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 10.0594;
    errDCR_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.00998189;
    CT_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.170992;
    DCR_Del_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 12.6367;
    errDCR_Del_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.0274834;
    CT_Del_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.226124;

    HV = 32.00;
    index = find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV);
    DCR_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 13.5249;
    errDCR_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.0115744;
    CT_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.216435;
    DCR_Del_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 17.1982;
    errDCR_Del_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.0278527;
    CT_Del_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.246158;

    HV = 33.00;
    index = find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV);
    DCR_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 16.5454;
    errDCR_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.0128018;
    CT_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.262708;
    DCR_Del_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 21.2371;
    errDCR_Del_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.0300278;
    CT_Del_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.270279;

    HV = 34.00;
    index = find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV);
    DCR_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 18.605;
    errDCR_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.0135753;
    CT_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.292313;
    DCR_Del_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 23.8468;
    errDCR_Del_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.0323347;
    CT_Del_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.292139;

    HV = 35.00;
    index = find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV);
    DCR_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 20.498;
    errDCR_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.0142493;
    CT_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.32891;
    DCR_Del_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 26.4365;
    errDCR_Del_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.0350154;
    CT_Del_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.319597;

    HV = 36.00;
    index = find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV);
    DCR_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 22.2579;
    errDCR_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.0148485;
    CT_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.354199;
    DCR_Del_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 28.9937;
    errDCR_Del_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.0379824;
    CT_Del_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.338926;

    HV = 37.00;
    index = find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV);
    DCR_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 23.97;
    errDCR_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.0154091;
    CT_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.410107;
    DCR_Del_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 31.5214;
    errDCR_Del_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.0413399;
    CT_Del_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.387422;


    ///////////////////////////////////////////////////////////////////////////
    //      SiPM3
    ///////////////////////////////////////////////////////////////////////////
    // pe_0_5_vect[7] = {8, 8, 8, 10, 11, 12, 13}
    // pe_1_5_vect[7] = {20, 23, 26, 31, 35, 40, 41}
    HV = 31.00;
    index = find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV);
    DCR_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 11.9937;
    errDCR_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 0.0108995;
    CT_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 0.205997;
    DCR_Del_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 15.2157;
    errDCR_Del_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 0.0271581;
    CT_Del_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 0.233567;

    HV = 32.00;
    index = find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV);
    DCR_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 15.5797;
    errDCR_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 0.0124226;
    CT_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 0.246363;
    DCR_Del_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 20.022;
    errDCR_Del_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 0.0290943;
    CT_Del_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 0.25344;

    HV = 33.00;
    index = find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV);
    DCR_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 18.4252;
    errDCR_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 0.0135096;
    CT_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 0.291183;
    DCR_Del_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 23.6899;
    errDCR_Del_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 0.0320873;
    CT_Del_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 0.28557;

    HV = 34.00;
    index = find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV);
    DCR_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 20.7195;
    errDCR_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 0.0143261;
    CT_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 0.322676;
    DCR_Del_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 26.7714;
    errDCR_Del_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 0.0353575;
    CT_Del_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 0.310956;

    HV = 35.00;
    index = find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV);
    DCR_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 22.8275;
    errDCR_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 0.0150373;
    CT_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 0.358968;
    DCR_Del_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 29.9022;
    errDCR_Del_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 0.0389383;
    CT_Del_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 0.337003;

    HV = 36.00;
    index = find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV);
    DCR_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 24.7271;
    errDCR_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 0.0156505;
    CT_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 0.387899;
    DCR_Del_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 32.9017;
    errDCR_Del_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 0.0430857;
    CT_Del_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 0.357481;

    HV = 37.00;
    index = find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV);
    DCR_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 26.5393;
    errDCR_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 0.0162139;
    CT_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 0.438103;
    DCR_Del_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 35.8947;
    errDCR_Del_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 0.0478439;
    CT_Del_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 0.396995;


    //------------------------------


    // PERCENTAGE ERROR
    if(percentage_error_bool){
        double err_rel = 0.05;
        for(int i=0; i<n_DCR_1; i++){
            errDCR_1[i] = err_rel * DCR_1[i];
        }
        for(int i=0; i<n_DCR_2; i++){
            errDCR_2[i] = err_rel * DCR_2[i];
        }
        for(int i=0; i<n_DCR_3; i++){
            errDCR_3[i] = err_rel * DCR_3[i];
        }
        for(int i=0; i<n_DCR_1; i++){
            errCT_1[i] = err_rel * CT_1[i];
        }
        for(int i=0; i<n_DCR_2; i++){
            errCT_2[i] = err_rel * CT_2[i];
        }
        for(int i=0; i<n_DCR_3; i++){
            errCT_3[i] = err_rel * CT_3[i];
        }
        for(int i=0; i<n_DCR_1; i++){
            errDCR_Del_1[i] = err_rel * DCR_Del_1[i];
        }
        for(int i=0; i<n_DCR_1; i++){
            errDCR_Del_2[i] = err_rel * DCR_Del_2[i];
        }
        for(int i=0; i<n_DCR_1; i++){
            errDCR_Del_3[i] = err_rel * DCR_Del_3[i];
        }
        for(int i=0; i<n_DCR_1; i++){
            errCT_Del_1[i] = err_rel * CT_Del_1[i];
        }
        for(int i=0; i<n_DCR_1; i++){
            errCT_Del_2[i] = err_rel * CT_Del_2[i];
        }
        for(int i=0; i<n_DCR_1; i++){
            errCT_Del_3[i] = err_rel * CT_Del_3[i];
        }
    }

    //------------------------------

    // FIXED ERROR
    if(fix_error_bool){
        double err_fix_DCR = 1;
        double err_fix_CT = 0.02;
        for(int i=0; i<n_DCR_1; i++){
            errDCR_1[i] = err_fix_DCR;
        }
        for(int i=0; i<n_DCR_2; i++){
            errDCR_2[i] = err_fix_DCR;
        }
        for(int i=0; i<n_DCR_3; i++){
            errDCR_3[i] = err_fix_DCR;
        }
        for(int i=0; i<n_DCR_1; i++){
            errCT_1[i] = err_fix_CT;
        }
        for(int i=0; i<n_DCR_2; i++){
            errCT_2[i] = err_fix_CT;
        }
        for(int i=0; i<n_DCR_3; i++){
            errCT_3[i] = err_fix_CT;
        }
        for(int i=0; i<n_DCR_1; i++){
            errDCR_Del_1[i] = err_fix_DCR;
        }
        for(int i=0; i<n_DCR_1; i++){
            errDCR_Del_2[i] = err_fix_DCR;
        }
        for(int i=0; i<n_DCR_1; i++){
            errDCR_Del_3[i] = err_fix_DCR;
        }
        for(int i=0; i<n_DCR_1; i++){
            errCT_Del_1[i] = err_fix_CT;
        }
        for(int i=0; i<n_DCR_1; i++){
            errCT_Del_2[i] = err_fix_CT;
        }
        for(int i=0; i<n_DCR_1; i++){
            errCT_Del_3[i] = err_fix_CT;
        }
    }


    //------------------------------

    if(dcr_area){
        for(int i=0; i< n_DCR_1; i++){
            DCR_1[i]  /= Area;
            errDCR_1[i] /= Area;
            DCR_1[i]  *= 1e3;
            errDCR_1[i] *= 1e3;
        }
        for(int i=0; i< n_DCR_2; i++){
            DCR_2[i]  /= Area;
            errDCR_2[i] /= Area;
            DCR_2[i]  *= 1e3;
            errDCR_2[i] *= 1e3;
        }
        for(int i=0; i< n_DCR_3; i++){
            DCR_3[i]  /= Area;
            errDCR_3[i] /= Area;
            DCR_3[i]  *= 1e3;
            errDCR_3[i] *= 1e3;
        }

        //------------------------------

        for(int i=0; i< n_DCR_1; i++){
            DCR_Del_1[i]  /= Area;
            errDCR_Del_1[i] /= Area;
            DCR_Del_1[i]  *= 1e3;
            errDCR_Del_1[i] *= 1e3;
        }
        for(int i=0; i< n_DCR_2; i++){
            DCR_Del_2[i]  /= Area;
            errDCR_Del_2[i] /= Area;
            DCR_Del_2[i]  *= 1e3;
            errDCR_Del_2[i] *= 1e3;
        }
        for(int i=0; i< n_DCR_3; i++){
            DCR_Del_3[i]  /= Area;
            errDCR_Del_3[i] /= Area;
            DCR_Del_3[i]  *= 1e3;
            errDCR_Del_3[i] *= 1e3;
        }

    }


    //------------------------------
    //------------------------------


    TGraphErrors *gV_DCR_1  = new TGraphErrors(n_DCR_1, HV_1, DCR_1, errHV_1, errDCR_1);
    TGraphErrors *gV_DCR_2  = new TGraphErrors(n_DCR_2, HV_2, DCR_2, errHV_2, errDCR_2);
    TGraphErrors *gV_DCR_3  = new TGraphErrors(n_DCR_3, HV_3, DCR_3, errHV_3, errDCR_3);

    TGraphErrors *gV_CT_1  = new TGraphErrors(n_DCR_1, HV_1, CT_1, errHV_1, errCT_1);
    TGraphErrors *gV_CT_2  = new TGraphErrors(n_DCR_2, HV_2, CT_2, errHV_2, errCT_2);
    TGraphErrors *gV_CT_3  = new TGraphErrors(n_DCR_3, HV_3, CT_3, errHV_3, errCT_3);

    TGraphErrors *gV_DCR_Del_1  = new TGraphErrors(n_DCR_1, HV_1, DCR_Del_1, errHV_1, errDCR_Del_1);
    TGraphErrors *gV_DCR_Del_2  = new TGraphErrors(n_DCR_2, HV_2, DCR_Del_2, errHV_2, errDCR_Del_2);
    TGraphErrors *gV_DCR_Del_3  = new TGraphErrors(n_DCR_3, HV_3, DCR_Del_3, errHV_3, errDCR_3);

    TGraphErrors *gV_CT_Del_1  = new TGraphErrors(n_DCR_1, HV_1, CT_Del_1, errHV_1, errCT_Del_1);
    TGraphErrors *gV_CT_Del_2  = new TGraphErrors(n_DCR_2, HV_2, CT_Del_2, errHV_2, errCT_Del_2);
    TGraphErrors *gV_CT_Del_3  = new TGraphErrors(n_DCR_3, HV_3, CT_Del_3, errHV_3, errCT_Del_3);


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

    // if(draw_all_bool){
    //     TCanvas *cV_DCR_1 = new TCanvas("cV_DCR_1", "cV_DCR_1",w,h);
    //     cV_DCR_1->SetGrid();
    //     gV_DCR_1->Draw("AP");
    //
    //     TCanvas *cV_DCR_2 = new TCanvas("cV_DCR_2", "cV_DCR_2",w,h);
    //     cV_DCR_2->SetGrid();
    //     gV_DCR_2->Draw("AP");
    //
    //     TCanvas *cV_DCR_3 = new TCanvas("cV_DCR_3", "cV_DCR_3",w,h);
    //     cV_DCR_3->SetGrid();
    //     gV_DCR_3->Draw("AP");
    // }

    //------------------------------

    // if(draw_all_bool){
    //     TCanvas *cV_CT_1 = new TCanvas("cV_CT_1", "cV_CT_1",w,h);
    //     cV_CT_1->SetGrid();
    //     gV_CT_1->Draw("AP");
    //
    //     TCanvas *cV_CT_2 = new TCanvas("cV_CT_2", "cV_CT_2",w,h);
    //     cV_CT_2->SetGrid();
    //     gV_CT_2->Draw("AP");
    //
    //     TCanvas *cV_CT_3 = new TCanvas("cV_CT_3", "cV_CT_3",w,h);
    //     cV_CT_3->SetGrid();
    //     gV_CT_3->Draw("AP");
    // }


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

    // if(draw_all_bool){
    //     TCanvas *cV_DCR_Del_1 = new TCanvas("cV_DCR_Del_1", "cV_DCR_Del_1",w,h);
    //     cV_DCR_Del_1->SetGrid();
    //     gV_DCR_Del_1->Draw("AP");
    //
    //     TCanvas *cV_DCR_Del_2 = new TCanvas("cV_DCR_Del_2", "cV_DCR_Del_2",w,h);
    //     cV_DCR_Del_2->SetGrid();
    //     gV_DCR_Del_2->Draw("AP");
    //
    //     TCanvas *cV_DCR_Del_3 = new TCanvas("cV_DCR_Del_3", "cV_DCR_Del_3",w,h);
    //     cV_DCR_Del_3->SetGrid();
    //     gV_DCR_Del_3->Draw("AP");
    // }


    //------------------------------

    // if(draw_all_bool){
    //     TCanvas *cV_CT_Del_1 = new TCanvas("cV_CT_Del_1", "cV_CT_Del_1",w,h);
    //     cV_CT_Del_1->SetGrid();
    //     gV_CT_Del_1->Draw("AP");
    //
    //     TCanvas *cV_CT_Del_2 = new TCanvas("cV_CT_Del_2", "cV_CT_Del_2",w,h);
    //     cV_CT_Del_2->SetGrid();
    //     gV_CT_Del_2->Draw("AP");
    //
    //     TCanvas *cV_CT_Del_3 = new TCanvas("cV_CT_Del_3", "cV_CT_Del_3",w,h);
    //     cV_CT_Del_3->SetGrid();
    //     gV_CT_Del_3->Draw("AP");
    // }


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
    int N_Z = n_DCR_1;
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
    // for(int i=0; i<n_DCR_1; i++){
    //     printf("$ %.2f \\pm %.2f $ & $ %.2f \\pm %.2f $ & $ %.2f \\pm %.2f $ & $ %.2f \\pm %.2f $ & $ %.2f \\pm %.2f $ & $ %.2f \\pm %.2f $ & $ %.2f \\pm %.2f $ \\\\ \n", HV_1[i], errHV_1[i], DCR_1[i], errDCR_1[i],DCR_2[i], errDCR_2[i],DCR_3[i], errDCR_3[i], CT_1[i], errCT_1[i],CT_2[i], errCT_2[i],CT_3[i], errCT_3[i]);
    // }
    cout<<endl<<endl;
    cout<<"-------------------------------------------------"<<endl;
    cout<<endl<<endl;
    cout<<"For LaTeX"<<endl<<endl;
    cout<<"\% From file DCR_CrossTalk_FBK_HD3_2_from_cnt_data_2018_07.C"<<endl;
    cout<<"\% From cnt"<<endl;
    cout<<"\% HV & DCR_1 (kHz/mm^2) & DCR_2 (kHz/mm^2) & DCR_3 (kHz/mm^2) \\\\"<<endl;
    for(int i=0; i<n_DCR_1; i++){
        printf("$ %.2f \\pm %.2f $ & $ %.2f \\pm %.2f $ & $ %.2f \\pm %.2f $ & $ %.2f \\pm %.2f $  \\\\ \n", HV_1[i], errHV_1[i], DCR_1[i], errDCR_1[i],DCR_2[i], errDCR_2[i],DCR_3[i], errDCR_3[i]);
    }
    cout<<endl<<endl;
    cout<<"\% From file DCR_CrossTalk_FBK_HD3_2_from_cnt_data_2018_07.C"<<endl;
    cout<<"\% From cnt"<<endl;
    cout<<"\% HV & CT_1 & CT_2  & CT_3  \\\\"<<endl;
    for(int i=0; i<n_DCR_1; i++){
        printf("$ %.2f \\pm %.2f $ & $ %.2f \\pm %.2f $ & $ %.2f \\pm %.2f $ & $ %.2f \\pm %.2f $  \\\\ \n", HV_1[i], errHV_1[i], CT_1[i], errCT_1[i],CT_2[i], errCT_2[i],CT_3[i], errCT_3[i]);
    }

    cout<<endl<<endl;
    cout<<"-------------------------------------------------"<<endl;
    cout<<endl<<endl;
    cout<<"\% From file DCR_CrossTalk_FBK_HD3_2_from_cnt_data_2018_07.C"<<endl;
    cout<<"\% From Del"<<endl;
    cout<<"\% HV & DCR_Del_1 (kHz/mm^2) & DCR_Del_2 (kHz/mm^2) & DCR_Del_3 (kHz/mm^2) \\\\"<<endl;
    for(int i=0; i<n_DCR_1; i++){
        printf("$ %.2f \\pm %.2f $ & $ %.2f \\pm %.2f $ & $ %.2f \\pm %.2f $ & $ %.2f \\pm %.2f $  \\\\ \n", HV_1[i], errHV_1[i], DCR_Del_1[i], errDCR_Del_1[i],DCR_Del_2[i], errDCR_Del_2[i],DCR_Del_3[i], errDCR_Del_3[i]);
    }
    cout<<endl<<endl;
    cout<<"\% From file DCR_CrossTalk_FBK_HD3_2_from_cnt_data_2018_07.C"<<endl;
    cout<<"\% From Del"<<endl;
    cout<<"\% HV & CT_Del_1 & CT_Del_2  & CT_Del_3  \\\\"<<endl;
    for(int i=0; i<n_DCR_1; i++){
        printf("$ %.2f \\pm %.2f $ & $ %.2f \\pm %.2f $ & $ %.2f \\pm %.2f $ & $ %.2f \\pm %.2f $  \\\\ \n", HV_1[i], errHV_1[i], CT_Del_1[i], errCT_Del_1[i],CT_Del_2[i], errCT_Del_2[i],CT_Del_3[i], errCT_Del_3[i]);
    }




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
