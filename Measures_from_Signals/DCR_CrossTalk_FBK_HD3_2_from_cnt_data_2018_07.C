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
 *  > for HV = 32 ... 37:
 *    minyhistDelays = 15;  maxyhistDelays = 100;
 *    expDelLow_max  = minyhistDelays*1.25; expDelHigh_max = maxyhistDelays;
 *
 *  > for HV = 31: (I need a bit more points)
 *    minyhistDelays = 15;  maxyhistDelays = 127;
 *    expDelLow_max  = minyhistDelays*1.25; expDelHigh_max = maxyhistDelays;
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
    bool fix_error_bool = true;

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



    ///////////////////////////////////////////////////////////////////////////
    //      SiPM1
    ///////////////////////////////////////////////////////////////////////////
    HV = 31.00;
    index = find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV);
    DCR_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 10.8468;
    errDCR_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.0110521;
    CT_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.176791;
    DCR_Del_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 9.75963;
    errDCR_Del_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.0398614;
    CT_Del_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.140421;



    HV = 32.00;
    index = find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV);
    DCR_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 14.7845;
    errDCR_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.0131704;
    CT_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.20516;
    DCR_Del_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 13.6261;
    errDCR_Del_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.0470757;
    CT_Del_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.17856;

    HV = 33.00;
    index = find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV);
    DCR_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 17.7624;
    errDCR_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.0146535;
    CT_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.240605;
    DCR_Del_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 17.1121;
    errDCR_Del_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.0428226;
    CT_Del_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.229532;

    HV = 34.00;
    index = find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV);
    DCR_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 20.2242;
    errDCR_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.0158254;
    CT_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.266764;
    DCR_Del_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 19.8707;
    errDCR_Del_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.0405474;
    CT_Del_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.26092;

    HV = 35.00;
    index = find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV);
    DCR_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 22.5769;
    errDCR_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.0169096;
    CT_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.295021;
    DCR_Del_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 22.4412;
    errDCR_Del_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.0390127;
    CT_Del_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.288113;

    HV = 36.00;
    index = find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV);
    DCR_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 24.7735;
    errDCR_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.017896;
    CT_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.320079;
    DCR_Del_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 24.7027;
    errDCR_Del_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.037829;
    CT_Del_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.306099;

    HV = 37.00;
    index = find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV);
    DCR_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 27.0179;
    errDCR_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.0188823;
    CT_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.358409;
    DCR_Del_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 27.065;
    errDCR_Del_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.0369479;
    CT_Del_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.356612;


    ///////////////////////////////////////////////////////////////////////////
    //      SiPM2
    ///////////////////////////////////////////////////////////////////////////
    HV = 31.00;
    index = find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV);
    DCR_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 10.9316;
    errDCR_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.0111002;
    CT_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.154444;
    DCR_Del_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 10.1127;
    errDCR_Del_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.0399677;
    CT_Del_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.111517;



    HV = 32.00;
    index = find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV);
    DCR_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 14.5568;
    errDCR_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.0130533;
    CT_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.182738;
    DCR_Del_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 14.2059;
    errDCR_Del_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.0502267;
    CT_Del_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.165761;



    HV = 33.00;
    index = find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV);
    DCR_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 18.4081;
    errDCR_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.014965;
    CT_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.241006;
    DCR_Del_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 17.6203;
    errDCR_Del_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.041782;
    CT_Del_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.233664;

    HV = 34.00;
    index = find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV);
    DCR_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 21.2632;
    errDCR_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.0163081;
    CT_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.259632;
    DCR_Del_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 20.8746;
    errDCR_Del_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.0393913;
    CT_Del_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.248935;

    HV = 35.00;
    index = find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV);
    DCR_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 23.3732;
    errDCR_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.0172698;
    CT_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.303938;
    DCR_Del_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 22.9998;
    errDCR_Del_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.0379174;
    CT_Del_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.308938;

    HV = 36.00;
    index = find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV);
    DCR_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 25.8183;
    errDCR_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.0183576;
    CT_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.334995;
    DCR_Del_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 25.5563;
    errDCR_Del_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.036735;
    CT_Del_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.335642;

    HV = 37.00;
    index = find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV);
    DCR_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 28.4578;
    errDCR_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.019505;
    CT_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.365338;
    DCR_Del_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 28.4297;
    errDCR_Del_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.0360842;
    CT_Del_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.365736;



    ///////////////////////////////////////////////////////////////////////////
    //      SiPM3
    ///////////////////////////////////////////////////////////////////////////
    HV = 31.00;
    index = find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV);
    DCR_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 13.2422;
    errDCR_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 0.012366;
    CT_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 0.183176;
    DCR_Del_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 12.1804;
    errDCR_Del_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 0.0353279;
    CT_Del_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 0.149937;



    HV = 32.00;
    index = find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV);
    DCR_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 17.6465;
    errDCR_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 0.0145972;
    CT_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 0.213261;
    DCR_Del_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 16.3993;
    errDCR_Del_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 0.0423881;
    CT_Del_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 0.19709;

    HV = 33.00;
    index = find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV);
    DCR_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 21.2565;
    errDCR_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 0.016305;
    CT_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 0.248154;
    DCR_Del_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 20.6269;
    errDCR_Del_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 0.0393206;
    CT_Del_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 0.243895;

    HV = 34.00;
    index = find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV);
    DCR_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 24.272;
    errDCR_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 0.0176727;
    CT_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 0.272539;
    DCR_Del_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 23.9895;
    errDCR_Del_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 0.0376377;
    CT_Del_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 0.266177;

    HV = 35.00;
    index = find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV);
    DCR_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 27.1553;
    errDCR_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 0.018942;
    CT_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 0.300822;
    DCR_Del_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 26.9617;
    errDCR_Del_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 0.0363447;
    CT_Del_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 0.289481;

    HV = 36.00;
    index = find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV);
    DCR_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 29.8447;
    errDCR_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 0.0200983;
    CT_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 0.323454;
    DCR_Del_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 29.8948;
    errDCR_Del_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 0.0358119;
    CT_Del_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 0.302173;

    HV = 37.00;
    index = find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV);
    DCR_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 32.4947;
    errDCR_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 0.0212159;
    CT_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 0.36373;
    DCR_Del_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 32.7241;
    errDCR_Del_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 0.0353436;
    CT_Del_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 0.350756;


    //------------------------------


    // PERCENTAGE ERROR
    if(percentage_error_bool){
        double err_rel = 0.05;
        for(int i=0; i<n_DCR; i++){
            errDCR_1[i] = err_rel * DCR_1[i];
        }
        for(int i=0; i<n_DCR; i++){
            errDCR_2[i] = err_rel * DCR_2[i];
        }
        for(int i=0; i<n_DCR; i++){
            errDCR_3[i] = err_rel * DCR_3[i];
        }
        for(int i=0; i<n_DCR; i++){
            errCT_1[i] = err_rel * CT_1[i];
        }
        for(int i=0; i<n_DCR; i++){
            errCT_2[i] = err_rel * CT_2[i];
        }
        for(int i=0; i<n_DCR; i++){
            errCT_3[i] = err_rel * CT_3[i];
        }
        for(int i=0; i<n_DCR; i++){
            errDCR_Del_1[i] = err_rel * DCR_Del_1[i];
        }
        for(int i=0; i<n_DCR; i++){
            errDCR_Del_2[i] = err_rel * DCR_Del_2[i];
        }
        for(int i=0; i<n_DCR; i++){
            errDCR_Del_3[i] = err_rel * DCR_Del_3[i];
        }
        for(int i=0; i<n_DCR; i++){
            errCT_Del_1[i] = err_rel * CT_Del_1[i];
        }
        for(int i=0; i<n_DCR; i++){
            errCT_Del_2[i] = err_rel * CT_Del_2[i];
        }
        for(int i=0; i<n_DCR; i++){
            errCT_Del_3[i] = err_rel * CT_Del_3[i];
        }
    }

    //------------------------------

    // FIXED ERROR
    if(fix_error_bool){
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
    }


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


    //##########################################################################
    //          FIND ERRORS
    //##########################################################################

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
    // pe_0_5_vect[7] = {8, 8, 8, 8, 8, 8, 10}
    // pe_1_5_vect[7] = {21, 25, 30, 35, 40, 45, 50}
    //minyhistDelays = 20;  maxyhistDelays = 200
    //expDelLow_max = 20;  expDelHigh_max = 200
    HV = 31.00;
    index = find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV);
    DCR_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 10.9316;
    errDCR_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 0.0111002;
    CT_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 0.131453;
    DCR_Del_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 11.602;
    errDCR_Del_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 0.0240823;
    CT_Del_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 0.171492;

    HV = 32.00;
    index = find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV);
    DCR_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 15.1275;
    errDCR_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 0.0133455;
    CT_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 0.158579;
    DCR_Del_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 16.2447;
    errDCR_Del_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 0.0221862;
    CT_Del_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 0.187555;

    HV = 33.00;
    index = find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV);
    DCR_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 18.8724;
    errDCR_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 0.0151872;
    CT_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 0.175978;
    DCR_Del_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 20.3384;
    errDCR_Del_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 0.0220846;
    CT_Del_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 0.19224;

    HV = 34.00;
    index = find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV);
    DCR_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 21.7974;
    errDCR_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 0.0165538;
    CT_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 0.198008;
    DCR_Del_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 23.5415;
    errDCR_Del_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 0.0226857;
    CT_Del_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 0.203376;

    HV = 35.00;
    index = find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV);
    DCR_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 24.3871;
    errDCR_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 0.0177241;
    CT_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 0.215636;
    DCR_Del_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 26.3411;
    errDCR_Del_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 0.0234478;
    CT_Del_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 0.212713;

    HV = 36.00;
    index = find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV);
    DCR_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 26.8525;
    errDCR_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 0.0188102;
    CT_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 0.232341;
    DCR_Del_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 29.0031;
    errDCR_Del_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 0.0242963;
    CT_Del_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 0.221849;

    HV = 37.00;
    index = find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV);
    DCR_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 29.0874;
    errDCR_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 0.0197751;
    CT_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 0.25166;
    DCR_Del_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 31.3833;
    errDCR_Del_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 0.0253658;
    CT_Del_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 0.238954;


    ///////////////////////////////////////////////////////////////////////////
    //      SiPM3 Ea
    ///////////////////////////////////////////////////////////////////////////
    // pe_0_5_vect[7] = {8, 8, 8, 8, 8, 8, 10}
    // pe_1_5_vect[7] = {21, 25, 30, 35, 40, 45, 50}
    //minyhistDelays = 20;  maxyhistDelays = 200
    //expDelLow_max = 20;  expDelHigh_max = 200
    HV = 31.00;
    index = find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV);
    DCR_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = 13.2422;
    errDCR_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = 0.012366;
    CT_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = 0.165705;
    DCR_Del_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = 14.0369;
    errDCR_Del_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = 0.0225116;
    CT_Del_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = 0.193459;

    HV = 32.00;
    index = find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV);
    DCR_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = 17.6465;
    errDCR_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = 0.0145972;
    CT_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = 0.194027;
    DCR_Del_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = 18.9897;
    errDCR_Del_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = 0.0219023;
    CT_Del_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = 0.214116;

    HV = 33.00;
    index = find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV);
    DCR_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = 21.2565;
    errDCR_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = 0.016305;
    CT_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = 0.219144;
    DCR_Del_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = 22.9336;
    errDCR_Del_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = 0.0225356;
    CT_Del_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = 0.227474;

    HV = 34.00;
    index = find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV);
    DCR_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = 24.3689;
    errDCR_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = 0.017716;
    CT_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = 0.244875;
    DCR_Del_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = 26.2717;
    errDCR_Del_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = 0.023497;
    CT_Del_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = 0.246811;

    HV = 35.00;
    index = find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV);
    DCR_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = 27.3095;
    errDCR_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = 0.019009;
    CT_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = 0.265673;
    DCR_Del_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = 29.457;
    errDCR_Del_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = 0.024497;
    CT_Del_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = 0.261155;

    HV = 36.00;
    index = find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV);
    DCR_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = 30.1011;
    errDCR_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = 0.0202073;
    CT_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = 0.283635;
    DCR_Del_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = 32.525;
    errDCR_Del_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = 0.0256689;
    CT_Del_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = 0.272759;

    HV = 37.00;
    index = find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV);
    DCR_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = 32.8109;
    errDCR_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = 0.021348;
    CT_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = 0.301727;
    DCR_Del_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = 35.5253;
    errDCR_Del_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = 0.0270075;
    CT_Del_3_Ea[find_index(HV_3_Ea,  sizeof(HV_3_Ea)/sizeof(HV_3_Ea[0]), HV)] = 0.283775;

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
    // pe_0_5_vect[7] = {9, 10, 11, 12, 13, 15, 18}
    // pe_1_5_vect[7] = {18, 22, 29, 30, 32, 35, 30}
    //minyhistDelays = 20;  maxyhistDelays = 200
    //expDelLow_max = 20;  expDelHigh_max = 200
    HV = 31.00;
    index = find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV);
    DCR_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 10.6895;
    errDCR_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 0.0109625;
    CT_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 0.210911;
    DCR_Del_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 11.353;
    errDCR_Del_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 0.0243102;
    CT_Del_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 0.254267;

    HV = 32.00;
    index = find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV);
    DCR_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 14.5568;
    errDCR_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 0.0130533;
    CT_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 0.212621;
    DCR_Del_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 15.642;
    errDCR_Del_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 0.0223518;
    CT_Del_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 0.252214;

    HV = 33.00;
    index = find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV);
    DCR_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 18.0494;
    errDCR_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 0.0147924;
    CT_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 0.200895;
    DCR_Del_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 19.4685;
    errDCR_Del_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 0.0221294;
    CT_Del_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 0.222964;

    HV = 34.00;
    index = find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV);
    DCR_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 21.015;
    errDCR_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 0.0161934;
    CT_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 0.262698;
    DCR_Del_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 22.7308;
    errDCR_Del_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 0.0226296;
    CT_Del_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 0.285649;

    HV = 35.00;
    index = find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV);
    DCR_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 23.6327;
    errDCR_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 0.0173866;
    CT_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 0.3006;
    DCR_Del_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 25.5558;
    errDCR_Del_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 0.0233162;
    CT_Del_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 0.320003;

    HV = 36.00;
    index = find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV);
    DCR_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 25.8183;
    errDCR_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 0.0183576;
    CT_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 0.334995;
    DCR_Del_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 27.9323;
    errDCR_Del_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 0.023866;
    CT_Del_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 0.346388;

    HV = 37.00;
    index = find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV);
    DCR_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 27.6054;
    errDCR_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 0.0191372;
    CT_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 0.440531;
    DCR_Del_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 29.8953;
    errDCR_Del_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 0.0241753;
    CT_Del_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 0.430482;


    ///////////////////////////////////////////////////////////////////////////
    //      SiPM3 Eb
    ///////////////////////////////////////////////////////////////////////////
    // pe_0_5_vect[7] = {9, 10, 11, 12, 13, 15, 18}
    // pe_1_5_vect[7] = {18, 22, 29, 30, 32, 35, 30}
    //minyhistDelays = 20;  maxyhistDelays = 200
    //expDelLow_max = 20;  expDelHigh_max = 200
    HV = 31.00;
    index = find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV);
    DCR_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 13.1021;
    errDCR_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 0.0122915;
    CT_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 0.232156;
    DCR_Del_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 13.88;
    errDCR_Del_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 0.0225726;
    CT_Del_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 0.253223;

    HV = 32.00;
    index = find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV);
    DCR_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 17.3636;
    errDCR_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 0.0144594;
    CT_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 0.228256;
    DCR_Del_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 18.6832;
    errDCR_Del_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 0.0218945;
    CT_Del_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 0.254529;

    HV = 33.00;
    index = find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV);
    DCR_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 20.9634;
    errDCR_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 0.0161695;
    CT_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 0.231824;
    DCR_Del_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 22.6361;
    errDCR_Del_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 0.0224599;
    CT_Del_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 0.246736;

    HV = 34.00;
    index = find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV);
    DCR_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 24.0479;
    errDCR_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 0.0175726;
    CT_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 0.279325;
    DCR_Del_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 25.9453;
    errDCR_Del_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 0.0233416;
    CT_Del_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 0.295508;

    HV = 35.00;
    index = find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV);
    DCR_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 26.898;
    errDCR_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 0.0188301;
    CT_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 0.314934;
    DCR_Del_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 29.0269;
    errDCR_Del_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 0.0242327;
    CT_Del_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 0.327042;

    HV = 36.00;
    index = find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV);
    DCR_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 29.279;
    errDCR_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 0.019857;
    CT_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 0.349254;
    DCR_Del_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 31.6833;
    errDCR_Del_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 0.0250331;
    CT_Del_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 0.35539;

    HV = 37.00;
    index = find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV);
    DCR_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 31.0595;
    errDCR_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 0.0206131;
    CT_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 0.473316;
    DCR_Del_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 33.7397;
    errDCR_Del_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 0.0253609;
    CT_Del_3_Eb[find_index(HV_3_Eb,  sizeof(HV_3_Eb)/sizeof(HV_3_Eb[0]), HV)] = 0.441148;


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

    TCanvas *cV_DCR_1_E = new TCanvas("cV_DCR_1_E", "cV_DCR_1_E",w,h);
    cV_DCR_1_E->SetGrid();
    TMultiGraph *mgV_DCR_1_E = new TMultiGraph("mgV_DCR_1_E", ";Bias Voltage (V); DCR");
    mgV_DCR_1_E->Add(gV_DCR_1_Ea);
    mgV_DCR_1_E->Add(gV_DCR_1_Eb);
    mgV_DCR_1_E->Draw("AP");

    TCanvas *cV_DCR_2_E = new TCanvas("cV_DCR_2_E", "cV_DCR_2_E",w,h);
    cV_DCR_2_E->SetGrid();
    TMultiGraph *mgV_DCR_2_E = new TMultiGraph("mgV_DCR_2_E", ";Bias Voltage (V); DCR");
    mgV_DCR_2_E->Add(gV_DCR_2_Ea);
    mgV_DCR_2_E->Add(gV_DCR_2_Eb);
    mgV_DCR_2_E->Draw("AP");

    TCanvas *cV_DCR_3_E = new TCanvas("cV_DCR_3_E", "cV_DCR_3_E",w,h);
    cV_DCR_3_E->SetGrid();
    TMultiGraph *mgV_DCR_3_E = new TMultiGraph("mgV_DCR_3_E", ";Bias Voltage (V); DCR");
    mgV_DCR_3_E->Add(gV_DCR_3_Ea);
    mgV_DCR_3_E->Add(gV_DCR_3_Eb);
    mgV_DCR_3_E->Draw("AP");

    //----------

    TCanvas *cV_CT_1_E = new TCanvas("cV_CT_1_E", "cV_CT_1_E",w,h);
    cV_CT_1_E->SetGrid();
    TMultiGraph *mgV_CT_1_E = new TMultiGraph("mgV_CT_1_E", ";Bias Voltage (V); CT");
    mgV_CT_1_E->Add(gV_CT_1_Ea);
    mgV_CT_1_E->Add(gV_CT_1_Eb);
    mgV_CT_1_E->Draw("AP");

    TCanvas *cV_CT_2_E = new TCanvas("cV_CT_2_E", "cV_CT_2_E",w,h);
    cV_CT_2_E->SetGrid();
    TMultiGraph *mgV_CT_2_E = new TMultiGraph("mgV_CT_2_E", ";Bias Voltage (V); CT");
    mgV_CT_2_E->Add(gV_CT_2_Ea);
    mgV_CT_2_E->Add(gV_CT_2_Eb);
    mgV_CT_2_E->Draw("AP");

    TCanvas *cV_CT_3_E = new TCanvas("cV_CT_3_E", "cV_CT_3_E",w,h);
    cV_CT_3_E->SetGrid();
    TMultiGraph *mgV_CT_3_E = new TMultiGraph("mgV_CT_3_E", ";Bias Voltage (V); CT");
    mgV_CT_3_E->Add(gV_CT_3_Ea);
    mgV_CT_3_E->Add(gV_CT_3_Eb);
    mgV_CT_3_E->Draw("AP");

    //----------

    TCanvas *cV_DCR_Del_1_E = new TCanvas("cV_DCR_Del_1_E", "cV_DCR_Del_1_E",w,h);
    cV_DCR_Del_1_E->SetGrid();
    TMultiGraph *mgV_DCR_Del_1_E = new TMultiGraph("mgV_DCR_Del_1_E", ";Bias Voltage (V); DCR_Del");
    mgV_DCR_Del_1_E->Add(gV_DCR_Del_1_Ea);
    mgV_DCR_Del_1_E->Add(gV_DCR_Del_1_Eb);
    mgV_DCR_Del_1_E->Draw("AP");

    TCanvas *cV_DCR_Del_2_E = new TCanvas("cV_DCR_Del_2_E", "cV_DCR_Del_2_E",w,h);
    cV_DCR_Del_2_E->SetGrid();
    TMultiGraph *mgV_DCR_Del_2_E = new TMultiGraph("mgV_DCR_Del_2_E", ";Bias Voltage (V); DCR_Del");
    mgV_DCR_Del_2_E->Add(gV_DCR_Del_2_Ea);
    mgV_DCR_Del_2_E->Add(gV_DCR_Del_2_Eb);
    mgV_DCR_Del_2_E->Draw("AP");

    TCanvas *cV_DCR_Del_3_E = new TCanvas("cV_DCR_Del_3_E", "cV_DCR_Del_3_E",w,h);
    cV_DCR_Del_3_E->SetGrid();
    TMultiGraph *mgV_DCR_Del_3_E = new TMultiGraph("mgV_DCR_Del_3_E", ";Bias Voltage (V); DCR_Del");
    mgV_DCR_Del_3_E->Add(gV_DCR_Del_3_Ea);
    mgV_DCR_Del_3_E->Add(gV_DCR_Del_3_Eb);
    mgV_DCR_Del_3_E->Draw("AP");

    //----------

    TCanvas *cV_CT_Del_1_E = new TCanvas("cV_CT_Del_1_E", "cV_CT_Del_1_E",w,h);
    cV_CT_Del_1_E->SetGrid();
    TMultiGraph *mgV_CT_Del_1_E = new TMultiGraph("mgV_CT_Del_1_E", ";Bias Voltage (V); CT_Del");
    mgV_CT_Del_1_E->Add(gV_CT_Del_1_Ea);
    mgV_CT_Del_1_E->Add(gV_CT_Del_1_Eb);
    mgV_CT_Del_1_E->Draw("AP");

    TCanvas *cV_CT_Del_2_E = new TCanvas("cV_CT_Del_2_E", "cV_CT_Del_2_E",w,h);
    cV_CT_Del_2_E->SetGrid();
    TMultiGraph *mgV_CT_Del_2_E = new TMultiGraph("mgV_CT_Del_2_E", ";Bias Voltage (V); CT_Del");
    mgV_CT_Del_2_E->Add(gV_CT_Del_2_Ea);
    mgV_CT_Del_2_E->Add(gV_CT_Del_2_Eb);
    mgV_CT_Del_2_E->Draw("AP");

    TCanvas *cV_CT_Del_3_E = new TCanvas("cV_CT_Del_3_E", "cV_CT_Del_3_E",w,h);
    cV_CT_Del_3_E->SetGrid();
    TMultiGraph *mgV_CT_Del_3_E = new TMultiGraph("mgV_CT_Del_3_E", ";Bias Voltage (V); CT_Del");
    mgV_CT_Del_3_E->Add(gV_CT_Del_3_Ea);
    mgV_CT_Del_3_E->Add(gV_CT_Del_3_Eb);
    mgV_CT_Del_3_E->Draw("AP");




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
