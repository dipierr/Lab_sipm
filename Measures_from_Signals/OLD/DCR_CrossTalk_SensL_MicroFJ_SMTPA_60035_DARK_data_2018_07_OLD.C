/******************************************************************************\
 *  DCR_CrossTalk_SensL_MicroFJ_SMTPA_60035_DARK_data_2018_07.C
 *
 *  Values obtained by Ana_Traces_SiPM.cxx (version of 18/08/2018, 1)
 *
 * KEY POINTS:
 *  > DCR_CT_1SiPM_nHVs(...)
 *  > dleddt = 8
 *  > NO Trace smoothing
 *  > thr at 0.5pe and 1.5 pe set manually
 *  > min_thr_to_find_peaks = 7.1;//first thr value in the DCR vs thr plot (mV)
 *  > max_thr_to_find_peaks = 60; //last thr value in the DCR vs thr plot (mV)
 *  > Area = 36.844900; // 6.07*6.07 mm^2
 *
 *  FILES ANALYZED:
 *
 * 20180725_SensL_MicroFJ-SMTPA-60035_01_STD_DARK_AgilentE3641A_28.00_
 *        Adapt.SensL.AS.02_AS_2_100000ev_01.dat
 * 20180725_SensL_MicroFJ-SMTPA-60035_01_STD_DARK_AgilentE3641A_29.00_
 *        Adapt.SensL.AS.02_AS_2_100000ev_01.dat
 * 20180725_SensL_MicroFJ-SMTPA-60035_01_STD_DARK_AgilentE3641A_30.00_
 *        Adapt.SensL.AS.02_AS_2_100000ev_01.dat
 * 20180725_SensL_MicroFJ-SMTPA-60035_01_STD_DARK_AgilentE3641A_31.00_
 *        Adapt.SensL.AS.02_AS_2_100000ev_01.dat
 *
 * 20180725_SensL_MicroFJ-SMTPA-60035_02_STD_DARK_AgilentE3641A_28.00_
 *        Adapt.SensL.AS.02_AS_2_100000ev_01.dat
 * 20180725_SensL_MicroFJ-SMTPA-60035_02_STD_DARK_AgilentE3641A_29.00_
 *        Adapt.SensL.AS.02_AS_2_100000ev_01.dat
 * 20180725_SensL_MicroFJ-SMTPA-60035_02_STD_DARK_AgilentE3641A_30.00_
 *        Adapt.SensL.AS.02_AS_2_100000ev_01.dat
 * 20180725_SensL_MicroFJ-SMTPA-60035_02_STD_DARK_AgilentE3641A_31.00_
 *        Adapt.SensL.AS.02_AS_2_100000ev_01.dat
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

#define n_DCR 4
#define n_DCR 4


#define h 600
#define w 1000

void DCR_CrossTalk_SensL_MicroFJ_SMTPA_60035_DARK_data_2018_07();
int find_index(double v[],int N, double value);
double GetMaxCheckPercentage(double num, double LimUp, double LimDown, double min_percentage, double max_percentage);

char title_DCR[80];
char title_DCR_mg[80];

void DCR_CrossTalk_SensL_MicroFJ_SMTPA_60035_DARK_data_2018_07(){

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

    double min_percentage = 0.;
    double max_percentage = 0.1;

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
    double Area;
    Area =  36.844900;

    // Initialization
    for(int i=0; i<n_DCR; i++){
        HV_1[i] = errHV_1[i] = DCR_1[i] = errDCR_1[i] = CT_1[i] = errCT_1[i] = DCR_Del_1[i] = errDCR_Del_1[i] = CT_Del_1[i] = errCT_Del_1[i] = 0.;
    }
    for(int i=0; i<n_DCR; i++){
        HV_2[i] = errHV_2[i] = DCR_2[i] = errDCR_2[i] = CT_2[i] = errCT_2[i] = DCR_Del_2[i] = errDCR_Del_2[i] = CT_Del_2[i] = errCT_Del_2[i] = 0.;
    }

    // HV
    HV_1[0]    = 28.00;
    errHV_1[0] =  0.01;
    for(int i=1; i<n_DCR; i++){
        HV_1[i]    = HV_1[i-1]+1.;
        errHV_1[i] = errHV_1[0];
    }

    HV_2[0]    = 28.00;
    errHV_2[0] =  0.01;
    for(int i=1; i<n_DCR; i++){
        HV_2[i]    = HV_2[i-1]+1.;
        errHV_2[i] = errHV_2[0];
    }

    for(int i=0; i<n_DCR; i++){
        HV_1_Ea[i] = HV_1_Eb[i] = HV_1[i];
        HV_2_Ea[i] = HV_2_Eb[i] = HV_2[i];
    }




    ///////////////////////////////////////////////////////////////////////////
    //      SiPM1
    ///////////////////////////////////////////////////////////////////////////
    // pe_0_5_vect[4] = {8, 9, 10, 12}
    // pe_1_5_vect[4] = {16, 20, 23, 26}
    //minyhistDelays = 20;  maxyhistDelays = 150
    //expDelLow_max = 20;  expDelHigh_max = 150
    HV = 28.00;
    index = find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV);
    DCR_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 3.82583;
    errDCR_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.00637308;
    CT_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.183342;
    DCR_Del_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 3.40796;
    errDCR_Del_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.0742061;
    CT_Del_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.518048;

    HV = 29.00;
    index = find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV);
    DCR_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 6.05288;
    errDCR_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.00814971;
    CT_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.217608;
    DCR_Del_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 6.45785;
    errDCR_Del_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.0524259;
    CT_Del_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.451736;

    HV = 30.00;
    index = find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV);
    DCR_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 8.10417;
    errDCR_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.00957019;
    CT_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.285169;
    DCR_Del_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 9.0646;
    errDCR_Del_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.0442492;
    CT_Del_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.503234;

    HV = 31.00;
    index = find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV);
    DCR_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 10.3189;
    errDCR_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.0109671;
    CT_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.370261;
    DCR_Del_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 12.0701;
    errDCR_Del_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.041697;
    CT_Del_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.673024;


    //////// Ea ////////
    HV = 28.00;
    index = find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV);
    DCR_1_Ea[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)] = DCR_1[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)];
    errDCR_1_Ea[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)] = errDCR_1[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)];
    CT_1_Ea[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)] = CT_1[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)];
    DCR_Del_1_Ea[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)] = DCR_Del_1[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)];
    errDCR_Del_1_Ea[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)] = errDCR_Del_1[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)];
    CT_Del_1_Ea[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)] = CT_Del_1[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)];

    HV = 29.00;
    index = find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV);
    DCR_1_Ea[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)] = 6.15159;
    errDCR_1_Ea[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)] = 0.00822181;
    CT_1_Ea[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)] = 0.227681;
    DCR_Del_1_Ea[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)] = 6.6401;
    errDCR_Del_1_Ea[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)] = 0.0518691;
    CT_Del_1_Ea[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)] = 0.477769;

    HV = 30.00;
    index = find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV);
    DCR_1_Ea[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)] = 8.19421;
    errDCR_1_Ea[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)] = 0.00962935;
    CT_1_Ea[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)] = 0.293682;
    DCR_Del_1_Ea[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)] = 9.30066;
    errDCR_Del_1_Ea[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)] = 0.0440958;
    CT_Del_1_Ea[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)] = 0.53053;

    HV = 31.00;
    index = find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV);
    DCR_1_Ea[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)] = 10.5312;
    errDCR_1_Ea[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)] = 0.0110955;
    CT_1_Ea[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)] = 0.378103;
    DCR_Del_1_Ea[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)] = 12.454;
    errDCR_Del_1_Ea[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)] = 0.0412398;
    CT_Del_1_Ea[find_index(HV_1_Ea,  sizeof(HV_1_Ea)/sizeof(HV_1_Ea[0]), HV)] = 0.688589;


    //////// Eb ////////
    HV = 28.00;
    index = find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV);
    DCR_1_Eb[find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV)] = 3.54204;
    errDCR_1_Eb[find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV)] = 0.00611903;
    CT_1_Eb[find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV)] = 0.171595;
    DCR_Del_1_Eb[find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV)] = 3.09716;
    errDCR_Del_1_Eb[find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV)] = 0.0793107;
    CT_Del_1_Eb[find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV)] = 0.485327;

    HV = 29.00;
    index = find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV);
    DCR_1_Eb[find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV)] = 5.8196;
    errDCR_1_Eb[find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV)] = 0.0079775;
    CT_1_Eb[find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV)] = 0.207349;
    DCR_Del_1_Eb[find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV)] = 6.03054;
    errDCR_Del_1_Eb[find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV)] = 0.0540351;
    CT_Del_1_Eb[find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV)] = 0.448159;

    HV = 30.00;
    index = find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV);
    DCR_1_Eb[find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV)] = 7.96154;
    errDCR_1_Eb[find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV)] = 0.00947601;
    CT_1_Eb[find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV)] = 0.277069;
    DCR_Del_1_Eb[find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV)] = 8.70437;
    errDCR_Del_1_Eb[find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV)] = 0.0447069;
    CT_Del_1_Eb[find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV)] = 0.4682;

    HV = 31.00;
    index = find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV);
    DCR_1_Eb[find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV)] = 10.0332;
    errDCR_1_Eb[find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV)] = 0.010793;
    CT_1_Eb[find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV)] = 0.364252;
    DCR_Del_1_Eb[find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV)] = 11.6481;
    errDCR_Del_1_Eb[find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV)] = 0.0424915;
    CT_Del_1_Eb[find_index(HV_1_Eb,  sizeof(HV_1_Eb)/sizeof(HV_1_Eb[0]), HV)] = 0.657914;




    ///////////////////////////////////////////////////////////////////////////
    //      SiPM2
    ///////////////////////////////////////////////////////////////////////////
    // pe_0_5_vect[4] = {8, 8, 9, 10}
    // pe_1_5_vect[4] = {22, 26, 30, 36}
    //minyhistDelays = 20;  maxyhistDelays = 150
    //expDelLow_max = 20;  expDelHigh_max = 150
    HV = 28.00;
    index = find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV);
    DCR_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 5.68451;
    errDCR_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.00787657;
    CT_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.122127;
    DCR_Del_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 5.08447;
    errDCR_Del_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.0551179;
    CT_Del_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.170022;

    HV = 29.00;
    index = find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV);
    DCR_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 7.46106;
    errDCR_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.00914068;
    CT_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.188607;
    DCR_Del_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 7.89427;
    errDCR_Del_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.0464612;
    CT_Del_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.313993;

    HV = 30.00;
    index = find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV);
    DCR_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 9.09404;
    errDCR_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.0102087;
    CT_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.263666;
    DCR_Del_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 10.343;
    errDCR_Del_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.0419293;
    CT_Del_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.409825;

    HV = 31.00;
    index = find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV);
    DCR_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 11.6313;
    errDCR_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.0117482;
    CT_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.332936;
    DCR_Del_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 13.4533;
    errDCR_Del_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.0386315;
    CT_Del_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.509661;


    //////// Ea ////////
    HV = 28.00;
    index = find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV);
    DCR_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = DCR_2[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)];
    errDCR_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = errDCR_2[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)];
    CT_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = CT_2[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)];
    DCR_Del_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = DCR_Del_2[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)];
    errDCR_Del_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = errDCR_Del_2[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)];
    CT_Del_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = CT_Del_2[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)];

    HV = 29.00;
    index = find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV);
    DCR_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = DCR_2[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)];
    errDCR_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = errDCR_2[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)];
    CT_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = CT_2[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)];
    DCR_Del_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = DCR_Del_2[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)];
    errDCR_Del_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = errDCR_Del_2[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)];
    CT_Del_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = CT_Del_2[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)];

    HV = 30.00;
    index = find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV);
    DCR_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 9.2344;
    errDCR_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 0.0102972;
    CT_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 0.266998;
    DCR_Del_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 10.4042;
    errDCR_Del_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 0.0416963;
    CT_Del_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 0.433206;

    HV = 31.00;
    index = find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV);
    DCR_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 11.7479;
    errDCR_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 0.0118162;
    CT_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 0.337408;
    DCR_Del_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 13.6162;
    errDCR_Del_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 0.0386876;
    CT_Del_2_Ea[find_index(HV_2_Ea,  sizeof(HV_2_Ea)/sizeof(HV_2_Ea[0]), HV)] = 0.52897;


    //////// Eb ////////
    HV = 28.00;
    index = find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV);
    DCR_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 5.39919;
    errDCR_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 0.00766027;
    CT_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 0.118032;
    DCR_Del_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 4.96559;
    errDCR_Del_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 0.0570655;
    CT_Del_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 0.162745;

    HV = 29.00;
    index = find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV);
    DCR_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 7.23414;
    errDCR_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 0.00898599;
    CT_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 0.187606;
    DCR_Del_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 7.84738;
    errDCR_Del_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 0.0472589;
    CT_Del_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 0.279413;

    HV = 30.00;
    index = find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV);
    DCR_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 8.96586;
    errDCR_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 0.0101274;
    CT_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 0.261656;
    DCR_Del_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 10.2109;
    errDCR_Del_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 0.0420182;
    CT_Del_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 0.384838;

    HV = 31.00;
    index = find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV);
    DCR_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 11.5344;
    errDCR_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 0.0116915;
    CT_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 0.327286;
    DCR_Del_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 13.2512;
    errDCR_Del_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 0.0385537;
    CT_Del_2_Eb[find_index(HV_2_Eb,  sizeof(HV_2_Eb)/sizeof(HV_2_Eb[0]), HV)] = 0.485201;






    // ERRORS
    for(int i=0; i<n_DCR; i++){
        errDCR_1[i] = GetMaxCheckPercentage(DCR_1[i], DCR_1_Ea[i], DCR_1_Eb[i], min_percentage, max_percentage);
        errDCR_2[i] = GetMaxCheckPercentage(DCR_2[i], DCR_2_Ea[i], DCR_2_Eb[i], min_percentage, max_percentage);
        errDCR_Del_1[i] = GetMaxCheckPercentage(DCR_Del_1[i], DCR_Del_1_Ea[i], DCR_Del_1_Eb[i], min_percentage, max_percentage);
        errDCR_Del_2[i] = GetMaxCheckPercentage(DCR_Del_2[i], DCR_Del_2_Ea[i], DCR_Del_2_Eb[i], min_percentage, max_percentage);
        errCT_1[i] = GetMaxCheckPercentage(CT_1[i], CT_1_Ea[i], CT_1_Eb[i], min_percentage, max_percentage);
        errCT_2[i] = GetMaxCheckPercentage(CT_2[i], CT_2_Ea[i], CT_2_Eb[i], min_percentage, max_percentage);
        errCT_Del_1[i] = GetMaxCheckPercentage(CT_Del_1[i], CT_Del_1_Ea[i], CT_Del_1_Eb[i], min_percentage, max_percentage);
        errCT_Del_2[i] = GetMaxCheckPercentage(CT_Del_2[i], CT_Del_2_Ea[i], CT_Del_2_Eb[i], min_percentage, max_percentage);
    } // end error from vect


    //------------------------------



    if(dcr_area){
        for(int i=0; i< n_DCR; i++){
            DCR_1[i]  /= Area;
            errDCR_1[i] /= Area;
            DCR_1[i]  *= 1e3;
            errDCR_1[i] *= 1e3;
            DCR_2[i]  /= Area;
            errDCR_2[i] /= Area;
            DCR_2[i]  *= 1e3;
            errDCR_2[i] *= 1e3;
            DCR_Del_1[i]  /= Area;
            errDCR_Del_1[i] /= Area;
            DCR_Del_1[i]  *= 1e3;
            errDCR_Del_1[i] *= 1e3;
            DCR_Del_2[i]  /= Area;
            errDCR_Del_2[i] /= Area;
            DCR_Del_2[i]  *= 1e3;
            errDCR_Del_2[i] *= 1e3;

            DCR_1_Ea[i]  /= Area;
            errDCR_1_Ea[i] /= Area;
            DCR_1_Ea[i]  *= 1e3;
            errDCR_1_Ea[i] *= 1e3;
            DCR_2_Ea[i]  /= Area;
            errDCR_2_Ea[i] /= Area;
            DCR_2_Ea[i]  *= 1e3;
            errDCR_2_Ea[i] *= 1e3;
            DCR_Del_1_Ea[i]  /= Area;
            errDCR_Del_1_Ea[i] /= Area;
            DCR_Del_1_Ea[i]  *= 1e3;
            errDCR_Del_1_Ea[i] *= 1e3;
            DCR_Del_2_Ea[i]  /= Area;
            errDCR_Del_2_Ea[i] /= Area;
            DCR_Del_2_Ea[i]  *= 1e3;
            errDCR_Del_2_Ea[i] *= 1e3;

            DCR_1_Eb[i]  /= Area;
            errDCR_1_Eb[i] /= Area;
            DCR_1_Eb[i]  *= 1e3;
            errDCR_1_Eb[i] *= 1e3;
            DCR_2_Eb[i]  /= Area;
            errDCR_2_Eb[i] /= Area;
            DCR_2_Eb[i]  *= 1e3;
            errDCR_2_Eb[i] *= 1e3;
            DCR_Del_1_Eb[i]  /= Area;
            errDCR_Del_1_Eb[i] /= Area;
            DCR_Del_1_Eb[i]  *= 1e3;
            errDCR_Del_1_Eb[i] *= 1e3;
            DCR_Del_2_Eb[i]  /= Area;
            errDCR_Del_2_Eb[i] /= Area;
            DCR_Del_2_Eb[i]  *= 1e3;
            errDCR_Del_2_Eb[i] *= 1e3;
        }


    }


    //------------------------------
    //------------------------------


    TGraphErrors *gV_DCR_1  = new TGraphErrors(n_DCR, HV_1, DCR_1, errHV_1, errDCR_1);
    TGraphErrors *gV_DCR_2  = new TGraphErrors(n_DCR, HV_2, DCR_2, errHV_2, errDCR_2);
    TGraphErrors *gV_CT_1  = new TGraphErrors(n_DCR, HV_1, CT_1, errHV_1, errCT_1);
    TGraphErrors *gV_CT_2  = new TGraphErrors(n_DCR, HV_2, CT_2, errHV_2, errCT_2);
    TGraphErrors *gV_DCR_Del_1  = new TGraphErrors(n_DCR, HV_1, DCR_Del_1, errHV_1, errDCR_Del_1);
    TGraphErrors *gV_DCR_Del_2  = new TGraphErrors(n_DCR, HV_2, DCR_Del_2, errHV_2, errDCR_Del_2);
    TGraphErrors *gV_CT_Del_1  = new TGraphErrors(n_DCR, HV_1, CT_Del_1, errHV_1, errCT_Del_1);
    TGraphErrors *gV_CT_Del_2  = new TGraphErrors(n_DCR, HV_2, CT_Del_2, errHV_2, errCT_Del_2);

    TGraphErrors *gV_DCR_1_Ea  = new TGraphErrors(n_DCR, HV_1, DCR_1_Ea, errHV_1, errDCR_1_Ea);
    TGraphErrors *gV_DCR_2_Ea  = new TGraphErrors(n_DCR, HV_2_Ea, DCR_2_Ea, errHV_2_Ea, errDCR_2_Ea);
    TGraphErrors *gV_CT_1_Ea  = new TGraphErrors(n_DCR, HV_1, CT_1_Ea, errHV_1, errCT_1_Ea);
    TGraphErrors *gV_CT_2_Ea  = new TGraphErrors(n_DCR, HV_2_Ea, CT_2_Ea, errHV_2_Ea, errCT_2_Ea);
    TGraphErrors *gV_DCR_Del_1_Ea  = new TGraphErrors(n_DCR, HV_1, DCR_Del_1_Ea, errHV_1, errDCR_Del_1_Ea);
    TGraphErrors *gV_DCR_Del_2_Ea  = new TGraphErrors(n_DCR, HV_2_Ea, DCR_Del_2_Ea, errHV_2_Ea, errDCR_Del_2_Ea);
    TGraphErrors *gV_CT_Del_1_Ea  = new TGraphErrors(n_DCR, HV_1, CT_Del_1_Ea, errHV_1, errCT_Del_1_Ea);
    TGraphErrors *gV_CT_Del_2_Ea  = new TGraphErrors(n_DCR, HV_2_Ea, CT_Del_2_Ea, errHV_2_Ea, errCT_Del_2_Ea);

    TGraphErrors *gV_DCR_1_Eb  = new TGraphErrors(n_DCR, HV_1, DCR_1_Eb, errHV_1, errDCR_1_Eb);
    TGraphErrors *gV_DCR_2_Eb  = new TGraphErrors(n_DCR, HV_2_Eb, DCR_2_Eb, errHV_2_Eb, errDCR_2_Eb);
    TGraphErrors *gV_CT_1_Eb  = new TGraphErrors(n_DCR, HV_1, CT_1_Eb, errHV_1, errCT_1_Eb);
    TGraphErrors *gV_CT_2_Eb  = new TGraphErrors(n_DCR, HV_2_Eb, CT_2_Eb, errHV_2_Eb, errCT_2_Eb);
    TGraphErrors *gV_DCR_Del_1_Eb  = new TGraphErrors(n_DCR, HV_1, DCR_Del_1_Eb, errHV_1, errDCR_Del_1_Eb);
    TGraphErrors *gV_DCR_Del_2_Eb  = new TGraphErrors(n_DCR, HV_2_Eb, DCR_Del_2_Eb, errHV_2_Eb, errDCR_Del_2_Eb);
    TGraphErrors *gV_CT_Del_1_Eb  = new TGraphErrors(n_DCR, HV_1, CT_Del_1_Eb, errHV_1, errCT_Del_1_Eb);
    TGraphErrors *gV_CT_Del_2_Eb  = new TGraphErrors(n_DCR, HV_2_Eb, CT_Del_2_Eb, errHV_2_Eb, errCT_Del_2_Eb);


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


    //------------------------------

    gV_CT_1->SetMarkerStyle(20);
    gV_CT_1->SetMarkerColor(kOrange+1);
    gV_CT_1->SetTitle();
    gV_CT_1->GetXaxis()->SetTitle("Bias Voltage (V)");
    gV_CT_1->GetYaxis()->SetTitle("Cross Talk");

    gV_CT_2->SetMarkerStyle(20);
    gV_CT_2->SetMarkerColor(kRed);
    gV_CT_2->SetTitle();
    gV_CT_2->GetXaxis()->SetTitle("Bias Voltage (V)");
    gV_CT_2->GetYaxis()->SetTitle("Cross Talk");


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


    //------------------------------

    gV_CT_Del_1->SetMarkerStyle(22);
    gV_CT_Del_1->SetMarkerColor(kOrange+1);
    gV_CT_Del_1->SetTitle();
    gV_CT_Del_1->GetXaxis()->SetTitle("Bias Voltage (V)");
    gV_CT_Del_1->GetYaxis()->SetTitle("Cross Talk");

    gV_CT_Del_2->SetMarkerStyle(22);
    gV_CT_Del_2->SetMarkerColor(kRed);
    gV_CT_Del_2->SetTitle();
    gV_CT_Del_2->GetXaxis()->SetTitle("Bias Voltage (V)");
    gV_CT_Del_2->GetYaxis()->SetTitle("Cross Talk");


    //------------------------------

    if(draw_all_bool){
        TCanvas *cV_DCR_1 = new TCanvas("cV_DCR_1", "cV_DCR_1",w,h);
        cV_DCR_1->SetGrid();
        gV_DCR_1->Draw("AP");

        TCanvas *cV_DCR_2 = new TCanvas("cV_DCR_2", "cV_DCR_2",w,h);
        cV_DCR_2->SetGrid();
        gV_DCR_2->Draw("AP");

    }

    //------------------------------

    if(draw_all_bool){
        TCanvas *cV_CT_1 = new TCanvas("cV_CT_1", "cV_CT_1",w,h);
        cV_CT_1->SetGrid();
        gV_CT_1->Draw("AP");

        TCanvas *cV_CT_2 = new TCanvas("cV_CT_2", "cV_CT_2",w,h);
        cV_CT_2->SetGrid();
        gV_CT_2->Draw("AP");

    }


    //------------------------------

    auto legendDCR = new TLegend(0.15,0.70,0.35,0.85);
    legendDCR->AddEntry(gV_DCR_1,"SensL MicroFJ (1)","p");
    legendDCR->AddEntry(gV_DCR_2,"SensL MicroFJ (2)","p");

    auto legendCT = new TLegend(0.15,0.70,0.35,0.85);
    legendCT->AddEntry(gV_CT_1,"SensL MicroFJ (1)","p");
    legendCT->AddEntry(gV_CT_2,"SensL MicroFJ (2)","p");



    //------------------------------


    TCanvas *cDCR = new TCanvas("cDCR", "cDCR",w,h);
    cDCR->SetGrid();
    TMultiGraph *mgDCR = new TMultiGraph("mgDCR", title_DCR_mg);
    mgDCR->Add(gV_DCR_1);
    mgDCR->Add(gV_DCR_2);
    mgDCR->Draw("AP");
    legendDCR->Draw();

    //------------------------------

    TCanvas *cCT = new TCanvas("cCT", "cCT",w,h);
    cCT->SetGrid();
    TMultiGraph *mgCT = new TMultiGraph("mgCT", ";Bias Voltage (V); P_{CT}");
    mgCT->Add(gV_CT_1);
    mgCT->Add(gV_CT_2);
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


    }


    //------------------------------

    if(draw_all_bool){
        TCanvas *cV_CT_Del_1 = new TCanvas("cV_CT_Del_1", "cV_CT_Del_1",w,h);
        cV_CT_Del_1->SetGrid();
        gV_CT_Del_1->Draw("AP");

        TCanvas *cV_CT_Del_2 = new TCanvas("cV_CT_Del_2", "cV_CT_Del_2",w,h);
        cV_CT_Del_2->SetGrid();
        gV_CT_Del_2->Draw("AP");


    }


    //------------------------------


    auto legendDCR_Del = new TLegend(0.15,0.70,0.35,0.85);
    legendDCR_Del->AddEntry(gV_DCR_Del_1,"SensL MicroFJ (1)","p");
    legendDCR_Del->AddEntry(gV_DCR_Del_2,"SensL MicroFJ (2)","p");

    auto legendCT_Del = new TLegend(0.15,0.70,0.35,0.85);
    legendCT_Del->AddEntry(gV_CT_Del_1,"SensL MicroFJ (1)","p");
    legendCT_Del->AddEntry(gV_CT_Del_2,"SensL MicroFJ (2)","p");

    //------------------------------

    auto legendDCR_CNT_Del = new TLegend(0.15,0.70,0.35,0.85);
    legendDCR_CNT_Del->AddEntry(gV_DCR_1,"","p");
    legendDCR_CNT_Del->AddEntry(gV_DCR_Del_1,"SensL MicroFJ (1)","p");
    legendDCR_CNT_Del->AddEntry(gV_DCR_2,"","p");
    legendDCR_CNT_Del->AddEntry(gV_DCR_Del_2,"SensL MicroFJ (2)","p");
    legendDCR_CNT_Del->SetNColumns(2);

    auto legendCT_CNT_Del = new TLegend(0.15,0.70,0.35,0.85);
    legendCT_CNT_Del->AddEntry(gV_CT_1,"","p");
    legendCT_CNT_Del->AddEntry(gV_CT_Del_1,"SensL MicroFJ (1)","p");
    legendCT_CNT_Del->AddEntry(gV_CT_2,"","p");
    legendCT_CNT_Del->AddEntry(gV_CT_Del_2,"SensL MicroFJ (2)","p");
    legendCT_CNT_Del->SetNColumns(2);


    //------------------------------

    TCanvas *cDCR_Del = new TCanvas("cDCR_Del", "cDCR_Del",w,h);
    cDCR_Del->SetGrid();
    TMultiGraph *mgDCR_Del = new TMultiGraph("mgDCR_Del", title_DCR_mg);
    mgDCR_Del->Add(gV_DCR_Del_1);
    mgDCR_Del->Add(gV_DCR_Del_2);
    mgDCR_Del->Draw("AP");
    legendDCR_Del->Draw();


    //------------------------------

    TCanvas *cCT_Del = new TCanvas("cCT_Del", "cCT_Del",w,h);
    cCT_Del->SetGrid();
    TMultiGraph *mgCT_Del = new TMultiGraph("mgCT_Del", ";Bias Voltage (V); P_{CT}");
    mgCT_Del->Add(gV_CT_Del_1);
    mgCT_Del->Add(gV_CT_Del_2);
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




    }


    //------------------------------
    //------------------------------

    TCanvas *cDCR_CNT_Del = new TCanvas("cDCR_CNT_Del", "cDCR_CNT_Del",w,h);
    cDCR_CNT_Del->SetGrid();
    TMultiGraph *mgDCR_CNT_Del = new TMultiGraph("mgDCR_CNT_Del", title_DCR_mg);
    mgDCR_CNT_Del->Add(gV_DCR_1);
    mgDCR_CNT_Del->Add(gV_DCR_2);
    mgDCR_CNT_Del->Add(gV_DCR_Del_1);
    mgDCR_CNT_Del->Add(gV_DCR_Del_2);
    mgDCR_CNT_Del->Draw("AP");
    legendDCR_CNT_Del->Draw();


    //------------------------------

    TCanvas *cCT_CNT_Del = new TCanvas("cCT_CNT_Del", "cCT_CNT_Del",w,h);
    cCT_CNT_Del->SetGrid();
    TMultiGraph *mgCT_CNT_Del = new TMultiGraph("mgCT_CNT_Del", ";Bias Voltage (V); P_{CT}");
    mgCT_CNT_Del->Add(gV_CT_1);
    mgCT_CNT_Del->Add(gV_CT_2);
    mgCT_CNT_Del->Add(gV_CT_Del_1);
    mgCT_CNT_Del->Add(gV_CT_Del_2);
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
    int n_SiPM_tot = 2;

    n_SiPM = 0;
    for(int i=0; i<N_Z; i++){
        testZ_DCR[n_SiPM][i] = TMath::Abs((DCR_1[i] - DCR_Del_1[i]) / (TMath::Sqrt( errDCR_1[i]*errDCR_1[i] + errDCR_Del_1[i]*errDCR_Del_1[i] )));
    }

    n_SiPM = 1;
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
    cout<<endl<<endl;
    cout<<"For LaTeX"<<endl<<endl;
    cout<<"\% From file DCR_CrossTalk_SensL_MicroFJ_SMTPA_60035_DARK_data_2018_07.C"<<endl;
    cout<<"\% From cnt"<<endl;
    cout<<"\% HV & DCR_1 (kHz/mm^2) & DCR_2 (kHz/mm^2) & CT_1  & CT_2   \\\\"<<endl;
    for(int i=0; i<n_DCR; i++){
        printf("$ %.2f \\pm %.2f $ & $ %.2f \\pm %.2f $ & $ %.2f \\pm %.2f $ & $ %.2f \\pm %.2f $ & $ %.2f \\pm %.2f $  \\\\ \n", HV_1[i], errHV_1[i], DCR_1[i], errDCR_1[i],DCR_2[i], errDCR_2[i], CT_1[i], errCT_1[i],CT_2[i], errCT_2[i]);
    }


    //----------

    gV_DCR_1_Ea->SetMarkerStyle(4);
    gV_DCR_2_Ea->SetMarkerStyle(4);
    gV_CT_1_Ea->SetMarkerStyle(4);
    gV_CT_2_Ea->SetMarkerStyle(4);
    gV_DCR_Del_1_Ea->SetMarkerStyle(4);
    gV_DCR_Del_2_Ea->SetMarkerStyle(4);
    gV_CT_Del_1_Ea->SetMarkerStyle(4);
    gV_CT_Del_2_Ea->SetMarkerStyle(4);
    gV_DCR_1_Eb->SetMarkerStyle(4);
    gV_DCR_2_Eb->SetMarkerStyle(4);
    gV_CT_1_Eb->SetMarkerStyle(4);
    gV_CT_2_Eb->SetMarkerStyle(4);
    gV_DCR_Del_1_Eb->SetMarkerStyle(4);
    gV_DCR_Del_2_Eb->SetMarkerStyle(4);
    gV_CT_Del_1_Eb->SetMarkerStyle(4);
    gV_CT_Del_2_Eb->SetMarkerStyle(4);


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

    }


    ////////////////////////////////////////////////////////////////////////////
    //  PRINT
    ////////////////////////////////////////////////////////////////////////////
    cout<<endl<<endl;
    cout<<"//-------------"<<endl;
    cout<<"// SensL - DARK"<<endl;
    cout<<"//-------------"<<endl;
    cout<<"double HV_SensL_DARK[] = {"; for(int i=0; i<n_DCR-1; i++) cout<<HV_1[i]<<", "; cout<<HV_1[n_DCR-1]<<"};"<<endl;
    cout<<"double errHV_SensL_DARK[] = {"; for(int i=0; i<n_DCR-1; i++) cout<<errHV_1[i]<<", "; cout<<errHV_1[n_DCR-1]<<"};"<<endl;

    cout<<"double DCR_1_SensL[] = {"; for(int i=0; i<n_DCR-1; i++) cout<<DCR_1[i]<<", "; cout<<DCR_1[n_DCR-1]<<"};"<<endl;
    cout<<"double errDCR_1_SensL[] = {"; for(int i=0; i<n_DCR-1; i++) cout<<errDCR_1[i]<<", "; cout<<errDCR_1[n_DCR-1]<<"};"<<endl;
    cout<<"double DCR_Del_1_SensL[] = {"; for(int i=0; i<n_DCR-1; i++) cout<<DCR_Del_1[i]<<", "; cout<<DCR_Del_1[n_DCR-1]<<"};"<<endl;
    cout<<"double errDCR_Del_1_SensL[] = {"; for(int i=0; i<n_DCR-1; i++) cout<<errDCR_Del_1[i]<<", "; cout<<errDCR_Del_1[n_DCR-1]<<"};"<<endl;

    cout<<"double DCR_2_SensL[] = {"; for(int i=0; i<n_DCR-1; i++) cout<<DCR_2[i]<<", "; cout<<DCR_2[n_DCR-1]<<"};"<<endl;
    cout<<"double errDCR_2_SensL[] = {"; for(int i=0; i<n_DCR-1; i++) cout<<errDCR_2[i]<<", "; cout<<errDCR_2[n_DCR-1]<<"};"<<endl;
    cout<<"double DCR_Del_2_SensL[] = {"; for(int i=0; i<n_DCR-1; i++) cout<<DCR_Del_2[i]<<", "; cout<<DCR_Del_2[n_DCR-1]<<"};"<<endl;
    cout<<"double errDCR_Del_2_SensL[] = {"; for(int i=0; i<n_DCR-1; i++) cout<<errDCR_Del_2[i]<<", "; cout<<errDCR_Del_2[n_DCR-1]<<"};"<<endl;

    cout<<"double CT_1_SensL[] = {"; for(int i=0; i<n_DCR-1; i++) cout<<CT_1[i]<<", "; cout<<CT_1[n_DCR-1]<<"};"<<endl;
    cout<<"double errCT_1_SensL[] = {"; for(int i=0; i<n_DCR-1; i++) cout<<errCT_1[i]<<", "; cout<<errCT_1[n_DCR-1]<<"};"<<endl;
    cout<<"double CT_Del_1_SensL[] = {"; for(int i=0; i<n_DCR-1; i++) cout<<CT_Del_1[i]<<", "; cout<<CT_Del_1[n_DCR-1]<<"};"<<endl;
    cout<<"double errCT_Del_1_SensL[] = {"; for(int i=0; i<n_DCR-1; i++) cout<<errCT_Del_1[i]<<", "; cout<<errCT_Del_1[n_DCR-1]<<"};"<<endl;

    cout<<"double CT_2_SensL[] = {"; for(int i=0; i<n_DCR-1; i++) cout<<CT_2[i]<<", "; cout<<CT_2[n_DCR-1]<<"};"<<endl;
    cout<<"double errCT_2_SensL[] = {"; for(int i=0; i<n_DCR-1; i++) cout<<errCT_2[i]<<", "; cout<<errCT_2[n_DCR-1]<<"};"<<endl;
    cout<<"double CT_Del_2_SensL[] = {"; for(int i=0; i<n_DCR-1; i++) cout<<CT_Del_2[i]<<", "; cout<<CT_Del_2[n_DCR-1]<<"};"<<endl;
    cout<<"double errCT_Del_2_SensL[] = {"; for(int i=0; i<n_DCR-1; i++) cout<<errCT_Del_2[i]<<", "; cout<<errCT_Del_2[n_DCR-1]<<"};"<<endl;



    ////////////////////////////////////////////////////////////////////////////
    //          OLD
    ////////////////////////////////////////////////////////////////////////////

    // PERCENTAGE ERROR
    // if(percentage_error_bool){
    //     double err_rel = 0.05;
    //     for(int i=0; i<n_DCR; i++){
    //         errDCR_1[i] = err_rel * DCR_1[i];
    //     }
    //     for(int i=0; i<n_DCR; i++){
    //         errDCR_2[i] = err_rel * DCR_2[i];
    //     }
    //
    //     for(int i=0; i<n_DCR; i++){
    //         errCT_1[i] = err_rel * CT_1[i];
    //     }
    //     for(int i=0; i<n_DCR; i++){
    //         errCT_2[i] = err_rel * CT_2[i];
    //     }
    //
    //     for(int i=0; i<n_DCR; i++){
    //         errDCR_Del_1[i] = err_rel * DCR_Del_1[i];
    //     }
    //     for(int i=0; i<n_DCR; i++){
    //         errDCR_Del_2[i] = err_rel * DCR_Del_2[i];
    //     }
    //
    //     for(int i=0; i<n_DCR; i++){
    //         errCT_Del_1[i] = err_rel * CT_Del_1[i];
    //     }
    //     for(int i=0; i<n_DCR; i++){
    //         errCT_Del_2[i] = err_rel * CT_Del_2[i];
    //     }
    //
    // }

    //------------------------------

    // FIXED ERROR
    // if(fix_error_bool){
    //     double err_fix_DCR = 1;
    //     double err_fix_CT = 0.02;
    //     for(int i=0; i<n_DCR; i++){
    //         errDCR_1[i] = err_fix_DCR;
    //     }
    //     for(int i=0; i<n_DCR; i++){
    //         errDCR_2[i] = err_fix_DCR;
    //     }
    //
    //     for(int i=0; i<n_DCR; i++){
    //         errCT_1[i] = err_fix_CT;
    //     }
    //     for(int i=0; i<n_DCR; i++){
    //         errCT_2[i] = err_fix_CT;
    //     }
    //
    //     for(int i=0; i<n_DCR; i++){
    //         errDCR_Del_1[i] = err_fix_DCR;
    //     }
    //     for(int i=0; i<n_DCR; i++){
    //         errDCR_Del_2[i] = err_fix_DCR;
    //     }
    //
    //     for(int i=0; i<n_DCR; i++){
    //         errCT_Del_1[i] = err_fix_CT;
    //     }
    //     for(int i=0; i<n_DCR; i++){
    //         errCT_Del_2[i] = err_fix_CT;
    //     }
    //
    // }


    //------------------------------




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


double GetMaxCheckPercentage(double num, double LimUp, double LimDown, double min_percentage, double max_percentage){
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
