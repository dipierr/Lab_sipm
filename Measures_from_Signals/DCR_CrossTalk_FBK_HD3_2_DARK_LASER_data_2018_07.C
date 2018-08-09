/******************************************************************************\
 * DCR_CrossTalk_FBK_HD3_2_DARK_LASER_data_2018_07.C
 *
 * DCR_CrossTalk_FBK_HD3_2_from_cnt_data_2018_07.C
 *      +
 * OperationPoint_LASER_PLS_FitGausSum_01.C
 *
 *  Update: 09/08/2018
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


#define h 600
#define w 1000

TGraphErrors *gV_CT_1;
TGraphErrors *gV_CT_2;
TGraphErrors *gV_CT_3;
TGraphErrors *gV_CT_Del_1;
TGraphErrors *gV_CT_Del_2;
TGraphErrors *gV_CT_Del_3;
TGraphErrors *gV_MS;



char title_DCR[80];
char title_DCR_mg[80];

void DCR_CrossTalk_FBK_HD3_2_DARK_LASER_data_2018_07();
// void DCR_CrossTalk_FBK_HD3_2_from_cnt_data_2018_07();
int find_index(double v[],int N, double value);
// void OperationPoint_LASER_PLS_FitGausSum_01();

void DCR_CrossTalk_FBK_HD3_2_DARK_LASER_data_2018_07(){

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
    bool draw_all_bool = false;

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
    gV_CT_1->GetYaxis()->SetTitle("Cross Talk");

    gV_CT_2->SetMarkerStyle(20);
    gV_CT_2->SetMarkerColor(kRed);
    gV_CT_2->SetTitle();
    gV_CT_2->GetXaxis()->SetTitle("Bias Voltage (V)");
    gV_CT_2->GetYaxis()->SetTitle("Cross Talk");

    gV_CT_3->SetMarkerStyle(20);
    gV_CT_3->SetMarkerColor(kMagenta);
    gV_CT_3->SetTitle();
    gV_CT_3->GetXaxis()->SetTitle("Bias Voltage (V)");
    gV_CT_3->GetYaxis()->SetTitle("Cross Talk");

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
    gV_CT_Del_1->GetYaxis()->SetTitle("Cross Talk");

    gV_CT_Del_2->SetMarkerStyle(22);
    gV_CT_Del_2->SetMarkerColor(kRed);
    gV_CT_Del_2->SetTitle();
    gV_CT_Del_2->GetXaxis()->SetTitle("Bias Voltage (V)");
    gV_CT_Del_2->GetYaxis()->SetTitle("Cross Talk");

    gV_CT_Del_3->SetMarkerStyle(22);
    gV_CT_Del_3->SetMarkerColor(kMagenta);
    gV_CT_Del_3->SetTitle();
    gV_CT_Del_3->GetXaxis()->SetTitle("Bias Voltage (V)");
    gV_CT_Del_3->GetYaxis()->SetTitle("Cross Talk");

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
    TMultiGraph *mgCT = new TMultiGraph("mgCT", ";Bias Voltage (V); Cross Talk");
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
    TMultiGraph *mgCT_Del = new TMultiGraph("mgCT_Del", ";Bias Voltage (V);Cross Talk");
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
        TMultiGraph *mgCT_CNT_Del_1 = new TMultiGraph("mgCT_CNT_Del_1", ";Bias Voltage (V);Cross Talk");
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
        TMultiGraph *mgCT_CNT_Del_2 = new TMultiGraph("mgCT_CNT_Del_2", ";Bias Voltage (V);Cross Talk");
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
        TMultiGraph *mgCT_CNT_Del_3 = new TMultiGraph("mgCT_CNT_Del_3", ";Bias Voltage (V);Cross Talk");
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
    TMultiGraph *mgCT_CNT_Del = new TMultiGraph("mgCT_CNT_Del", ";Bias Voltage (V);Cross Talk");
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




//*****************************************************************************
//*****************************************************************************
//*****************************************************************************
//*****************************************************************************
//*****************************************************************************
//*****************************************************************************
//*****************************************************************************
//*****************************************************************************

// void OperationPoint_LASER_PLS_FitGausSum_01(){
// ERRORS
// bool fix_error_bool = true;

// HV
double HV_LASER[] =  {29.00,30.00,    31.00,        32.00,        33.00,        34.00,               35.00,               36.00};
double errHV_LASER[]={0.01, 0.01,     0.01,         0.01,         0.01,         0.01,                0.01,                0.01};
// double HV_LASER[] =  {29.00,30.00,    31.00,        32.00,        33.00,        34.00,        34.50,        35.00,        35.50,        36.00};
// double errHV_LASER[]={0.01, 0.01,     0.01,         0.01,         0.01,         0.01,         0.01,         0.01,         0.01,         0.01};



// PEAK 0
// H_peak_0
double H_peak_0[n_GAIN] = {0.};
double errH_peak_0[n_GAIN] = {0.};

// Sigma_peak_0
double Sigma_peak_0[n_GAIN] = {0.};
double errSigma_peak_0[n_GAIN] = {0.};

// PEAK 1
// H_peak_1
double H_peak_1[n_GAIN] = {0.};
double errH_peak_1[n_GAIN] = {0.};

// Sigma_peak_1
double Sigma_peak_1[n_GAIN] = {0.};
double errSigma_peak_1[n_GAIN] = {0.};

// Mean_peak_1
double Mean_peak_1[n_GAIN] = {0.};
double errMean_peak_1[n_GAIN] = {0.};

// PEAK 2
// Mean_peak_2
double Mean_peak_2[n_GAIN] = {0.};
double errMean_peak_2[n_GAIN] = {0.};

// GLOBAL HIST
// Mean hist global
double Mean_hg[n_MEAN] = {0.};
double errMean_hg[n_MEAN] = {0.};

// Standard_dev_hist_global
double Std_hg[n_MEAN] = {0.};
double errStd_hg[n_MEAN] = {0.};

 // Integral of the average
double Integral[n_GAIN] = {0.};
double errIntegral[n_GAIN] = {0.};

// Entries
double Entries[n_GAIN] = {0.};

// GAIN
double GAIN[n_GAIN] = {0.};
double errGAIN[n_GAIN] = {0.};

// Other variables
double GainWeighted[n_GAIN] = {0.};
double errGainWeighted[n_GAIN] = {0.};
double Mean_hist_St_dev[n_MEAN] = {0.};
double errMean_hist_St_dev[n_MEAN] = {0.};
double Integral_GAIN[n_GAIN] = {0.};
double errIntegral_GAIN[n_GAIN] = {0.};
double HV_LASER_GAIN[n_GAIN] = {0.};
double errHV_LASER_GAIN[n_GAIN] = {0.};
double Area0[n_GAIN] = {0.};
double Area1[n_GAIN] = {0.};
double errArea0[n_GAIN] = {0.};
double errArea1[n_GAIN] = {0.};
double Prob_0pe[n_GAIN] = {0.};
double errProb_0pe[n_GAIN] = {0.};
double Prob_1pe[n_GAIN] = {0.};
double errProb_1pe[n_GAIN] = {0.};
double Prob_1peS[n_GAIN] = {0.};
double errProb_1peS[n_GAIN] = {0.};
double Mu[n_GAIN] = {0.};
double errMu[n_GAIN] = {0.};
double Prob_Cross_Talk[n_GAIN] = {0.};
double errProb_Cross_Talk[n_GAIN] = {0.};
double errEntries[n_GAIN] = {0.};



// HV_LASER = 29 V


// HV_LASER = 30 V



// HV_LASER = 31 V [0]
// Window for LED peak: (168, 177) ns  (minLED_amp, maxLED_amp)
// Other: dleddt = 9; smooth_trace_bool = 0;
//        histLED_low = -50; histLED_high = 700; histLED_binw = 0.7;
//        fit_low = -10; fit_high = 29;
H_peak_0[0]          = 798.911;
errH_peak_0[0]       = 12.4535;
Sigma_peak_0[0]      = 2.45416;
errSigma_peak_0[0]   = 0.0310133;
H_peak_1[0]          = 1180.98;
errH_peak_1[0]       = 16.6779;
Sigma_peak_1[0]      = 2.96697;
errSigma_peak_1[0]   = 0.0387114;
Mean_peak_1[0]       = 5.9277;
errMean_peak_1[0]    = 0.0543554;
Mean_peak_2[0]       = 8.38187;
errMean_peak_2[0]    = 0.0625806;
GAIN[0]              = 2.45416;
errGAIN[0]           = 0.0310133;
Mean_hg[2]           = 23.9805;
errMean_hg[2]        = 0.0713393;
Std_hg[2]            = 15.9518;
errStd_hg[2]         = 0.0504445;
Entries[0]           = 49999;


// HV_LASER = 32 V [1]
// Window for LED peak: (168, 177) ns  (minLED_amp, maxLED_amp)
// Other: dleddt = 9; smooth_trace_bool = 0;
//        histLED_low = -50; histLED_high = 700; histLED_binw = 0.7;
//        fit_low = -10; fit_high = 34;
H_peak_0[1]          = 602.069;
errH_peak_0[1]       = 9.94509;
Sigma_peak_0[1]      = 2.80037;
errSigma_peak_0[1]   = 0.0322954;
H_peak_1[1]          = 955.733;
errH_peak_1[1]       = 11.8414;
Sigma_peak_1[1]      = 3.274;
errSigma_peak_1[1]   = 0.0365093;
Mean_peak_1[1]       = 5.91894;
errMean_peak_1[1]    = 0.0551977;
Mean_peak_2[1]       = 8.71931;
errMean_peak_2[1]    = 0.0639514;
GAIN[1]              = 2.80037;
errGAIN[1]           = 0.0322954;
Mean_hg[3]           = 32.2026;
errMean_hg[3]        = 0.0987055;
Std_hg[3]            = 22.071;
errStd_hg[3]         = 0.0697953;
Entries[1]           = 49999;



// HV_LASER = 33 V [2]
// Window for LED peak: (168, 177) ns  (minLED_amp, maxLED_amp)
// Other: dleddt = 9; smooth_trace_bool = 0;
//        histLED_low = -50; histLED_high = 700; histLED_binw = 0.7;
//        fit_low = -10; fit_high = 38;
H_peak_0[2]          = 465.179;
errH_peak_0[2]       = 8.49261;
Sigma_peak_0[2]      = 3.20809;
errSigma_peak_0[2]   = 0.0409653;
H_peak_1[2]          = 782.935;
errH_peak_1[2]       = 9.37982;
Sigma_peak_1[2]      = 3.66242;
errSigma_peak_1[2]   = 0.0447216;
Mean_peak_1[2]       = 5.98651;
errMean_peak_1[2]    = 0.065181;
Mean_peak_2[2]       = 9.19459;
errMean_peak_2[2]    = 0.0769852;
GAIN[2]              = 3.20809;
errGAIN[2]           = 0.0409653;
Mean_hg[4]           = 41.6162;
errMean_hg[4]        = 0.131494;
Std_hg[4]            = 29.4028;
errStd_hg[4]         = 0.0929806;
Entries[2]           = 49999;





// HV_LASER = 34 V [3]
// Window for LED peak: (168, 177) ns  (minLED_amp, maxLED_amp)
// Other: dleddt = 9; smooth_trace_bool = 0;
//        histLED_low = -50; histLED_high = 700; histLED_binw = 0.7;
//        fit_low = -10; fit_high = 46;
H_peak_0[3]          = 367.064;
errH_peak_0[3]       = 7.19149;
Sigma_peak_0[3]      = 3.55882;
errSigma_peak_0[3]   = 0.0455046;
H_peak_1[3]          = 619.135;
errH_peak_1[3]       = 7.96463;
Sigma_peak_1[3]      = 4.17185;
errSigma_peak_1[3]   = 0.0486095;
Mean_peak_1[3]       = 5.51912;
errMean_peak_1[3]    = 0.0708883;
Mean_peak_2[3]       = 9.07794;
errMean_peak_2[3]    = 0.0842367;
GAIN[3]              = 3.55882;
errGAIN[3]           = 0.0455046;
Mean_hg[5]           = 52.3904;
errMean_hg[5]        = 0.167837;
Std_hg[5]            = 37.529;
errStd_hg[5]         = 0.118678;
Entries[3]           = 49999;




// HV_LASER = 35.00 V [4]
// Window for LED peak: (168, 177) ns  (minLED_amp, maxLED_amp)
// Other: dleddt = 9; smooth_trace_bool = 0;
//        histLED_low = -50; histLED_high = 700; histLED_binw = 0.7;
//        fit_low = -10; fit_high = 51;
H_peak_0[4]          = 310.268;
errH_peak_0[4]       = 6.52916;
Sigma_peak_0[4]      = 3.78661;
errSigma_peak_0[4]   = 0.054867;
H_peak_1[4]          = 517.979;
errH_peak_1[4]       = 6.82964;
Sigma_peak_1[4]      = 4.54551;
errSigma_peak_1[4]   = 0.0561603;
Mean_peak_1[4]       = 4.9444;
errMean_peak_1[4]    = 0.0799529;
Mean_peak_2[4]       = 8.73101;
errMean_peak_2[4]    = 0.0969682;
GAIN[4]              = 3.78661;
errGAIN[4]           = 0.054867;
Mean_hg[6]           = 63.8248;
errMean_hg[6]        = 0.20694;
Std_hg[6]            = 46.2728;
errStd_hg[6]         = 0.146329;
Entries[4]           = 49999;







// HV_LASER = 36 V ([5])
// Window for LED peak: (168, 177) ns  (minLED_amp, maxLED_amp)
// Other: dleddt = 9; smooth_trace_bool = 0;
//        histLED_low = -50; histLED_high = 700; histLED_binw = 0.7;
//        fit_low = -10; fit_high = 60;
H_peak_0[5]          = 253.065;
errH_peak_0[5]       = 5.72925;
Sigma_peak_0[5]      = 4.17716;
errSigma_peak_0[5]   = 0.0685689;
H_peak_1[5]          = 394.896;
errH_peak_1[5]       = 5.7841;
Sigma_peak_1[5]      = 5.1988;
errSigma_peak_1[5]   = 0.0663333;
Mean_peak_1[5]       = 4.39273;
errMean_peak_1[5]    = 0.0982355;
Mean_peak_2[5]       = 8.56989;
errMean_peak_2[5]    = 0.119799;
GAIN[5]              = 4.17716;
errGAIN[5]           = 0.0685689;
Mean_hg[7]           = 76.9744;
errMean_hg[7]        = 0.252346;
Std_hg[7]            = 56.4257;
errStd_hg[7]         = 0.178435;
Entries[5]           = 49999;







//------------------------------

// select only points for n_GAIN
for (int i = 0; i < n_GAIN; i++) {
    HV_LASER_GAIN[i] = HV_LASER[i+n_MEAN-n_GAIN];
    errHV_LASER_GAIN[i] = errHV_LASER[i+n_MEAN-n_GAIN];
}

for (int i=0; i<n_GAIN; i++){

    // GainWeighted
    GainWeighted[i] = GAIN[i]/TMath::Sqrt(Sigma_peak_1[i]*Sigma_peak_1[i]-Sigma_peak_0[i]*Sigma_peak_0[i]);
    errGainWeighted[i] = errGAIN[i]*errGAIN[i]/(Sigma_peak_1[i]*Sigma_peak_1[i]-Sigma_peak_0[i]*Sigma_peak_0[i]);
    errGainWeighted[i] += errSigma_peak_1[i]*errSigma_peak_1[i]*Sigma_peak_1[i]*Sigma_peak_1[i]*GAIN[i]*GAIN[i]/TMath::Power(Sigma_peak_1[i]*Sigma_peak_1[i]-Sigma_peak_0[i]*Sigma_peak_0[i],3);
    errGainWeighted[i] += errSigma_peak_0[i]*errSigma_peak_0[i]*Sigma_peak_0[i]*Sigma_peak_0[i]*GAIN[i]*GAIN[i]/TMath::Power(Sigma_peak_1[i]*Sigma_peak_1[i]-Sigma_peak_0[i]*Sigma_peak_0[i],3);
    errGainWeighted[i] = TMath::Sqrt(errGainWeighted[i]);

    // Integral_GAIN
    // Integral_GAIN[i]=TMath::Abs(Integral[i])/(GAIN[i]*Cross_Talk[i]);
    // errIntegral_GAIN[i]=Integral[i]/(GAIN[i]*Cross_Talk[i])*TMath::Power(errGAIN[i]*errGAIN[i]/(GAIN[i]*GAIN[i])+errCross_Talk[i]*errCross_Talk[i]/(Cross_Talk[i]*Cross_Talk[i]),0.5);
    //Integral_GAIN[i]=TMath::Abs(Integral[i])/(GAIN[i]);
    //errIntegral_GAIN[i]=Integral[i]*errGAIN[i]/(GAIN[i]*GAIN[i]);

    //Probabilità CrossTalk LED
    errEntries[i]=TMath::Power(Entries[i],0.5);
    Area0[i]=H_peak_0[i]*Sigma_peak_0[i]*TMath::Power(2*TMath::Pi(),0.5)/1;
    Prob_0pe[i]=Area0[i]/Entries[i];
    Mu[i]=-TMath::Log(Prob_0pe[i]);
    Prob_1pe[i]=Mu[i]*TMath::Exp(-Mu[i]);
    Area1[i]=H_peak_1[i]*Sigma_peak_1[i]*TMath::Power(2*TMath::Pi(),0.5)/1;
    Prob_1peS[i]=Area1[i]/Entries[i];
    Prob_Cross_Talk[i]=1-(Prob_1peS[i]/Prob_1pe[i]);
    cout << "Cross Talk " << Prob_Cross_Talk[i] << endl;
    cout << "Prob_1peS[i]/Prob_1pe[i]\t" << Prob_1peS[i]/Prob_1pe[i] << endl;
    cout << "p 0 \t" << Prob_0pe[i] << endl;
    cout << "p 1 s\t" << Prob_1peS[i] << endl;
    cout << "p 1\t" << Prob_1pe[i] << endl<<endl    ;


    errArea1[i]=Area1[i]*TMath::Power(TMath::Power(errH_peak_1[i]/H_peak_1[i],2)+TMath::Power(errSigma_peak_1[i]/Sigma_peak_1[i],2),0.5);
    errProb_1peS[i]=errArea1[i]/Entries[i];
    errMu[i]=errArea0[i]/(Prob_0pe[i]*Entries[i]);
    errProb_1pe[i]=TMath::Exp(-Mu[i])*TMath::Abs(1-Mu[i])*errMu[i];
    errArea0[i]=Area0[i]*TMath::Power(TMath::Power(errH_peak_0[i]/H_peak_0[i],2)+TMath::Power(errSigma_peak_0[i]/Sigma_peak_0[i],2),0.5);
    errProb_0pe[i]=errArea0[i]/Entries[i];
    errProb_Cross_Talk[i]=TMath::Power(TMath::Power(errProb_1peS[i]/Prob_1pe[i],2)+TMath::Power(Prob_1peS[i]*errProb_1pe[i]/Prob_1pe[i],2),0.5);
}
double temp;
for (int i = 0; i < n_MEAN; i++) {
  // Mean_hg/St_dev
  Mean_hist_St_dev[i]= Mean_hg[i]/Std_hg[i];


  errMean_hist_St_dev[i] = TMath::Power(Std_hg[i],-1)*TMath::Power(errMean_hg[i]*errMean_hg[i]+Mean_hg[i]*Mean_hg[i]*errStd_hg[i]*errStd_hg[i]/(Std_hg[i]*Std_hg[i]),0.5);

  // temp = errMean_hg[i]*errMean_hg[i]+Mean_hg[i]*Mean_hg[i]*errStd_hg[i]*errStd_hg[i]/(Std_hg[i]*Std_hg[i]);
  // temp = TMath::Power(temp, 0.5);
  // temp *= 1/Std_hg[i];
  //


}

///////////////////////////////////////////////////////////////////////////
//          FIXED ERROR
///////////////////////////////////////////////////////////////////////////

if(fix_error_bool){
    // double err_fix_CT = 0.015;
    double err_fix_CT = 0.02;
    for(int i=0; i<n_GAIN; i++){
        errProb_Cross_Talk[i] = err_fix_CT;
    }
}

///////////////////////////////////////////////////////////////////////////

for(int i=0; i<n_GAIN; i++){
        cout<<"Error CT = "<<errProb_Cross_Talk[i]<<endl;
    }


//------------------------------
// Gain_Weighted
//------------------------------

TGraphErrors *gV_GW  = new TGraphErrors(n_GAIN, HV_LASER_GAIN, GainWeighted, errHV_LASER_GAIN, errGainWeighted);


//------------------------------

gV_GW->SetMarkerStyle(47);
gV_GW->SetMarkerColor(kOrange+2);
gV_GW->SetTitle();
gV_GW->GetXaxis()->SetTitle("Voltage (V)");
gV_GW->GetYaxis()->SetTitle("Weighted Gain");

//------------------------------

TCanvas *cV_GW = new TCanvas("cV_GW", "cV_GW",w,h);
cV_GW->SetGrid();
gV_GW->Draw("AP");

 //------------------------------
 // Mean/Dev_st
 //------------------------------

TGraphErrors *gV_MS  = new TGraphErrors(n_MEAN, HV_LASER, Mean_hist_St_dev, errHV_LASER, errMean_hist_St_dev);


//------------------------------

gV_MS->SetMarkerStyle(47);
gV_MS->SetMarkerColor(kOrange+2);
gV_MS->SetTitle();
gV_MS->GetXaxis()->SetTitle("Voltage (V)");
gV_MS->GetYaxis()->SetTitle("Mean/St_Dev ()");

//------------------------------

TCanvas *cV_MS = new TCanvas("cV_MS", "cV_MS",w,h);
cV_MS->SetGrid();
gV_MS->Draw("AP");


//------------------------------
// Cross Talk
//------------------------------
//for(int i=0; i<n_GAIN; i++){
  //  cout<<HV_LASER_GAIN[i]<<"\t"<<Prob_Cross_Talk[i]<<endl;
//}
TGraphErrors *gV_PCT  = new TGraphErrors(n_GAIN, HV_LASER_GAIN, Prob_Cross_Talk, errHV_LASER_GAIN, errProb_Cross_Talk);


//------------------------------

gV_PCT->SetMarkerStyle(21);
// gV_PCT->SetMarkerSize(2);
gV_PCT->SetMarkerColor(kOrange+2);
gV_PCT->SetTitle();
gV_PCT->GetXaxis()->SetTitle("Voltage (V)");
gV_PCT->GetYaxis()->SetTitle("Cross Talk");

//------------------------------

TCanvas *cV_PCT = new TCanvas("cV_PCT", "cV_PCT",w,h);
cV_PCT->SetGrid();
gV_PCT->Draw("AP");






//*****************************************************************************
//*****************************************************************************
//*****************************************************************************
//*****************************************************************************
//*****************************************************************************
//*****************************************************************************
//*****************************************************************************
//*****************************************************************************



  /////////////////////////////////////////////////////////////////////////////
  //    MERGE
  /////////////////////////////////////////////////////////////////////////////

      TCanvas *cV_CT_DARK_LASER = new TCanvas("cV_CT_DARK_LASER", "cV_CT_DARK_LASER",w,h);
      cV_CT_DARK_LASER->SetGrid();
      TMultiGraph *mgCT_DARK_LASER = new TMultiGraph("mgCT_DARK_LASER", ";Bias Voltage (V);Cross Talk");
      mgCT_DARK_LASER->Add(gV_CT_1);
      // mgCT_DARK_LASER->Add(gV_CT_2);
      // mgCT_DARK_LASER->Add(gV_CT_3);
      mgCT_DARK_LASER->Add(gV_CT_Del_1);
      // mgCT_DARK_LASER->Add(gV_CT_Del_2);
      // mgCT_DARK_LASER->Add(gV_CT_Del_3);
      mgCT_DARK_LASER->Add(gV_PCT);
      mgCT_DARK_LASER->Draw("AP");

      auto legendCT_DARK_LASER = new TLegend(0.15,0.70,0.35,0.85);
      legendCT_DARK_LASER->AddEntry(gV_CT_1,    "HD3-2 (1), from DCR","p");
      legendCT_DARK_LASER->AddEntry(gV_CT_Del_1,"HD3-2 (1), from DELAYS","p");
      legendCT_DARK_LASER->AddEntry(gV_PCT,     "HD3-2 (1), from LASER","p");
      legendCT_DARK_LASER->Draw();
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
