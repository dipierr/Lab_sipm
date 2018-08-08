/******************************************************************************\
 * DCR_CrossTalk_FBK_HD3_2_from_cnt_data_2018_07.C
 *
 * GAIN values obtained by Ana_Traces_SiPM.cxx (version of 07/08/2018, 1)
 *
 * KEY POINTS:
 *  > DCR_CT_1SiPM_nHVs(...)
 *  > dleddt = 6
 *  > NO trace smoothing
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

#define n_DCR_1 6
#define n_DCR_2 6
#define n_DCR_3 6


#define h 600
#define w 1000

void DCR_CrossTalk_FBK_HD3_2_from_cnt_data_2018_07();
int find_index(double v[],int N, double value);

char title_DCR[80];
char title_DCR_mg[80];

void DCR_CrossTalk_FBK_HD3_2_from_cnt_data_2018_07(){

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

    bool dcr_area = true;

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
    HV_1[0]    = 32.00;
    errHV_1[0] =  0.01;
    for(int i=1; i<n_DCR_1; i++){
        HV_1[i]    = HV_1[i-1]+1.;
        errHV_1[i] = errHV_1[0];
    }

    HV_2[0]    = 32.00;
    errHV_2[0] =  0.01;
    for(int i=1; i<n_DCR_2; i++){
        HV_2[i]    = HV_2[i-1]+1.;
        errHV_2[i] = errHV_2[0];
    }

    HV_3[0]    = 32.00;
    errHV_3[0] =  0.01;
    for(int i=1; i<n_DCR_3; i++){
        HV_3[i]    = HV_3[i-1]+1.;
        errHV_3[i] = errHV_3[0];
    }



    ///////////////////////////////////////////////////////////////////////////
    //      SiPM1
    ///////////////////////////////////////////////////////////////////////////
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

    // PERCENTAGE ERROR
    if(fix_error_bool){
        double err_fix_DCR = 1;
        double err_fix_CT = 0.015;
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
