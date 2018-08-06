/******************************************************************************\
 * DCR_CrossTalk_FBK_HD3_2_from_cnt_data_2018_07.cxx
 *
 * GAIN values obtained by Ana_Traces_SiPM.cxx (version of 07/08/2018, 1)
 *
 * KEY POINTS:
 *  > DCR_CT_1SiPM_nHVs(string filelist, int nfile_in_list, int last_event_n)
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
#define w 800

void DCR_CrossTalk_FBK_HD3_2_from_cnt_data_2018_07();
int find_index(double v[],int N, double value);

void DCR_CrossTalk_FBK_HD3_2_from_cnt_data_2018_07(){

    // HV:
    // SiPM1:
    double HV_1[n_DCR_1], errHV_1[n_DCR_1];
    double DCR_1[n_DCR_1], errDCR_1[n_DCR_1];
    double CT_1[n_DCR_1], errCT_1[n_DCR_1];
    // SiPM2:
    double HV_2[n_DCR_2], errHV_2[n_DCR_2];
    double DCR_2[n_DCR_2], errDCR_2[n_DCR_2];
    double CT_2[n_DCR_2], errCT_2[n_DCR_2];
    // SiPM3:
    double HV_3[n_DCR_3], errHV_3[n_DCR_3];
    double DCR_3[n_DCR_3], errDCR_3[n_DCR_3];
    double CT_3[n_DCR_3], errCT_3[n_DCR_3];

    double HV = 0.;

    // index:
    int index = 0;

    // Initialization
    for(int i=0; i<n_DCR_1; i++){
        HV_1[i] = errHV_1[i] = DCR_1[i] = errDCR_1[i] = CT_1[i] = errCT_1[i] = 0.;
    }
    for(int i=0; i<n_DCR_2; i++){
        HV_2[i] = errHV_2[i] = DCR_2[i] = errDCR_2[i] = CT_2[i] = errCT_2[i] = 0.;
    }
    for(int i=0; i<n_DCR_3; i++){
        HV_3[i] = errHV_3[i] = DCR_3[i] = errDCR_3[i] = CT_3[i] = errCT_3[i] = 0.;
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
    DCR_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 14.7862;
    errDCR_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.0416529;
    CT_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.205871;

    HV = 33.00;
    index = find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV);
    DCR_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 17.8257;
    errDCR_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.0464375;
    CT_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.240943;

    HV = 34.00;
    index = find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV);
    DCR_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 20.2136;
    errDCR_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.0500308;
    CT_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.26685;

    HV = 35.00;
    index = find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV);
    DCR_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 22.5155;
    errDCR_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.0533869;
    CT_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.294154;

    HV = 36.00;
    index = find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV);
    DCR_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 24.7567;
    errDCR_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.056571;
    CT_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.320874;

    HV = 37.00;
    index = find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV);
    DCR_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 27.0207;
    errDCR_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.0597175;
    CT_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.356573;


    ///////////////////////////////////////////////////////////////////////////
    //      SiPM2
    ///////////////////////////////////////////////////////////////////////////
    HV = 32.00;
    index = find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV);
    DCR_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 15.057;
    errDCR_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.0420906;
    CT_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.188846;

    HV = 33.00;
    index = find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV);
    DCR_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 18.7635;
    errDCR_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.0478639;
    CT_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.227392;

    HV = 34.00;
    index = find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV);
    DCR_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 21.4827;
    errDCR_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.0518928;
    CT_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.250666;

    HV = 35.00;
    index = find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV);
    DCR_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 23.9344;
    errDCR_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.0554115;
    CT_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.279799;

    HV = 36.00;
    index = find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV);
    DCR_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 26.3345;
    errDCR_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.0587705;
    CT_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.300818;

    HV = 37.00;
    index = find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV);
    DCR_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 28.6496;
    errDCR_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.0619437;
    CT_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.346385;



    ///////////////////////////////////////////////////////////////////////////
    //      SiPM3
    ///////////////////////////////////////////////////////////////////////////
    HV = 32.00;
    index = find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV);
    DCR_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 17.6168;
    errDCR_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 0.046117;
    CT_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 0.211385;

    HV = 33.00;
    index = find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV);
    DCR_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 21.23;
    errDCR_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 0.0515244;
    CT_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 0.248883;

    HV = 34.00;
    index = find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV);
    DCR_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 24.2307;
    errDCR_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 0.0558303;
    CT_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 0.275;

    HV = 35.00;
    index = find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV);
    DCR_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 27.119;
    errDCR_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 0.0598526;
    CT_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 0.300411;

    HV = 36.00;
    index = find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV);
    DCR_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 29.8312;
    errDCR_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 0.0635411;
    CT_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 0.323539;

    HV = 37.00;
    index = find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV);
    DCR_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 32.4055;
    errDCR_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 0.0669758;
    CT_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 0.363071;


    //------------------------------


    // PERCENTAGE ERROR
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


    //------------------------------

    TGraphErrors *gV_DCR_1  = new TGraphErrors(n_DCR_1, HV_1, DCR_1, errHV_1, errDCR_1);
    TGraphErrors *gV_DCR_2  = new TGraphErrors(n_DCR_2, HV_2, DCR_2, errHV_2, errDCR_2);
    TGraphErrors *gV_DCR_3  = new TGraphErrors(n_DCR_3, HV_3, DCR_3, errHV_3, errDCR_3);

    TGraphErrors *gV_CT_1  = new TGraphErrors(n_DCR_1, HV_1, CT_1, errHV_1, errCT_1);
    TGraphErrors *gV_CT_2  = new TGraphErrors(n_DCR_2, HV_2, CT_2, errHV_2, errCT_2);
    TGraphErrors *gV_CT_3  = new TGraphErrors(n_DCR_3, HV_3, CT_3, errHV_3, errCT_3);


    //------------------------------

    gV_DCR_1->SetMarkerStyle(20);
    gV_DCR_1->SetMarkerColor(kOrange+1);
    gV_DCR_1->SetTitle();
    gV_DCR_1->GetXaxis()->SetTitle("Bias Voltage (V)");
    gV_DCR_1->GetYaxis()->SetTitle("DCR (MHz)");

    gV_DCR_2->SetMarkerStyle(20);
    gV_DCR_2->SetMarkerColor(kRed);
    gV_DCR_2->SetTitle();
    gV_DCR_2->GetXaxis()->SetTitle("Bias Voltage (V)");
    gV_DCR_2->GetYaxis()->SetTitle("DCR (MHz)");

    gV_DCR_3->SetMarkerStyle(20);
    gV_DCR_3->SetMarkerColor(kMagenta);
    gV_DCR_3->SetTitle();
    gV_DCR_3->GetXaxis()->SetTitle("Bias Voltage (V)");
    gV_DCR_3->GetYaxis()->SetTitle("DCR (MHz)");

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

    TCanvas *cV_DCR_1 = new TCanvas("cV_DCR_1", "cV_DCR_1",w,h);
    cV_DCR_1->SetGrid();
    gV_DCR_1->Draw("AP");

    TCanvas *cV_DCR_2 = new TCanvas("cV_DCR_2", "cV_DCR_2",w,h);
    cV_DCR_2->SetGrid();
    gV_DCR_2->Draw("AP");

    TCanvas *cV_DCR_3 = new TCanvas("cV_DCR_3", "cV_DCR_3",w,h);
    cV_DCR_3->SetGrid();
    gV_DCR_3->Draw("AP");

    //------------------------------

    TCanvas *cV_CT_1 = new TCanvas("cV_CT_1", "cV_CT_1",w,h);
    cV_CT_1->SetGrid();
    gV_CT_1->Draw("AP");

    TCanvas *cV_CT_2 = new TCanvas("cV_CT_2", "cV_CT_2",w,h);
    cV_CT_2->SetGrid();
    gV_CT_2->Draw("AP");

    TCanvas *cV_CT_3 = new TCanvas("cV_CT_3", "cV_CT_3",w,h);
    cV_CT_3->SetGrid();
    gV_CT_3->Draw("AP");

    //------------------------------

    TCanvas *cDCR = new TCanvas("cDCR", "cDCR",w,h);
    cDCR->SetGrid();
    TMultiGraph *mgDCR = new TMultiGraph("mgDCR", ";Bias Voltage (V);DCR (MHz)");
    mgDCR->Add(gV_DCR_1);
    mgDCR->Add(gV_DCR_2);
    mgDCR->Add(gV_DCR_3);
    mgDCR->Draw("AP");

    //------------------------------

    TCanvas *cCT = new TCanvas("cCT", "cCT",w,h);
    cCT->SetGrid();
    TMultiGraph *mgCT = new TMultiGraph("mgCT", ";Bias Voltage (V);Cross Talk");
    mgCT->Add(gV_CT_1);
    mgCT->Add(gV_CT_2);
    mgCT->Add(gV_CT_3);
    mgCT->Draw("AP");

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
