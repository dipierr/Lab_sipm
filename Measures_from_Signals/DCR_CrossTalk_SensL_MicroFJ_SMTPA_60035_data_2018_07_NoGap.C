/******************************************************************************\
 *  DCR_CrossTalk_SensL_MicroFJ_SMTPA_60035_data_2018_07_NoGap.C
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
 *  > No gap between peaks (when I found a peak I don't jump of a certain value
 *    but I only do i++, i.e. gap_between_peaks=1)
 *  > trace_time = SUM(trace_DLED[0][trace_DLED_length-1] - trace_DLED[0][0])/
 *                 / n_ev_tot
 *
 *  > minyhistDelays = 50;  maxyhistDelays = 250; bins_Delays = 100;
 *    expDelLow_max  = minyhistDelays; expDelHigh_max = maxyhistDelays;
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

#define n_DCR_1 4
#define n_DCR_2 4


#define h 600
#define w 1000

void DCR_CrossTalk_SensL_MicroFJ_SMTPA_60035_data_2018_07_NoGap();
int find_index(double v[],int N, double value);

char title_DCR[80];
char title_DCR_mg[80];

void DCR_CrossTalk_SensL_MicroFJ_SMTPA_60035_data_2018_07_NoGap(){

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
    Area = 6.07*6.07;

    // Initialization
    for(int i=0; i<n_DCR_1; i++){
        HV_1[i] = errHV_1[i] = DCR_1[i] = errDCR_1[i] = CT_1[i] = errCT_1[i] = DCR_Del_1[i] = errDCR_Del_1[i] = CT_Del_1[i] = errCT_Del_1[i] = 0.;
    }
    for(int i=0; i<n_DCR_2; i++){
        HV_2[i] = errHV_2[i] = DCR_2[i] = errDCR_2[i] = CT_2[i] = errCT_2[i] = DCR_Del_2[i] = errDCR_Del_2[i] = CT_Del_2[i] = errCT_Del_2[i] = 0.;
    }

    // HV
    HV_1[0]    = 28.00;
    errHV_1[0] =  0.01;
    for(int i=1; i<n_DCR_1; i++){
        HV_1[i]    = HV_1[i-1]+1.;
        errHV_1[i] = errHV_1[0];
    }

    HV_2[0]    = 28.00;
    errHV_2[0] =  0.01;
    for(int i=1; i<n_DCR_2; i++){
        HV_2[i]    = HV_2[i-1]+1.;
        errHV_2[i] = errHV_2[0];
    }




    ///////////////////////////////////////////////////////////////////////////
    //      SiPM1
    ///////////////////////////////////////////////////////////////////////////
    // pe_0_5_vect[4] = {8, 8, 9, 9}
    // pe_1_5_vect[4] = {17, 20, 22, 25}
    HV = 28.00;
    index = find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV);
    DCR_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 3.65626;
    errDCR_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.00602383;
    CT_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.165949;
    DCR_Del_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 4.45971;
    errDCR_Del_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.0457989;
    CT_Del_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.318728;

    HV = 29.00;
    index = find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV);
    DCR_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 5.74691;
    errDCR_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.0075522;
    CT_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.228337;
    DCR_Del_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 6.9172;
    errDCR_Del_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.0345151;
    CT_Del_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.309471;

    HV = 30.00;
    index = find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV);
    DCR_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 7.5311;
    errDCR_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.00864545;
    CT_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.312339;
    DCR_Del_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 8.74777;
    errDCR_Del_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.0307576;
    CT_Del_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.371864;

    HV = 31.00;
    index = find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV);
    DCR_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 9.75207;
    errDCR_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.00983805;
    CT_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.394809;
    DCR_Del_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 10.5298;
    errDCR_Del_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.0289919;
    CT_Del_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0.43609;





    ///////////////////////////////////////////////////////////////////////////
    //      SiPM2
    ///////////////////////////////////////////////////////////////////////////
    // pe_0_5_vect[4] = {10, 10, 12, 15}
    // pe_1_5_vect[4] = {22, 26, 30, 35}
    HV = 28.00;
    index = find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV);
    DCR_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 4.82297;
    errDCR_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.00691851;
    CT_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.144258;
    DCR_Del_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 5.77328;
    errDCR_Del_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.0377755;
    CT_Del_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.292451;

    HV = 29.00;
    index = find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV);
    DCR_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 6.51595;
    errDCR_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.00804167;
    CT_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.21519;
    DCR_Del_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 7.76909;
    errDCR_Del_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.0325642;
    CT_Del_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.301591;

    HV = 30.00;
    index = find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV);
    DCR_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 8.01951;
    errDCR_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.00892139;
    CT_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.295959;
    DCR_Del_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 9.51106;
    errDCR_Del_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.0300056;
    CT_Del_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.351001;

    HV = 31.00;
    index = find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV);
    DCR_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 9.91351;
    errDCR_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.00991915;
    CT_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.387812;
    DCR_Del_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 11.0824;
    errDCR_Del_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.02822;
    CT_Del_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0.412221;



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

        for(int i=0; i<n_DCR_1; i++){
            errCT_1[i] = err_rel * CT_1[i];
        }
        for(int i=0; i<n_DCR_2; i++){
            errCT_2[i] = err_rel * CT_2[i];
        }

        for(int i=0; i<n_DCR_1; i++){
            errDCR_Del_1[i] = err_rel * DCR_Del_1[i];
        }
        for(int i=0; i<n_DCR_1; i++){
            errDCR_Del_2[i] = err_rel * DCR_Del_2[i];
        }

        for(int i=0; i<n_DCR_1; i++){
            errCT_Del_1[i] = err_rel * CT_Del_1[i];
        }
        for(int i=0; i<n_DCR_1; i++){
            errCT_Del_2[i] = err_rel * CT_Del_2[i];
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

        for(int i=0; i<n_DCR_1; i++){
            errCT_1[i] = err_fix_CT;
        }
        for(int i=0; i<n_DCR_2; i++){
            errCT_2[i] = err_fix_CT;
        }

        for(int i=0; i<n_DCR_1; i++){
            errDCR_Del_1[i] = err_fix_DCR;
        }
        for(int i=0; i<n_DCR_1; i++){
            errDCR_Del_2[i] = err_fix_DCR;
        }

        for(int i=0; i<n_DCR_1; i++){
            errCT_Del_1[i] = err_fix_CT;
        }
        for(int i=0; i<n_DCR_1; i++){
            errCT_Del_2[i] = err_fix_CT;
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


    }


    //------------------------------
    //------------------------------


    TGraphErrors *gV_DCR_1  = new TGraphErrors(n_DCR_1, HV_1, DCR_1, errHV_1, errDCR_1);
    TGraphErrors *gV_DCR_2  = new TGraphErrors(n_DCR_2, HV_2, DCR_2, errHV_2, errDCR_2);

    TGraphErrors *gV_CT_1  = new TGraphErrors(n_DCR_1, HV_1, CT_1, errHV_1, errCT_1);
    TGraphErrors *gV_CT_2  = new TGraphErrors(n_DCR_2, HV_2, CT_2, errHV_2, errCT_2);

    TGraphErrors *gV_DCR_Del_1  = new TGraphErrors(n_DCR_1, HV_1, DCR_Del_1, errHV_1, errDCR_Del_1);
    TGraphErrors *gV_DCR_Del_2  = new TGraphErrors(n_DCR_2, HV_2, DCR_Del_2, errHV_2, errDCR_Del_2);

    TGraphErrors *gV_CT_Del_1  = new TGraphErrors(n_DCR_1, HV_1, CT_Del_1, errHV_1, errCT_Del_1);
    TGraphErrors *gV_CT_Del_2  = new TGraphErrors(n_DCR_2, HV_2, CT_Del_2, errHV_2, errCT_Del_2);


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
    int N_Z = n_DCR_1;
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
    cout<<"\% From file DCR_CrossTalk_SensL_MicroFJ_SMTPA_60035_data_2018_07.C"<<endl;
    cout<<"\% From cnt"<<endl;
    cout<<"\% HV & DCR_1 (kHz/mm^2) & DCR_2 (kHz/mm^2) & CT_1  & CT_2   \\\\"<<endl;
    for(int i=0; i<n_DCR_1; i++){
        printf("$ %.2f \\pm %.2f $ & $ %.2f \\pm %.2f $ & $ %.2f \\pm %.2f $ & $ %.2f \\pm %.2f $ & $ %.2f \\pm %.2f $  \\\\ \n", HV_1[i], errHV_1[i], DCR_1[i], errDCR_1[i],DCR_2[i], errDCR_2[i], CT_1[i], errCT_1[i],CT_2[i], errCT_2[i]);
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
