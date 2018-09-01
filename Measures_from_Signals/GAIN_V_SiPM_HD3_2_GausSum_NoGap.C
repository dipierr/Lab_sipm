/******************************************************************************\
 * GAIN_V_SiPM_HD3_2_GausSum_NoGap.C
 *
 * GAIN values obtained by Ana_Traces_SiPM.cxx (version of 06/08/2018, 4)
 *
 * KEY POINTS:
 *  > Ana1(...)
 *  > dleddt = 6
 *  > NO trace smoothing
 *  > No gap between peaks (when I found a peak I don't jump of a certain value
 *    but I only do i++, i.e. gap_between_peaks=1)
 *  > trace_time = SUM(trace_DLED[0][trace_DLED_length-1] - trace_DLED[0][0])/
 *                 / n_ev_tot
 *  > thr (parameter of Ana1(...)) set depending on the situation
 *  > fitted using gaus_sum_12, i.e.:
 *      "[0]*TMath::Exp( - (x-[2]-[3])*(x-[2]-[3])/( 2*([4]*[4] + [5]*[5] )) )
 *       + [1]*TMath::Exp( - (x-[2]-2*[3])*(x-[2]-2*[3])/( 2*([4]*[4] +
 *      4*[5]*[5] ) ) ) "
 *  > fit range set manually
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
#include <vector>
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

#define n_GAIN_1 6
#define n_GAIN_2 6
#define n_GAIN_3 6


#define h 600
#define w 1000

void GAIN_V_SiPM_HD3_2_GausSum_NoGap();
int find_index(double vect[],int N, double value);

void GAIN_V_SiPM_HD3_2_GausSum_NoGap(){

    // HV:
    // SiPM1:
    double HV_1[n_GAIN_1], errHV_1[n_GAIN_1];
    double GAIN_1[n_GAIN_1], errGAIN_1[n_GAIN_1];
    // SiPM2:
    double HV_2[n_GAIN_2], errHV_2[n_GAIN_2];
    double GAIN_2[n_GAIN_2], errGAIN_2[n_GAIN_2];
    // SiPM3:
    double HV_3[n_GAIN_3], errHV_3[n_GAIN_3];
    double GAIN_3[n_GAIN_3], errGAIN_3[n_GAIN_3];

    double HV = 0.;

    // index:
    int index = 0;

    // ERRORS
    bool percentage_error_bool = false;
    bool fix_error_bool = true;

    // Initialization
    for(int i=0; i<n_GAIN_1; i++){
        HV_1[i] = errHV_1[i] = GAIN_1[i] = errGAIN_1[i] = 0.;
    }
    for(int i=0; i<n_GAIN_2; i++){
        HV_2[i] = errHV_2[i] = GAIN_2[i] = errGAIN_2[i] = 0.;
    }
    for(int i=0; i<n_GAIN_3; i++){
        HV_3[i] = errHV_3[i] = GAIN_3[i] = errGAIN_3[i] = 0.;
    }

    // HV
    HV_1[0]    = 32.00;
    errHV_1[0] =  0.01;
    for(int i=1; i<n_GAIN_1; i++){
        HV_1[i]    = HV_1[i-1]+1.;
        errHV_1[i] = errHV_1[0];
    }

    HV_2[0]    = 32.00;
    errHV_2[0] =  0.01;
    for(int i=1; i<n_GAIN_2; i++){
        HV_2[i]    = HV_2[i-1]+1.;
        errHV_2[i] = errHV_2[0];
    }

    HV_3[0]    = 32.00;
    errHV_3[0] =  0.01;
    for(int i=1; i<n_GAIN_3; i++){
        HV_3[i]    = HV_3[i-1]+1.;
        errHV_3[i] = errHV_3[0];
    }



    ///////////////////////////////////////////////////////////////////////////
    //      SiPM1
    ///////////////////////////////////////////////////////////////////////////
    HV = 32.;
    index = find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV);
    // cout<< index <<endl;
    GAIN_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 1.36493e+01; //1.36383e+01;
    errGAIN_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 1.10267e-02; //1.12055e-02;

    HV = 33;
    GAIN_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 1.64070e+01; //1.59077e+01;
    errGAIN_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 8.67319e-03; //9.03540e-03;

    HV = 34;
    GAIN_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 1.90921e+01; //1.89570e+01;
    errGAIN_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 8.31971e-03; //8.47637e-03;

    HV = 35;
    GAIN_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 2.19225e+01; //2.17825e+01;
    errGAIN_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 1.10134e-02; //8.68068e-03;

    HV = 36;
    GAIN_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 2.45017e+01; //2.45054e+01;
    errGAIN_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 8.98663e-03; //8.84089e-03;

    HV = 37;
    GAIN_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 2.71682e+01; //2.68621e+01;
    errGAIN_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 9.74418e-03; //2.91218e-02;



    ///////////////////////////////////////////////////////////////////////////
    //      SiPM2
    ///////////////////////////////////////////////////////////////////////////
    HV = 32.;
    index = find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV);
    // cout<< index <<endl;
    GAIN_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 1.26260e+01; //1.31614e+01;
    errGAIN_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 1.14298e-02; //1.19338e-02;

    HV = 33;
    GAIN_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 1.49182e+01; //1.50984e+01;
    errGAIN_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 9.93221e-03; //1.09344e-02;

    HV = 34;
    GAIN_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 1.77442e+01; // 1.74047e+01;
    errGAIN_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 9.17827e-03; //9.84946e-03;

    HV = 35;
    GAIN_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 1.99976e+01; //2.02401e+01;
    errGAIN_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 8.82832e-03; //8.91698e-03;

    HV = 36;
    GAIN_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 2.24596e+01; //2.25491e+01  ;
    errGAIN_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 8.99406e-03; //8.93731e-03;

    HV = 37;
    GAIN_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 2.46676e+01; //2.51216e+01;
    errGAIN_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 9.50977e-03; //9.41230e-03;


    ///////////////////////////////////////////////////////////////////////////
    //      SiPM3
    ///////////////////////////////////////////////////////////////////////////
    HV = 32.;
    index = find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV);
    // cout<< index <<endl;
    GAIN_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 1.32656e+01; // 1.32845e+01;
    errGAIN_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 9.13034e-03; // 1.09552e-02;

    HV = 33;
    GAIN_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 1.64862e+01; // 1.60972e+01;
    errGAIN_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 8.19336e-03; // 8.65181e-03;

    HV = 34;
    GAIN_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 1.92959e+01; // 1.92261e+01;
    errGAIN_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 7.85349e-03; //8.09490e-03 ;

    HV = 35;
    GAIN_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 2.20319e+01; // 2.20918e+01;
    errGAIN_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 8.14659e-03; // 2.59910e-02;

    HV = 36;
    GAIN_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 2.45662e+01; // 2.40640e+01;
    errGAIN_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 8.84217e-03; // 2.93357e-02;

    HV = 37;
    GAIN_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 2.72791e+01; // 2.69107e+01;
    errGAIN_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 9.65506e-03; // 3.13209e-02;


    //------------------------------


    // PERCENTAGE ERROR
    if(percentage_error_bool){
        double err_rel = 0.05;
        for(int i=0; i<n_GAIN_1; i++){
            errGAIN_1[i] = err_rel * GAIN_1[i];
        }
        for(int i=0; i<n_GAIN_2; i++){
            errGAIN_2[i] = err_rel * GAIN_2[i];
        }
        for(int i=0; i<n_GAIN_3; i++){
            errGAIN_3[i] = err_rel * GAIN_3[i];
        }
    }
    if(fix_error_bool){
        double err_fix = 1;
        for(int i=0; i<n_GAIN_1; i++){
            errGAIN_1[i] = err_fix;
        }
        for(int i=0; i<n_GAIN_2; i++){
            errGAIN_2[i] = err_fix;
        }
        for(int i=0; i<n_GAIN_3; i++){
            errGAIN_3[i] = err_fix;
        }
    }



    //------------------------------


    TGraphErrors *gV_GAIN_1  = new TGraphErrors(n_GAIN_1, HV_1, GAIN_1, errHV_1, errGAIN_1);
    TGraphErrors *gV_GAIN_2  = new TGraphErrors(n_GAIN_2, HV_2, GAIN_2, errHV_2, errGAIN_2);
    TGraphErrors *gV_GAIN_3  = new TGraphErrors(n_GAIN_3, HV_3, GAIN_3, errHV_3, errGAIN_3);


    //------------------------------

    gV_GAIN_1->SetMarkerStyle(20);
    gV_GAIN_1->SetMarkerColor(kOrange+1);
    gV_GAIN_1->SetTitle();
    gV_GAIN_1->GetXaxis()->SetTitle("Bias Voltage (V)");
    gV_GAIN_1->GetYaxis()->SetTitle("GAIN (mV)");

    gV_GAIN_2->SetMarkerStyle(20);
    gV_GAIN_2->SetMarkerColor(kRed);
    gV_GAIN_2->SetTitle();
    gV_GAIN_2->GetXaxis()->SetTitle("Bias Voltage (V)");
    gV_GAIN_2->GetYaxis()->SetTitle("GAIN (mV)");

    gV_GAIN_3->SetMarkerStyle(20);
    gV_GAIN_3->SetMarkerColor(kMagenta);
    gV_GAIN_3->SetTitle();
    gV_GAIN_3->GetXaxis()->SetTitle("Bias Voltage (V)");
    gV_GAIN_3->GetYaxis()->SetTitle("GAIN (mV)");

    //------------------------------

    TCanvas *cV_GAIN_1 = new TCanvas("cV_GAIN_1", "cV_GAIN_1",w,h);
    cV_GAIN_1->SetGrid();
    gV_GAIN_1->Draw("AP");

    TCanvas *cV_GAIN_2 = new TCanvas("cV_GAIN_2", "cV_GAIN_2",w,h);
    cV_GAIN_2->SetGrid();
    gV_GAIN_2->Draw("AP");

    TCanvas *cV_GAIN_3 = new TCanvas("cV_GAIN_3", "cV_GAIN_3",w,h);
    cV_GAIN_3->SetGrid();
    gV_GAIN_3->Draw("AP");

    //------------------------------


    double V_bd_1,V_bd_2,V_bd_3;

    TF1 *linearFit1 = new TF1("linearFit1","pol1",-100,100);
    TF1 *linearFit2 = new TF1("linearFit2","pol1",-100,100);
    TF1 *linearFit3 = new TF1("linearFit3","pol1",-100,100);


    gV_GAIN_1->Fit("linearFit1");
    V_bd_1  = - linearFit1->GetParameter(0)/linearFit1->GetParameter(1);

    gV_GAIN_2->Fit("linearFit2");
    V_bd_2  = - linearFit2->GetParameter(0)/linearFit2->GetParameter(1);

    gV_GAIN_3->Fit("linearFit3");
    V_bd_3  = - linearFit3->GetParameter(0)/linearFit3->GetParameter(1);




    cout<<"V_bd_1 = "<<V_bd_1<<endl;
    cout<<"V_bd_2 = "<<V_bd_2<<endl;
    cout<<"V_bd_3 = "<<V_bd_3<<endl;

    auto legendGAIN = new TLegend(0.15,0.70,0.35,0.85);
    legendGAIN->AddEntry(gV_GAIN_1,"HD3-2 (1)","p");
    legendGAIN->AddEntry(gV_GAIN_2,"HD3-2 (2)","p");
    legendGAIN->AddEntry(gV_GAIN_3,"HD3-2 (3)","p");


    TCanvas *cGAIN = new TCanvas("cGAIN", "cGAIN",w,h);
    cGAIN->SetGrid();
    TMultiGraph *mgGAIN = new TMultiGraph("mgGAIN", ";Bias Voltage (V);GAIN (mV)");
    mgGAIN->Add(gV_GAIN_1);
    mgGAIN->Add(gV_GAIN_2);
    mgGAIN->Add(gV_GAIN_3);
    mgGAIN->Draw("AP");
    legendGAIN->Draw();


    //------------------------------

    ///////////////////////////////////////////////////////////////////////////

    //------------------------------

    double HV_Vbd_1[n_GAIN_1+1],HV_Vbd_2[n_GAIN_2+1],HV_Vbd_3[n_GAIN_3+1],errHV_Vbd_2[n_GAIN_2+1],errHV_Vbd_3[n_GAIN_3+1],errHV_Vbd_1[n_GAIN_1+1];
    double GAIN_Vbd_1[n_GAIN_1+1],GAIN_Vbd_2[n_GAIN_2+1],GAIN_Vbd_3[n_GAIN_3+1],errGAIN_Vbd_2[n_GAIN_2+1],errGAIN_Vbd_3[n_GAIN_3+1],errGAIN_Vbd_1[n_GAIN_1+1];

    HV_Vbd_1[0]      = V_bd_1;
    errHV_Vbd_1[0]   = 0.1;
    GAIN_Vbd_1[0]    = 0;
    errGAIN_Vbd_1[0] = 0;

    HV_Vbd_2[0]      = V_bd_2;
    errHV_Vbd_2[0]   = 0.1;
    GAIN_Vbd_2[0]    = 0;
    errGAIN_Vbd_2[0] = 0;

    HV_Vbd_3[0]      = V_bd_3;
    errHV_Vbd_3[0]   = 0.1;
    GAIN_Vbd_3[0]    = 0;
    errGAIN_Vbd_3[0] = 0;

    for(int i=0; i<n_GAIN_1; i++){
        HV_Vbd_1[i+1]      = HV_1[i];
        errHV_Vbd_1[i+1]   = errHV_1[i];
        GAIN_Vbd_1[i+1]    = GAIN_1[i];
        errGAIN_Vbd_1[i+1] = errGAIN_1[i];

        HV_Vbd_2[i+1]      = HV_2[i];
        errHV_Vbd_2[i+1]   = errHV_2[i];
        GAIN_Vbd_2[i+1]    = GAIN_2[i];
        errGAIN_Vbd_2[i+1] = errGAIN_2[i];

        HV_Vbd_3[i+1]      = HV_3[i];
        errHV_Vbd_3[i+1]   = errHV_3[i];
        GAIN_Vbd_3[i+1]    = GAIN_3[i];
        errGAIN_Vbd_3[i+1] = errGAIN_3[i];
    }


    TGraphErrors *gV_GAIN_Vbd_1  = new TGraphErrors(n_GAIN_1+1, HV_Vbd_1, GAIN_Vbd_1, errHV_Vbd_1, errGAIN_Vbd_1);
    TGraphErrors *gV_GAIN_Vbd_2  = new TGraphErrors(n_GAIN_2+1, HV_Vbd_2, GAIN_Vbd_2, errHV_Vbd_2, errGAIN_Vbd_2);
    TGraphErrors *gV_GAIN_Vbd_3  = new TGraphErrors(n_GAIN_3+1, HV_Vbd_3, GAIN_Vbd_3, errHV_Vbd_3, errGAIN_Vbd_3);


    //------------------------------

    gV_GAIN_Vbd_1->SetMarkerStyle(20);
    gV_GAIN_Vbd_1->SetMarkerColor(kOrange+1);
    gV_GAIN_Vbd_1->SetTitle();
    gV_GAIN_Vbd_1->GetXaxis()->SetTitle("Bias Voltage (V)");
    gV_GAIN_Vbd_1->GetYaxis()->SetTitle("GAIN (mV)");

    gV_GAIN_Vbd_2->SetMarkerStyle(20);
    gV_GAIN_Vbd_2->SetMarkerColor(kRed);
    gV_GAIN_Vbd_2->SetTitle();
    gV_GAIN_Vbd_2->GetXaxis()->SetTitle("Bias Voltage (V)");
    gV_GAIN_Vbd_2->GetYaxis()->SetTitle("GAIN (mV)");

    gV_GAIN_Vbd_3->SetMarkerStyle(20);
    gV_GAIN_Vbd_3->SetMarkerColor(kMagenta);
    gV_GAIN_Vbd_3->SetTitle();
    gV_GAIN_Vbd_3->GetXaxis()->SetTitle("Bias Voltage (V)");
    gV_GAIN_Vbd_3->GetYaxis()->SetTitle("GAIN (mV)");

    //------------------------------

    TCanvas *cV_GAIN_Vbd_1 = new TCanvas("cV_GAIN_Vbd_1", "cV_GAIN_Vbd_1",w,h);
    cV_GAIN_Vbd_1->SetGrid();
    gV_GAIN_Vbd_1->Draw("AP");

    TCanvas *cV_GAIN_Vbd_2 = new TCanvas("cV_GAIN_Vbd_2", "cV_GAIN_Vbd_2",w,h);
    cV_GAIN_Vbd_2->SetGrid();
    gV_GAIN_Vbd_2->Draw("AP");

    TCanvas *cV_GAIN_Vbd_3 = new TCanvas("cV_GAIN_Vbd_3", "cV_GAIN_Vbd_3",w,h);
    cV_GAIN_Vbd_3->SetGrid();
    gV_GAIN_Vbd_3->Draw("AP");

    //------------------------------



    TCanvas *cGAIN_Vbd = new TCanvas("cGAIN_Vbd", "cGAIN_Vbd",w,h);
    cGAIN_Vbd->SetGrid();
    TMultiGraph *mgGAIN_Vdb = new TMultiGraph("mgGAIN_Vdb", ";Bias Voltage (V);GAIN (mV)");
    mgGAIN_Vdb->Add(gV_GAIN_Vbd_1);
    mgGAIN_Vdb->Add(gV_GAIN_Vbd_2);
    mgGAIN_Vdb->Add(gV_GAIN_Vbd_3);
    mgGAIN_Vdb->Draw("AP");


    //------------------------------
    // For LaTeX
    //------------------------------
    cout<<endl<<endl;
    cout<<"For LaTeX"<<endl<<endl;
    cout<<"\% From file GAIN_V_SiPM_HD3_2_GausSum.C"<<endl;
    cout<<"\% HV & GAIN_1 (mV) & GAIN_2 (mV) & GAIN_3 (mV) \\\\"<<endl;
    for(int i=0; i<n_GAIN_1; i++){
        printf("$ %.2f \\pm %.2f $ & $ %.2f \\pm %.2f $ & $ %.2f \\pm %.2f $ & $ %.2f \\pm %.2f $ \\\\ \n", HV_1[i], errHV_1[i], GAIN_1[i], errGAIN_1[i],GAIN_2[i], errGAIN_2[i],GAIN_3[i], errGAIN_3[i]);
    }



}


int find_index(double vect[], int N, double value){
    int index =0;


    double epsilon = 0.1;

    for(int i=0; i<N; i++){
        // cout<<vect[i]<<endl;
        if((vect[i]>value-epsilon) && (vect[i]<value+epsilon)){
            index = i;
            // break;
        }
    }

    return index;

}
