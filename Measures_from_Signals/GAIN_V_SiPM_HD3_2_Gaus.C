/******************************************************************************\
 * GAIN_V_SiPM_HD3_2_Gaus.cxx
 *
 * GAIN values obtained by Ana_Traces_SiPM.cxx (version of 06/08/2018, 4)
 *
 * KEY POINTS:
 *  > Ana1(...)
 *  > dleddt = 6
 *  > NO trace smoothing
 *  > thr (parameter of Ana1(...)) set depending on the situation
 *  > fitted using two different gaussians using FitPanel
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
#define w 800

void GAIN_V_SiPM_HD3_2_Gaus();
int find_index(double vect[],int N, double value);

void GAIN_V_SiPM_HD3_2_Gaus(){

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
    cout<< index <<endl;
    GAIN_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 9999;
    errGAIN_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0;

    HV = 33;
    GAIN_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 9999;
    errGAIN_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0;

    HV = 34;
    GAIN_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 9999;
    errGAIN_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0;

    HV = 35;
    GAIN_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 9999;
    errGAIN_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0;

    HV = 36;
    GAIN_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 9999;
    errGAIN_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0;

    HV = 37;
    GAIN_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 9999;
    errGAIN_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 0;



    ///////////////////////////////////////////////////////////////////////////
    //      SiPM2
    ///////////////////////////////////////////////////////////////////////////
    HV = 32.;
    index = find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV);
    cout<< index <<endl;
    GAIN_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 2.63456e+01-1.43432e+01;
    errGAIN_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0;

    HV = 33;
    GAIN_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 3.15294e+01-1.72352e+01;
    errGAIN_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0;

    HV = 34;
    GAIN_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] =  3.66303e+01-1.94682e+01;
    errGAIN_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0;

    HV = 35;
    GAIN_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 4.11993e+01-2.13305e+01;
    errGAIN_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0;

    HV = 36;
    GAIN_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 4.57549e+01 - 2.33086e+01;
    errGAIN_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0;

    HV = 37;
    GAIN_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 5.01583e+01 - 2.51468e+01;
    errGAIN_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 0;


    ///////////////////////////////////////////////////////////////////////////
    //      SiPM3
    ///////////////////////////////////////////////////////////////////////////
    HV = 32.;
    index = find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV);
    cout<< index <<endl;
    GAIN_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] =  9999;
    errGAIN_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] =  0.;

    HV = 33;
    GAIN_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] =  9999;
    errGAIN_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] =  0.;

    HV = 34;
    GAIN_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] =  9999;
    errGAIN_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 0.;

    HV = 35;
    GAIN_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] =  9999;
    errGAIN_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] =  0.;

    HV = 36;
    GAIN_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] =  9999;
    errGAIN_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] =  0.;

    HV = 37;
    GAIN_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] = 9999;
    errGAIN_3[find_index(HV_3,  sizeof(HV_3)/sizeof(HV_3[0]), HV)] =  0.;



    //------------------------------


    // PERCENTAGE ERROR
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


    //------------------------------

    TGraphErrors *gV_GAIN_1  = new TGraphErrors(n_GAIN_1, HV_1, GAIN_1, errHV_1, errGAIN_1);
    TGraphErrors *gV_GAIN_2  = new TGraphErrors(n_GAIN_2, HV_2, GAIN_2, errHV_2, errGAIN_2);
    TGraphErrors *gV_GAIN_3  = new TGraphErrors(n_GAIN_3, HV_3, GAIN_3, errHV_3, errGAIN_3);


    //------------------------------

    gV_GAIN_1->SetMarkerStyle(20);
    gV_GAIN_1->SetMarkerColor(kOrange+2);
    gV_GAIN_1->SetTitle();
    gV_GAIN_1->GetXaxis()->SetTitle("Bias Voltage (V)");
    gV_GAIN_1->GetYaxis()->SetTitle("GAIN (mV)");

    gV_GAIN_2->SetMarkerStyle(20);
    gV_GAIN_2->SetMarkerColor(kOrange+2);
    gV_GAIN_2->SetTitle();
    gV_GAIN_2->GetXaxis()->SetTitle("Bias Voltage (V)");
    gV_GAIN_2->GetYaxis()->SetTitle("GAIN (mV)");

    gV_GAIN_3->SetMarkerStyle(20);
    gV_GAIN_3->SetMarkerColor(kOrange+2);
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
