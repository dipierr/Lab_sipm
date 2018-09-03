/******************************************************************************\
 * GAIN_V_SiPM_SensL_MicroFJ_SMTPA_GausSum_NoGap.C
 *
 * GAIN values obtained by Ana_Traces_SiPM.cxx (version of 02/09/2018)
 *
 * KEY POINTS:
 *  > Ana1(...)
 *  > dleddt = 8
 *  > NO trace smoothing
 *  > thr (parameter of Ana1(...)) set depending on the situation
 *  > fitted using gaus_sum_12, i.e.:
 *      "[0]*TMath::Exp( - (x-[2]-[3])*(x-[2]-[3])/( 2*([4]*[4] + [5]*[5] )) )
 *       + [1]*TMath::Exp( - (x-[2]-2*[3])*(x-[2]-2*[3])/( 2*([4]*[4] +
 *      4*[5]*[5] ) ) ) "
 *  > fit range set manually
 *
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

#define n_GAIN_1 4
#define n_GAIN_2 4
#define n_GAIN_3 1 // not used


#define h 600
#define w 1000

void GAIN_V_SiPM_SensL_MicroFJ_SMTPA_GausSum_NoGap();
int find_index(double vect[],int N, double value);

void GAIN_V_SiPM_SensL_MicroFJ_SMTPA_GausSum_NoGap(){

    // HV:
    // SiPM1:
    double HV_1[n_GAIN_1], errHV_1[n_GAIN_1];
    double GAIN_1[n_GAIN_1], errGAIN_1[n_GAIN_1];
    // SiPM2:
    double HV_2[n_GAIN_2], errHV_2[n_GAIN_2];
    double GAIN_2[n_GAIN_2], errGAIN_2[n_GAIN_2];

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


    // HV
    HV_1[0]    = 28.00;
    errHV_1[0] =  0.01;
    for(int i=1; i<n_GAIN_1; i++){
        HV_1[i]    = HV_1[i-1]+1.;
        errHV_1[i] = errHV_1[0];
    }

    HV_2[0]    = 28.00;
    errHV_2[0] =  0.01;
    for(int i=1; i<n_GAIN_2; i++){
        HV_2[i]    = HV_2[i-1]+1.;
        errHV_2[i] = errHV_2[0];
    }

    ///////////////////////////////////////////////////////////////////////////
    //      SiPM1
    ///////////////////////////////////////////////////////////////////////////
    HV = 28.;
    GAIN_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 7.04621e+00;//7.05501e+00;
    errGAIN_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 2.88865e-02;//2.78957e-02;

    HV = 29;
    GAIN_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 1.05541e+01;//1.05494e+01;
    errGAIN_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 1.13549e-02;//1.17454e-02;

    HV = 30;
    GAIN_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 1.35555e+01;//1.31027e+01;
    errGAIN_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 1.02567e-02;//3.20972e-02;

    HV = 31;
    GAIN_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 1.64390e+01;//1.60932e+01;
    errGAIN_1[find_index(HV_1,  sizeof(HV_1)/sizeof(HV_1[0]), HV)] = 1.15212e-02;//1.19408e-02;




    ///////////////////////////////////////////////////////////////////////////
    //      SiPM2
    ///////////////////////////////////////////////////////////////////////////
    HV = 28.;
    index = find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV);
    // cout<< index <<endl;
    GAIN_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 9.97040e+00;//1.05731e+01;
    errGAIN_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 1.85964e-02;//2.29943e-02;

    HV = 29;
    GAIN_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 1.50591e+01;//1.50629e+01;
    errGAIN_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 1.38029e-02;//1.38306e-02 ;

    HV = 30;
    GAIN_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 1.92912e+01;//1.92838e+01;
    errGAIN_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 1.26277e-02;//1.28872e-02;

    HV = 31;
    GAIN_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 2.35660e+01;//2.32833e+01;
    errGAIN_2[find_index(HV_2,  sizeof(HV_2)/sizeof(HV_2[0]), HV)] = 1.46330e-02;//1.52857e-02;




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

    }
    if(fix_error_bool){
        double err_fix = 1;
        for(int i=0; i<n_GAIN_1; i++){
            errGAIN_1[i] = err_fix;
        }
        for(int i=0; i<n_GAIN_2; i++){
            errGAIN_2[i] = err_fix;
        }

    }



    //------------------------------


    TGraphErrors *gV_GAIN_1  = new TGraphErrors(n_GAIN_1, HV_1, GAIN_1, errHV_1, errGAIN_1);
    TGraphErrors *gV_GAIN_2  = new TGraphErrors(n_GAIN_2, HV_2, GAIN_2, errHV_2, errGAIN_2);


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



    //------------------------------

    TCanvas *cV_GAIN_1 = new TCanvas("cV_GAIN_1", "cV_GAIN_1",w,h);
    cV_GAIN_1->SetGrid();
    gV_GAIN_1->Draw("AP");

    TCanvas *cV_GAIN_2 = new TCanvas("cV_GAIN_2", "cV_GAIN_2",w,h);
    cV_GAIN_2->SetGrid();
    gV_GAIN_2->Draw("AP");



    //------------------------------


    double V_bd_1,V_bd_2,V_bd_3;

    TF1 *linearFit1 = new TF1("linearFit1","pol1",-100,100);
    TF1 *linearFit2 = new TF1("linearFit2","pol1",-100,100);


    gV_GAIN_1->Fit("linearFit1");
    V_bd_1  = - linearFit1->GetParameter(0)/linearFit1->GetParameter(1);

    gV_GAIN_2->Fit("linearFit2");
    V_bd_2  = - linearFit2->GetParameter(0)/linearFit2->GetParameter(1);





    cout<<"V_bd_1 = "<<V_bd_1<<endl;
    cout<<"V_bd_2 = "<<V_bd_2<<endl;

    auto legendGAIN = new TLegend(0.15,0.70,0.35,0.85);
    legendGAIN->AddEntry(gV_GAIN_1,"SensL MicroFJ (1)","p");
    legendGAIN->AddEntry(gV_GAIN_2,"SensL MicroFJ (2)","p");


    TCanvas *cGAIN = new TCanvas("cGAIN", "cGAIN",w,h);
    cGAIN->SetGrid();
    TMultiGraph *mgGAIN = new TMultiGraph("mgGAIN", ";Bias Voltage (V);GAIN (mV)");
    mgGAIN->Add(gV_GAIN_1);
    mgGAIN->Add(gV_GAIN_2);
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


    for(int i=0; i<n_GAIN_1; i++){
        HV_Vbd_1[i+1]      = HV_1[i];
        errHV_Vbd_1[i+1]   = errHV_1[i];
        GAIN_Vbd_1[i+1]    = GAIN_1[i];
        errGAIN_Vbd_1[i+1] = errGAIN_1[i];

        HV_Vbd_2[i+1]      = HV_2[i];
        errHV_Vbd_2[i+1]   = errHV_2[i];
        GAIN_Vbd_2[i+1]    = GAIN_2[i];
        errGAIN_Vbd_2[i+1] = errGAIN_2[i];


    }


    TGraphErrors *gV_GAIN_Vbd_1  = new TGraphErrors(n_GAIN_1+1, HV_Vbd_1, GAIN_Vbd_1, errHV_Vbd_1, errGAIN_Vbd_1);
    TGraphErrors *gV_GAIN_Vbd_2  = new TGraphErrors(n_GAIN_2+1, HV_Vbd_2, GAIN_Vbd_2, errHV_Vbd_2, errGAIN_Vbd_2);


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


    //------------------------------

    TCanvas *cV_GAIN_Vbd_1 = new TCanvas("cV_GAIN_Vbd_1", "cV_GAIN_Vbd_1",w,h);
    cV_GAIN_Vbd_1->SetGrid();
    gV_GAIN_Vbd_1->Draw("AP");

    TCanvas *cV_GAIN_Vbd_2 = new TCanvas("cV_GAIN_Vbd_2", "cV_GAIN_Vbd_2",w,h);
    cV_GAIN_Vbd_2->SetGrid();
    gV_GAIN_Vbd_2->Draw("AP");



    //------------------------------



    TCanvas *cGAIN_Vbd = new TCanvas("cGAIN_Vbd", "cGAIN_Vbd",w,h);
    cGAIN_Vbd->SetGrid();
    TMultiGraph *mgGAIN_Vdb = new TMultiGraph("mgGAIN_Vdb", ";Bias Voltage (V);GAIN (mV)");
    mgGAIN_Vdb->Add(gV_GAIN_Vbd_1);
    mgGAIN_Vdb->Add(gV_GAIN_Vbd_2);
    mgGAIN_Vdb->Draw("AP");


    //------------------------------
    // For LaTeX
    //------------------------------
    cout<<endl<<endl;
    cout<<"For LaTeX"<<endl<<endl;
    cout<<"\% From file GAIN_V_SiPM_SensL_MicroFJ_SMTPA_GausSum.C"<<endl;
    cout<<"\% HV & GAIN_1 (mV) & GAIN_2 (mV) \\\\ "<<endl;
    for(int i=0; i<n_GAIN_1; i++){
        printf("$ %.2f \\pm %.2f $ & $ %.2f \\pm %.2f $ & $ %.2f \\pm %.2f $ \\\\ \n", HV_1[i], errHV_1[i], GAIN_1[i], errGAIN_1[i],GAIN_2[i], errGAIN_2[i]);
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
