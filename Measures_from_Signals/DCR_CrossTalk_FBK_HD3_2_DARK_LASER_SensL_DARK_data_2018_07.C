/******************************************************************************\
 *  DCR_CrossTalk_FBK_HD3_2_from_cnt_data_2018_07.C
 *             +
 *  DCR_CrossTalk_FBK_HD3_2_LASER_data_2018_07.C
 *             +
 *  DCR_CrossTalk_SensL_MicroFJ_SMTPA_60035_data_2018_07.C
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

#define w 1000
#define h 600


void DCR_CrossTalk_FBK_HD3_2_DARK_LASER_SensL_DARK_data_2018_07(){

    double HV_FBK_DARK[] = {31, 32, 33, 34, 35, 36, 37};
    double errHV_FBK_DARK[] = {0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01};
    double DCR_1_FBK_DARK[] = {301.3, 410.681, 493.403, 563.636, 630.558, 692.942, 756.225};
    double errDCR_1_FBK_DARK[] = {6.026, 8.21361, 9.86806, 11.2727, 12.6112, 13.8588, 26.2333};
    double DCR_Del_1_FBK_DARK[] = {296.506, 421.242, 512.722, 587.747, 658.806, 722.903, 789.519};
    double errDCR_Del_1_FBK_DARK[] = {5.93011, 8.42483, 10.2544, 11.7549, 13.1761, 15.3333, 28.925};
    double DCR_2_FBK_DARK[] = {303.656, 420.208, 524.233, 602.572, 675.086, 743.881, 807.983};
    double errDCR_2_FBK_DARK[] = {6.725, 15.8528, 22.8611, 18.8222, 18.6222, 26.7056, 41.1667};
    double DCR_Del_2_FBK_DARK[] = {297.931, 427.972, 541.625, 627.806, 704.772, 777.228, 844.7};
    double errDCR_Del_2_FBK_DARK[] = {6.72778, 16.6806, 23.4778, 19.2167, 19.025, 27.8667, 44.8611};
    double DCR_3_FBK_DARK[] = {367.839, 490.181, 589.292, 674.222, 754.314, 829.019, 897.669};
    double errDCR_3_FBK_DARK[] = {7.35678, 9.80361, 11.7858, 13.4844, 15.0863, 16.5804, 34.9056};
    double DCR_Del_3_FBK_DARK[] = {364.106, 502.686, 611.697, 703.564, 787.081, 866.672, 940.539};
    double errDCR_Del_3_FBK_DARK[] = {7.28211, 10.0537, 12.2339, 14.0713, 15.7416, 17.3334, 37.8861};

    double CT_1_FBK_LASER[] = {0.182358, 0.197717, 0.212695, 0.226478, 0.246893, 0.291059};
    double errCT_1_FBK_LASER[] = {0.02, 0.02, 0.02, 0.02, 0.02, 0.02};

    double HV_SensL[] = {28, 29, 30, 31};
    double DCR_1_SensL[] = {103.836, 164.281, 219.954, 280.063};
    double errDCR_1_SensL[] = {7.70229, 6.33222, 12.2104, 20.2025};
    double DCR_Del_1_SensL[] = {92.4948, 175.263, 246.02, 327.592};
    double errDCR_Del_1_SensL[] = {8.43455, 11.5954, 24.602, 28.2834};
    double DCR_2_SensL[] = {154.282, 202.5, 246.82, 315.683};
    double errDCR_2_SensL[] = {3.08564, 4.05, 11.9116, 14.3331};
    double DCR_Del_2_SensL[] = {137.997, 214.255, 280.717, 365.133};
    double errDCR_Del_2_SensL[] = {2.75993, 4.28511, 28.0717, 34.4145};

    int n_DCR_DARK = 7;

    TGraphErrors *gDCR_1_FBK_DARK  = new TGraphErrors(n_DCR_DARK, HV_FBK_DARK, DCR_1_FBK_DARK, errHV_FBK_DARK, errDCR_1_FBK_DARK);
    TGraphErrors *gDCR_2_FBK_DARK  = new TGraphErrors(n_DCR_DARK, HV_FBK_DARK, DCR_2_FBK_DARK, errHV_FBK_DARK, errDCR_2_FBK_DARK);
    TGraphErrors *gDCR_3_FBK_DARK  = new TGraphErrors(n_DCR_DARK, HV_FBK_DARK, DCR_3_FBK_DARK, errHV_FBK_DARK, errDCR_3_FBK_DARK);
    TGraphErrors *gDCR_Del_1_FBK_DARK  = new TGraphErrors(n_DCR_DARK, HV_FBK_DARK, DCR_Del_1_FBK_DARK, errHV_FBK_DARK, errDCR_Del_1_FBK_DARK);
    TGraphErrors *gDCR_Del_2_FBK_DARK  = new TGraphErrors(n_DCR_DARK, HV_FBK_DARK, DCR_Del_2_FBK_DARK, errHV_FBK_DARK, errDCR_Del_2_FBK_DARK);
    TGraphErrors *gDCR_Del_3_FBK_DARK  = new TGraphErrors(n_DCR_DARK, HV_FBK_DARK, DCR_Del_3_FBK_DARK, errHV_FBK_DARK, errDCR_Del_3_FBK_DARK);
    // TGraphErrors *gCT_1_FBK_DARK  = new TGraphErrors(n_DCR_DARK, HV_FBK_DARK, CT_1_FBK_DARK, errHV_FBK_DARK, errCT_1_FBK_DARK);
    // TGraphErrors *gCT_2_FBK_DARK  = new TGraphErrors(n_DCR_DARK, HV_FBK_DARK, CT_2_FBK_DARK, errHV_FBK_DARK, errCT_2_FBK_DARK);
    // TGraphErrors *gCT_3_FBK_DARK  = new TGraphErrors(n_DCR_DARK, HV_FBK_DARK, CT_3_FBK_DARK, errHV_FBK_DARK, errCT_3_FBK_DARK);
    // TGraphErrors *gCT_Del_1_FBK_DARK  = new TGraphErrors(n_DCR_DARK, HV_FBK_DARK, CT_Del_1_FBK_DARK, errHV_FBK_DARK, errCT_Del_1_FBK_DARK);
    // TGraphErrors *gCT_Del_2_FBK_DARK  = new TGraphErrors(n_DCR_DARK, HV_FBK_DARK, CT_Del_2_FBK_DARK, errHV_FBK_DARK, errCT_Del_2_FBK_DARK);
    // TGraphErrors *gCT_Del_3_FBK_DARK  = new TGraphErrors(n_DCR_DARK, HV_FBK_DARK, CT_Del_3_FBK_DARK, errHV_FBK_DARK, errCT_Del_3_FBK_DARK);


    gDCR_1_FBK_DARK->SetMarkerStyle(20);
    gDCR_1_FBK_DARK->SetMarkerColor(kOrange+1);
    gDCR_1_FBK_DARK->SetTitle();
    gDCR_1_FBK_DARK->GetXaxis()->SetTitle("Bias Voltage (V)");
    gDCR_1_FBK_DARK->GetYaxis()->SetTitle("title_DCR");

    gDCR_2_FBK_DARK->SetMarkerStyle(20);
    gDCR_2_FBK_DARK->SetMarkerColor(kRed);
    gDCR_2_FBK_DARK->SetTitle();
    gDCR_2_FBK_DARK->GetXaxis()->SetTitle("Bias Voltage (V)");
    gDCR_2_FBK_DARK->GetYaxis()->SetTitle("title_DCR");

    gDCR_3_FBK_DARK->SetMarkerStyle(20);
    gDCR_3_FBK_DARK->SetMarkerColor(kMagenta);
    gDCR_3_FBK_DARK->SetTitle();
    gDCR_3_FBK_DARK->GetXaxis()->SetTitle("Bias Voltage (V)");
    gDCR_3_FBK_DARK->GetYaxis()->SetTitle("title_DCR");

    //------------------------------

    // gCT_1_FBK_DARK->SetMarkerStyle(20);
    // gCT_1_FBK_DARK->SetMarkerColor(kOrange+1);
    // gCT_1_FBK_DARK->SetTitle();
    // gCT_1_FBK_DARK->GetXaxis()->SetTitle("Bias Voltage (V)");
    // gCT_1_FBK_DARK->GetYaxis()->SetTitle("P_{CT}");
    //
    // gCT_2_FBK_DARK->SetMarkerStyle(20);
    // gCT_2_FBK_DARK->SetMarkerColor(kRed);
    // gCT_2_FBK_DARK->SetTitle();
    // gCT_2_FBK_DARK->GetXaxis()->SetTitle("Bias Voltage (V)");
    // gCT_2_FBK_DARK->GetYaxis()->SetTitle("P_{CT}");
    //
    // gCT_3_FBK_DARK->SetMarkerStyle(20);
    // gCT_3_FBK_DARK->SetMarkerColor(kMagenta);
    // gCT_3_FBK_DARK->SetTitle();
    // gCT_3_FBK_DARK->GetXaxis()->SetTitle("Bias Voltage (V)");
    // gCT_3_FBK_DARK->GetYaxis()->SetTitle("P_{CT}");





    TCanvas *cV_CT_FBK = new TCanvas("cV_CT_FBK", "cV_CT_FBK",w,h);
    cV_CT_FBK->SetGrid();
    TMultiGraph *mgV_CT_Del_1_E = new TMultiGraph("mgV_CT_Del_1_E", ";Bias Voltage (V); CT_Del");
    mgV_CT_Del_1_E->Add(gDCR_1_FBK_DARK);
    mgV_CT_Del_1_E->Draw("AP");

}
