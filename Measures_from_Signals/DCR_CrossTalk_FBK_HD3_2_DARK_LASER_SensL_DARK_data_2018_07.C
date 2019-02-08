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

    //---------------
    // FBK HD3-2 DARK
    //---------------
    double HV_FBK_DARK[] = {31, 32, 33, 34, 35, 36, 37};
    double errHV_FBK_DARK[] = {0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01};

    // DCR
    double DCR_1_FBK_DARK[] = {301.3, 410.681, 493.403, 563.636, 630.558, 692.942, 756.225};
    double errDCR_1_FBK_DARK[] = {3.50833, 2.15556, 1.21389, 1.85278, 1.73611, 1.27778, 1.49444};
    double DCR_Del_1_FBK_DARK[] = {296.506, 421.242, 512.722, 587.747, 658.806, 722.903, 789.519};
    double errDCR_Del_1_FBK_DARK[] = {3.65, 1.86111, 0.891667, 1.41667, 1.38611, 1.28889, 1.71667};
    double DCR_2_FBK_DARK[] = {303.656, 416.183, 520.406, 596.392, 670.1, 739.65, 805.164};
    double errDCR_2_FBK_DARK[] = {6.725, 11.8278, 9.06944, 6.18056, 4.98611, 4.23056, 3.08056};
    double DCR_Del_2_FBK_DARK[] = {297.928, 423.736, 537.936, 621.978, 699.983, 773.25, 841.906};
    double errDCR_Del_2_FBK_DARK[] = {6.725, 12.4444, 9.06944, 5.86389, 4.78889, 3.97778, 3.31389};
    double DCR_3_FBK_DARK[] = {367.839, 490.181, 589.292, 674.222, 754.314, 829.019, 897.669};
    double errDCR_3_FBK_DARK[] = {3.89167, 1.88889, 2.89167, 2.37222, 3.18611, 3.66667, 6.18889};
    double DCR_Del_3_FBK_DARK[] = {364.106, 502.686, 611.697, 703.564, 787.081, 866.672, 940.539};
    double errDCR_Del_3_FBK_DARK[] = {4.41389, 1.74722, 2.84167, 2.61667, 3.625, 3.99444, 6.625};

    // CT
    double CT_1_FBK_DARK[] = {0.194892, 0.20516, 0.240606, 0.269749, 0.29342, 0.328475, 0.358377};
    double errCT_1_FBK_DARK[] = {0.016018, 0.013109, 0.011335, 0.01174, 0.014037, 0.012211, 0.012549};
    double CT_Del_1_FBK_DARK[] = {0.197997, 0.22506, 0.263293, 0.28729, 0.306855, 0.34085, 0.366474};
    double errCT_Del_1_FBK_DARK[] = {0.017453, 0.023314, 0.015238, 0.018588, 0.026278, 0.018788, 0.018547};
    double CT_2_FBK_DARK[] = {0.176788, 0.206578, 0.236805, 0.263922, 0.289775, 0.328639, 0.354963};
    double errCT_2_FBK_DARK[] = {0.018846, 0.021822, 0.022803, 0.016043, 0.018487, 0.019276, 0.015744};
    double CT_Del_2_FBK_DARK[] = {0.185645, 0.226305, 0.25734, 0.277028, 0.299984, 0.337352, 0.359528};
    double errCT_Del_2_FBK_DARK[] = {0.024479, 0.024405, 0.023906, 0.026155, 0.03034, 0.01872, 0.021716};
    double CT_3_FBK_DARK[] = {0.192589, 0.208325, 0.24549, 0.276746, 0.306715, 0.339257, 0.369078};
    double errCT_3_FBK_DARK[] = {0.016406, 0.013547, 0.014306, 0.011511, 0.012366, 0.011901, 0.013996};
    double CT_Del_3_FBK_DARK[] = {0.182632, 0.223419, 0.260092, 0.286346, 0.31063, 0.341269, 0.366328};
    double errCT_Del_3_FBK_DARK[] = {0.021971, 0.028105, 0.020238, 0.018358, 0.019569, 0.017404, 0.021229};

    //----------------
    // FBK HD3-2 LASER
    //----------------
    double HV_FBK_LASER[] = {31, 32, 33, 34, 35, 36};
    double errHV_FBK_LASER[] = {0.01, 0.01, 0.01, 0.01, 0.01, 0.01};
    double CT_1_FBK_LASER[] = {0.182358, 0.197717, 0.212695, 0.226478, 0.246893, 0.291059};
    double errCT_1_FBK_LASER[] = {0.02, 0.02, 0.02, 0.02, 0.02, 0.02};



    //-------------
    // SensL - DARK
    //-------------
    double HV_SensL[] = {28, 29, 30, 31};
    double errHV_SensL[] = {0.01, 0.01, 0.01, 0.01};
    double DCR_1_SensL[] = {103.836, 164.281, 219.954, 280.063};
    double errDCR_1_SensL[] = {7.70229, 6.33168, 3.87109, 7.75413};
    double DCR_Del_1_SensL[] = {91.5668, 175.066, 246.022, 327.595};
    double errDCR_Del_1_SensL[] = {8.83107, 11.5544, 9.78263, 11.4561};
    double DCR_2_SensL[] = {154.282, 202.5, 246.82, 315.683};
    double errDCR_2_SensL[] = {7.74381, 6.15852, 3.80948, 3.16462};
    double DCR_Del_2_SensL[] = {137.687, 214.212, 280.723, 365.136};
    double errDCR_Del_2_SensL[] = {3.24902, 1.2824, 3.60158, 5.48787};
    double CT_1_SensL[] = {0.183342, 0.217608, 0.285169, 0.370261};
    double errCT_1_SensL[] = {0.011747, 0.01026, 0.008513, 0.007842};
    double CT_Del_1_SensL[] = {0.502512, 0.464136, 0.523566, 0.681726};
    double errCT_Del_1_SensL[] = {0.017802, 0.041348, 0.038044, 0.01485};
    double CT_2_SensL[] = {0.122127, 0.188606, 0.263666, 0.332936};
    double errCT_2_SensL[] = {0.004095, 0.001001, 0.003332, 0.00565};
    double CT_Del_2_SensL[] = {0.0635249, 0.304844, 0.422022, 0.515758};
    double errCT_Del_2_SensL[] = {0.0297699, 0.039856, 0.026187, 0.024853};



    int n_DCR_DARK = 7;
    int n_DCR_LASER = 6;
    int n_DCR_SensL = 4;

    double V_bd_FBK = 27;
    double V_bd_SensL = 25.5;

    double HV_FBK_DARK_OV[7];
    double HV_FBK_LASER_OV[6];
    double HV_SensL_OV[4];

    for(int i=0; i<n_DCR_DARK; i++){
        HV_FBK_DARK_OV[i]=HV_FBK_DARK[i]-V_bd_FBK;
    }

    for(int i=0; i<n_DCR_LASER; i++){
        HV_FBK_LASER_OV[i]=HV_FBK_LASER[i]-V_bd_FBK;
    }

    for(int i=0; i<n_DCR_SensL; i++){
        HV_SensL_OV[i]=HV_SensL[i]-V_bd_SensL;
    }


    TGraphErrors *gDCR_1_FBK_DARK  = new TGraphErrors(n_DCR_DARK, HV_FBK_DARK, DCR_1_FBK_DARK, errHV_FBK_DARK, errDCR_1_FBK_DARK);
    TGraphErrors *gDCR_2_FBK_DARK  = new TGraphErrors(n_DCR_DARK, HV_FBK_DARK, DCR_2_FBK_DARK, errHV_FBK_DARK, errDCR_2_FBK_DARK);
    TGraphErrors *gDCR_3_FBK_DARK  = new TGraphErrors(n_DCR_DARK, HV_FBK_DARK, DCR_3_FBK_DARK, errHV_FBK_DARK, errDCR_3_FBK_DARK);
    TGraphErrors *gDCR_Del_1_FBK_DARK  = new TGraphErrors(n_DCR_DARK, HV_FBK_DARK, DCR_Del_1_FBK_DARK, errHV_FBK_DARK, errDCR_Del_1_FBK_DARK);
    TGraphErrors *gDCR_Del_2_FBK_DARK  = new TGraphErrors(n_DCR_DARK, HV_FBK_DARK, DCR_Del_2_FBK_DARK, errHV_FBK_DARK, errDCR_Del_2_FBK_DARK);
    TGraphErrors *gDCR_Del_3_FBK_DARK  = new TGraphErrors(n_DCR_DARK, HV_FBK_DARK, DCR_Del_3_FBK_DARK, errHV_FBK_DARK, errDCR_Del_3_FBK_DARK);
    TGraphErrors *gCT_1_FBK_DARK  = new TGraphErrors(n_DCR_DARK, HV_FBK_DARK, CT_1_FBK_DARK, errHV_FBK_DARK, errCT_1_FBK_DARK);
    TGraphErrors *gCT_2_FBK_DARK  = new TGraphErrors(n_DCR_DARK, HV_FBK_DARK, CT_2_FBK_DARK, errHV_FBK_DARK, errCT_2_FBK_DARK);
    TGraphErrors *gCT_3_FBK_DARK  = new TGraphErrors(n_DCR_DARK, HV_FBK_DARK, CT_3_FBK_DARK, errHV_FBK_DARK, errCT_3_FBK_DARK);
    TGraphErrors *gCT_Del_1_FBK_DARK  = new TGraphErrors(n_DCR_DARK, HV_FBK_DARK, CT_Del_1_FBK_DARK, errHV_FBK_DARK, errCT_Del_1_FBK_DARK);
    TGraphErrors *gCT_Del_2_FBK_DARK  = new TGraphErrors(n_DCR_DARK, HV_FBK_DARK, CT_Del_2_FBK_DARK, errHV_FBK_DARK, errCT_Del_2_FBK_DARK);
    TGraphErrors *gCT_Del_3_FBK_DARK  = new TGraphErrors(n_DCR_DARK, HV_FBK_DARK, CT_Del_3_FBK_DARK, errHV_FBK_DARK, errCT_Del_3_FBK_DARK);

    TGraphErrors *gCT_1_FBK_LASER  = new TGraphErrors(n_DCR_LASER, HV_FBK_LASER, CT_1_FBK_LASER, errHV_FBK_LASER, errCT_1_FBK_LASER);

    TGraphErrors *gDCR_1_SensL  = new TGraphErrors(n_DCR_SensL, HV_SensL, DCR_1_SensL, errHV_SensL, errDCR_1_SensL);
    TGraphErrors *gDCR_2_SensL  = new TGraphErrors(n_DCR_SensL, HV_SensL, DCR_2_SensL, errHV_SensL, errDCR_2_SensL);
    TGraphErrors *gDCR_Del_1_SensL  = new TGraphErrors(n_DCR_SensL, HV_SensL, DCR_Del_1_SensL, errHV_SensL, errDCR_Del_1_SensL);
    TGraphErrors *gDCR_Del_2_SensL  = new TGraphErrors(n_DCR_SensL, HV_SensL, DCR_Del_2_SensL, errHV_SensL, errDCR_Del_2_SensL);
    TGraphErrors *gCT_1_SensL  = new TGraphErrors(n_DCR_SensL, HV_SensL, CT_1_SensL, errHV_SensL, errCT_1_SensL);
    TGraphErrors *gCT_2_SensL  = new TGraphErrors(n_DCR_SensL, HV_SensL, CT_2_SensL, errHV_SensL, errCT_2_SensL);
    TGraphErrors *gCT_Del_1_SensL  = new TGraphErrors(n_DCR_SensL, HV_SensL, CT_Del_1_SensL, errHV_SensL, errCT_Del_1_SensL);
    TGraphErrors *gCT_Del_2_SensL  = new TGraphErrors(n_DCR_SensL, HV_SensL, CT_Del_2_SensL, errHV_SensL, errCT_Del_2_SensL);

    TGraphErrors *gDCR_1_FBK_DARK_OV  = new TGraphErrors(n_DCR_DARK, HV_FBK_DARK_OV, DCR_1_FBK_DARK, errHV_FBK_DARK, errDCR_1_FBK_DARK);
    TGraphErrors *gDCR_2_FBK_DARK_OV  = new TGraphErrors(n_DCR_DARK, HV_FBK_DARK_OV, DCR_2_FBK_DARK, errHV_FBK_DARK, errDCR_2_FBK_DARK);
    TGraphErrors *gDCR_3_FBK_DARK_OV  = new TGraphErrors(n_DCR_DARK, HV_FBK_DARK_OV, DCR_3_FBK_DARK, errHV_FBK_DARK, errDCR_3_FBK_DARK);
    TGraphErrors *gDCR_Del_1_FBK_DARK_OV  = new TGraphErrors(n_DCR_DARK, HV_FBK_DARK_OV, DCR_Del_1_FBK_DARK, errHV_FBK_DARK, errDCR_Del_1_FBK_DARK);
    TGraphErrors *gDCR_Del_2_FBK_DARK_OV  = new TGraphErrors(n_DCR_DARK, HV_FBK_DARK_OV, DCR_Del_2_FBK_DARK, errHV_FBK_DARK, errDCR_Del_2_FBK_DARK);
    TGraphErrors *gDCR_Del_3_FBK_DARK_OV  = new TGraphErrors(n_DCR_DARK, HV_FBK_DARK_OV, DCR_Del_3_FBK_DARK, errHV_FBK_DARK, errDCR_Del_3_FBK_DARK);
    TGraphErrors *gCT_1_FBK_DARK_OV  = new TGraphErrors(n_DCR_DARK, HV_FBK_DARK_OV, CT_1_FBK_DARK, errHV_FBK_DARK, errCT_1_FBK_DARK);
    TGraphErrors *gCT_2_FBK_DARK_OV  = new TGraphErrors(n_DCR_DARK, HV_FBK_DARK_OV, CT_2_FBK_DARK, errHV_FBK_DARK, errCT_2_FBK_DARK);
    TGraphErrors *gCT_3_FBK_DARK_OV  = new TGraphErrors(n_DCR_DARK, HV_FBK_DARK_OV, CT_3_FBK_DARK, errHV_FBK_DARK, errCT_3_FBK_DARK);
    TGraphErrors *gCT_Del_1_FBK_DARK_OV  = new TGraphErrors(n_DCR_DARK, HV_FBK_DARK_OV, CT_Del_1_FBK_DARK, errHV_FBK_DARK, errCT_Del_1_FBK_DARK);
    TGraphErrors *gCT_Del_2_FBK_DARK_OV  = new TGraphErrors(n_DCR_DARK, HV_FBK_DARK_OV, CT_Del_2_FBK_DARK, errHV_FBK_DARK, errCT_Del_2_FBK_DARK);
    TGraphErrors *gCT_Del_3_FBK_DARK_OV  = new TGraphErrors(n_DCR_DARK, HV_FBK_DARK_OV, CT_Del_3_FBK_DARK, errHV_FBK_DARK, errCT_Del_3_FBK_DARK);

    TGraphErrors *gCT_1_FBK_LASER_OV  = new TGraphErrors(n_DCR_LASER, HV_FBK_LASER_OV, CT_1_FBK_LASER, errHV_FBK_LASER, errCT_1_FBK_LASER);

    TGraphErrors *gDCR_1_SensL_OV  = new TGraphErrors(n_DCR_SensL, HV_SensL_OV, DCR_1_SensL, errHV_SensL, errDCR_1_SensL);
    TGraphErrors *gDCR_2_SensL_OV  = new TGraphErrors(n_DCR_SensL, HV_SensL_OV, DCR_2_SensL, errHV_SensL, errDCR_2_SensL);
    TGraphErrors *gDCR_Del_1_SensL_OV  = new TGraphErrors(n_DCR_SensL, HV_SensL_OV, DCR_Del_1_SensL, errHV_SensL, errDCR_Del_1_SensL);
    TGraphErrors *gDCR_Del_2_SensL_OV  = new TGraphErrors(n_DCR_SensL, HV_SensL_OV, DCR_Del_2_SensL, errHV_SensL, errDCR_Del_2_SensL);
    TGraphErrors *gCT_1_SensL_OV  = new TGraphErrors(n_DCR_SensL, HV_SensL_OV, CT_1_SensL, errHV_SensL, errCT_1_SensL);
    TGraphErrors *gCT_2_SensL_OV  = new TGraphErrors(n_DCR_SensL, HV_SensL_OV, CT_2_SensL, errHV_SensL, errCT_2_SensL);
    TGraphErrors *gCT_Del_1_SensL_OV  = new TGraphErrors(n_DCR_SensL, HV_SensL_OV, CT_Del_1_SensL, errHV_SensL, errCT_Del_1_SensL);
    TGraphErrors *gCT_Del_2_SensL_OV  = new TGraphErrors(n_DCR_SensL, HV_SensL_OV, CT_Del_2_SensL, errHV_SensL, errCT_Del_2_SensL);


    //------------------------------

    char title_P_CT_mg[] = ";Bias Voltage (V); P_{CT}";
    char title_P_CT[] = "P_{CT}";
    char title_DCR_Area[] = "\\frac{DCR}{Area} \\left(\\frac{kHz}{mm^2}\\right)";
    char title_DCR_Area_mg[] = ";Bias Voltage (V); \\frac{DCR}{Area} \\left(\\frac{kHz}{mm^2}\\right)";

    //------------------------------

    gDCR_1_FBK_DARK->SetMarkerStyle(20);
    gDCR_1_FBK_DARK->SetMarkerColor(kOrange+1);
    gDCR_1_FBK_DARK->SetTitle();
    gDCR_1_FBK_DARK->GetXaxis()->SetTitle("Bias Voltage (V)");
    gDCR_1_FBK_DARK->GetYaxis()->SetTitle(title_DCR_Area);

    gDCR_2_FBK_DARK->SetMarkerStyle(20);
    gDCR_2_FBK_DARK->SetMarkerColor(kRed);
    gDCR_2_FBK_DARK->SetTitle();
    gDCR_2_FBK_DARK->GetXaxis()->SetTitle("Bias Voltage (V)");
    gDCR_2_FBK_DARK->GetYaxis()->SetTitle(title_DCR_Area);

    gDCR_3_FBK_DARK->SetMarkerStyle(20);
    gDCR_3_FBK_DARK->SetMarkerColor(kMagenta);
    gDCR_3_FBK_DARK->SetTitle();
    gDCR_3_FBK_DARK->GetXaxis()->SetTitle("Bias Voltage (V)");
    gDCR_3_FBK_DARK->GetYaxis()->SetTitle(title_DCR_Area);

    //------------------------------

    gCT_1_FBK_DARK->SetMarkerStyle(20);
    gCT_1_FBK_DARK->SetMarkerColor(kOrange+1);
    gCT_1_FBK_DARK->SetTitle();
    gCT_1_FBK_DARK->GetXaxis()->SetTitle("Bias Voltage (V)");
    gCT_1_FBK_DARK->GetYaxis()->SetTitle(title_P_CT);

    gCT_2_FBK_DARK->SetMarkerStyle(20);
    gCT_2_FBK_DARK->SetMarkerColor(kRed);
    gCT_2_FBK_DARK->SetTitle();
    gCT_2_FBK_DARK->GetXaxis()->SetTitle("Bias Voltage (V)");
    gCT_2_FBK_DARK->GetYaxis()->SetTitle(title_P_CT);

    gCT_3_FBK_DARK->SetMarkerStyle(20);
    gCT_3_FBK_DARK->SetMarkerColor(kMagenta);
    gCT_3_FBK_DARK->SetTitle();
    gCT_3_FBK_DARK->GetXaxis()->SetTitle("Bias Voltage (V)");
    gCT_3_FBK_DARK->GetYaxis()->SetTitle(title_P_CT);

    //------------------------------

    gDCR_Del_1_FBK_DARK->SetMarkerStyle(22);
    gDCR_Del_1_FBK_DARK->SetMarkerColor(kOrange+1);
    gDCR_Del_1_FBK_DARK->SetTitle();
    gDCR_Del_1_FBK_DARK->GetXaxis()->SetTitle("Bias Voltage (V)");
    gDCR_Del_1_FBK_DARK->GetYaxis()->SetTitle(title_DCR_Area);

    gDCR_Del_2_FBK_DARK->SetMarkerStyle(22);
    gDCR_Del_2_FBK_DARK->SetMarkerColor(kRed);
    gDCR_Del_2_FBK_DARK->SetTitle();
    gDCR_Del_2_FBK_DARK->GetXaxis()->SetTitle("Bias Voltage (V)");
    gDCR_Del_2_FBK_DARK->GetYaxis()->SetTitle(title_DCR_Area);

    gDCR_Del_3_FBK_DARK->SetMarkerStyle(22);
    gDCR_Del_3_FBK_DARK->SetMarkerColor(kMagenta);
    gDCR_Del_3_FBK_DARK->SetTitle();
    gDCR_Del_3_FBK_DARK->GetXaxis()->SetTitle("Bias Voltage (V)");
    gDCR_Del_3_FBK_DARK->GetYaxis()->SetTitle(title_DCR_Area);

    //------------------------------

    gCT_Del_1_FBK_DARK->SetMarkerStyle(22);
    gCT_Del_1_FBK_DARK->SetMarkerColor(kOrange+1);
    gCT_Del_1_FBK_DARK->SetTitle();
    gCT_Del_1_FBK_DARK->GetXaxis()->SetTitle("Bias Voltage (V)");
    gCT_Del_1_FBK_DARK->GetYaxis()->SetTitle(title_P_CT);

    gCT_Del_2_FBK_DARK->SetMarkerStyle(22);
    gCT_Del_2_FBK_DARK->SetMarkerColor(kRed);
    gCT_Del_2_FBK_DARK->SetTitle();
    gCT_Del_2_FBK_DARK->GetXaxis()->SetTitle("Bias Voltage (V)");
    gCT_Del_2_FBK_DARK->GetYaxis()->SetTitle(title_P_CT);

    gCT_Del_3_FBK_DARK->SetMarkerStyle(22);
    gCT_Del_3_FBK_DARK->SetMarkerColor(kMagenta);
    gCT_Del_3_FBK_DARK->SetTitle();
    gCT_Del_3_FBK_DARK->GetXaxis()->SetTitle("Bias Voltage (V)");
    gCT_Del_3_FBK_DARK->GetYaxis()->SetTitle(title_P_CT);

    //------------------------------
    //------------------------------

    //------------------------------

    gDCR_1_FBK_DARK_OV->SetMarkerStyle(20);
    gDCR_1_FBK_DARK_OV->SetMarkerColor(kOrange+1);
    gDCR_1_FBK_DARK_OV->SetTitle();
    gDCR_1_FBK_DARK_OV->GetXaxis()->SetTitle("Bias Voltage (V)");
    gDCR_1_FBK_DARK_OV->GetYaxis()->SetTitle(title_DCR_Area);

    gDCR_2_FBK_DARK_OV->SetMarkerStyle(20);
    gDCR_2_FBK_DARK_OV->SetMarkerColor(kRed);
    gDCR_2_FBK_DARK_OV->SetTitle();
    gDCR_2_FBK_DARK_OV->GetXaxis()->SetTitle("Bias Voltage (V)");
    gDCR_2_FBK_DARK_OV->GetYaxis()->SetTitle(title_DCR_Area);

    gDCR_3_FBK_DARK_OV->SetMarkerStyle(20);
    gDCR_3_FBK_DARK_OV->SetMarkerColor(kMagenta);
    gDCR_3_FBK_DARK_OV->SetTitle();
    gDCR_3_FBK_DARK_OV->GetXaxis()->SetTitle("Bias Voltage (V)");
    gDCR_3_FBK_DARK_OV->GetYaxis()->SetTitle(title_DCR_Area);

    //------------------------------

    gCT_1_FBK_DARK_OV->SetMarkerStyle(20);
    gCT_1_FBK_DARK_OV->SetMarkerColor(kOrange+1);
    gCT_1_FBK_DARK_OV->SetTitle();
    gCT_1_FBK_DARK_OV->GetXaxis()->SetTitle("Bias Voltage (V)");
    gCT_1_FBK_DARK_OV->GetYaxis()->SetTitle(title_P_CT);

    gCT_2_FBK_DARK_OV->SetMarkerStyle(20);
    gCT_2_FBK_DARK_OV->SetMarkerColor(kRed);
    gCT_2_FBK_DARK_OV->SetTitle();
    gCT_2_FBK_DARK_OV->GetXaxis()->SetTitle("Bias Voltage (V)");
    gCT_2_FBK_DARK_OV->GetYaxis()->SetTitle(title_P_CT);

    gCT_3_FBK_DARK_OV->SetMarkerStyle(20);
    gCT_3_FBK_DARK_OV->SetMarkerColor(kMagenta);
    gCT_3_FBK_DARK_OV->SetTitle();
    gCT_3_FBK_DARK_OV->GetXaxis()->SetTitle("Bias Voltage (V)");
    gCT_3_FBK_DARK_OV->GetYaxis()->SetTitle(title_P_CT);

    //------------------------------

    gDCR_Del_1_FBK_DARK_OV->SetMarkerStyle(22);
    gDCR_Del_1_FBK_DARK_OV->SetMarkerColor(kOrange+1);
    gDCR_Del_1_FBK_DARK_OV->SetTitle();
    gDCR_Del_1_FBK_DARK_OV->GetXaxis()->SetTitle("Bias Voltage (V)");
    gDCR_Del_1_FBK_DARK_OV->GetYaxis()->SetTitle(title_DCR_Area);

    gDCR_Del_2_FBK_DARK_OV->SetMarkerStyle(22);
    gDCR_Del_2_FBK_DARK_OV->SetMarkerColor(kRed);
    gDCR_Del_2_FBK_DARK_OV->SetTitle();
    gDCR_Del_2_FBK_DARK_OV->GetXaxis()->SetTitle("Bias Voltage (V)");
    gDCR_Del_2_FBK_DARK_OV->GetYaxis()->SetTitle(title_DCR_Area);

    gDCR_Del_3_FBK_DARK_OV->SetMarkerStyle(22);
    gDCR_Del_3_FBK_DARK_OV->SetMarkerColor(kMagenta);
    gDCR_Del_3_FBK_DARK_OV->SetTitle();
    gDCR_Del_3_FBK_DARK_OV->GetXaxis()->SetTitle("Bias Voltage (V)");
    gDCR_Del_3_FBK_DARK_OV->GetYaxis()->SetTitle(title_DCR_Area);

    //------------------------------

    gCT_Del_1_FBK_DARK_OV->SetMarkerStyle(22);
    gCT_Del_1_FBK_DARK_OV->SetMarkerColor(kOrange+1);
    gCT_Del_1_FBK_DARK_OV->SetTitle();
    gCT_Del_1_FBK_DARK_OV->GetXaxis()->SetTitle("Bias Voltage (V)");
    gCT_Del_1_FBK_DARK_OV->GetYaxis()->SetTitle(title_P_CT);

    gCT_Del_2_FBK_DARK_OV->SetMarkerStyle(22);
    gCT_Del_2_FBK_DARK_OV->SetMarkerColor(kRed);
    gCT_Del_2_FBK_DARK_OV->SetTitle();
    gCT_Del_2_FBK_DARK_OV->GetXaxis()->SetTitle("Bias Voltage (V)");
    gCT_Del_2_FBK_DARK_OV->GetYaxis()->SetTitle(title_P_CT);

    gCT_Del_3_FBK_DARK_OV->SetMarkerStyle(22);
    gCT_Del_3_FBK_DARK_OV->SetMarkerColor(kMagenta);
    gCT_Del_3_FBK_DARK_OV->SetTitle();
    gCT_Del_3_FBK_DARK_OV->GetXaxis()->SetTitle("Bias Voltage (V)");
    gCT_Del_3_FBK_DARK_OV->GetYaxis()->SetTitle(title_P_CT);

    //------------------------------


    gCT_1_FBK_LASER->SetMarkerStyle(21);
    gCT_1_FBK_LASER->SetMarkerColor(kOrange+2);
    gCT_1_FBK_LASER->SetTitle();
    gCT_1_FBK_LASER->GetXaxis()->SetTitle("Bias Voltage (V)");
    gCT_1_FBK_LASER->GetYaxis()->SetTitle(title_P_CT);

    //------------------------------

    gDCR_1_SensL_OV->SetMarkerStyle(20);
    gDCR_1_SensL_OV->SetMarkerColor(kGreen+1);
    gDCR_1_SensL_OV->SetTitle();
    gDCR_1_SensL_OV->GetXaxis()->SetTitle("Overvoltage (V)");
    gDCR_1_SensL_OV->GetYaxis()->SetTitle(title_DCR_Area);

    gDCR_2_SensL_OV->SetMarkerStyle(20);
    gDCR_2_SensL_OV->SetMarkerColor(kBlue);
    gDCR_2_SensL_OV->SetTitle();
    gDCR_2_SensL_OV->GetXaxis()->SetTitle("Overvoltage (V)");
    gDCR_2_SensL_OV->GetYaxis()->SetTitle(title_DCR_Area);


    //------------------------------

    gCT_1_SensL_OV->SetMarkerStyle(20);
    gCT_1_SensL_OV->SetMarkerColor(kGreen+1);
    gCT_1_SensL_OV->SetTitle();
    gCT_1_SensL_OV->GetXaxis()->SetTitle("Overvoltage (V)");
    gCT_1_SensL_OV->GetYaxis()->SetTitle(title_P_CT);

    gCT_2_SensL_OV->SetMarkerStyle(20);
    gCT_2_SensL_OV->SetMarkerColor(kBlue);
    gCT_2_SensL_OV->SetTitle();
    gCT_2_SensL_OV->GetXaxis()->SetTitle("Overvoltage (V)");
    gCT_2_SensL_OV->GetYaxis()->SetTitle(title_P_CT);

    //------------------------------



    TCanvas *cV_CT_FBK_1 = new TCanvas("cV_CT_FBK_1", "cV_CT_FBK_1",w,h);
    cV_CT_FBK_1->SetGrid();
    auto legendCT_FBK_1 = new TLegend(0.15,0.70,0.45,0.85);
    legendCT_FBK_1->AddEntry(gCT_1_FBK_DARK,    "HD3-2 (1), from CNT (DARK)","p");
    legendCT_FBK_1->AddEntry(gCT_Del_1_FBK_DARK,"HD3-2 (1), from DELAYS (DARK)","p");
    legendCT_FBK_1->AddEntry(gCT_1_FBK_LASER,     "HD3-2 (1), from LASER","p");
    legendCT_FBK_1->Draw();
    TMultiGraph *mgV_CT_FBK_1 = new TMultiGraph("mgV_CT_FBK_1", title_P_CT_mg);
    mgV_CT_FBK_1->Add(gCT_1_FBK_DARK);
    mgV_CT_FBK_1->Add(gCT_Del_1_FBK_DARK);
    mgV_CT_FBK_1->Add(gCT_1_FBK_LASER);
    mgV_CT_FBK_1->Draw("AP");
    legendCT_FBK_1->Draw();

    //------------------------------

    TCanvas *cV_DCR_FBK_SensL_1 = new TCanvas("cV_DCR_FBK_SensL_1", "cV_DCR_FBK_SensL_1",w,h);
    cV_DCR_FBK_SensL_1->SetGrid();
    auto legendDCR_FBK_SensL = new TLegend(0.15,0.70,0.4,0.9);
    legendDCR_FBK_SensL->AddEntry(gDCR_1_FBK_DARK_OV, "FBK NUV HD3-2 (1)","p");
    legendDCR_FBK_SensL->AddEntry(gDCR_2_FBK_DARK_OV, "FBK NUV HD3-2 (2)","p");
    legendDCR_FBK_SensL->AddEntry(gDCR_3_FBK_DARK_OV, "FBK NUV HD3-2 (3)","p");
    legendDCR_FBK_SensL->AddEntry(gDCR_1_SensL_OV,    "SensL MicroFJ (1)","p");
    legendDCR_FBK_SensL->AddEntry(gDCR_2_SensL_OV,    "SensL MicroFJ (2)","p");
    TMultiGraph *mgV_DCR_FBK_SensL_1 = new TMultiGraph("mgV_DCR_FBK_SensL_1", title_DCR_Area_mg);
    mgV_DCR_FBK_SensL_1->Add(gDCR_1_FBK_DARK_OV);
    mgV_DCR_FBK_SensL_1->Add(gDCR_2_FBK_DARK_OV);
    mgV_DCR_FBK_SensL_1->Add(gDCR_3_FBK_DARK_OV);
    mgV_DCR_FBK_SensL_1->Add(gDCR_1_SensL_OV);
    mgV_DCR_FBK_SensL_1->Add(gDCR_2_SensL_OV);
    mgV_DCR_FBK_SensL_1->Draw("AP");
    legendDCR_FBK_SensL->Draw();


    //------------------------------

    TCanvas *cV_CT_FBK_SensL_1 = new TCanvas("cV_CT_FBK_SensL_1", "cV_CT_FBK_SensL_1",w,h);
    cV_CT_FBK_SensL_1->SetGrid();
    auto legendCT_FBK_SensL = new TLegend(0.15,0.70,0.4,0.9);
    legendCT_FBK_SensL->AddEntry(gCT_1_FBK_DARK_OV, "FBK NUV HD3-2 (1)","p");
    legendCT_FBK_SensL->AddEntry(gCT_2_FBK_DARK_OV, "FBK NUV HD3-2 (2)","p");
    legendCT_FBK_SensL->AddEntry(gCT_3_FBK_DARK_OV, "FBK NUV HD3-2 (3)","p");
    legendCT_FBK_SensL->AddEntry(gCT_1_SensL_OV,    "SensL MicroFJ (1)","p");
    legendCT_FBK_SensL->AddEntry(gCT_2_SensL_OV,    "SensL MicroFJ (2)","p");
    TMultiGraph *mgV_CT_FBK_SensL_1 = new TMultiGraph("mgV_CT_FBK_SensL_1", title_P_CT_mg);
    mgV_CT_FBK_SensL_1->Add(gCT_1_FBK_DARK_OV);
    mgV_CT_FBK_SensL_1->Add(gCT_2_FBK_DARK_OV);
    mgV_CT_FBK_SensL_1->Add(gCT_3_FBK_DARK_OV);
    mgV_CT_FBK_SensL_1->Add(gCT_1_SensL_OV);
    mgV_CT_FBK_SensL_1->Add(gCT_2_SensL_OV);
    mgV_CT_FBK_SensL_1->Draw("AP");
    legendCT_FBK_SensL->Draw();








}
