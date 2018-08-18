/******************************************************************************\
 *  MEASURES FROM PCB OF INTERFACE BOARD, MADE BY SILVANO GALLIAN
 *  FILTER 3 STAGES (92dB @ 100kHz)
 *
 *  MODEL based on
 *              > 20180608_SiPM_BOARD_10.pdf
 *
 *  SCHEMATIC:
 *               > NUMI AL-3T (OUT = 24.00V (*)) + Dummy Back Plane + Dragon + SCB + SW_Silvano + Filtro 3s + ch1 + R_LOAD
 *               > PCB (SW_Silvano + Filtro 3s + ch1 + R_LOAD) connected ch4 of SCB
 *               > Dummy Back Plane connected with Ethernet cable to PC
 *               > Ch enabled by PC (using ApplHV IP HV)
 *               > CONNECTION WITH SCB (as reported in Progetto_SiPM_BOARD_SW_1ch_StoD_02.txt)
 *                 # Connector:
 *                     20180622_Connettore_SCB.
 *                 # +3.3V
 *                     Connector: pin 2 and 4 (5.189 V (*))
 *                 # -3.3V
 *                     Connector: pin 37 (-2.955 V (*))
 *                 # +6V
 *                     Connector: pin 6 and 8 (5.941 V (*))
 *                 # HV_set (HV_set)
 *                     Connector: pin 38
 *                 # ENABLE
 *                     Connector: pin 10 (with the R = 10 kOhm as reported in 20180622_Connettore_SCB)
 *
 *  MEASURES TAKEN WITH
 *              > Fluke 179 True RMS Multimeter (also for the (*) above)
 *
 *  LOADs USED (R_LOAD):
 *              > MAXi      -> 39  kOhm
 *              > MINi      -> 220 kOhm (code: 224)
 *              > SuperMAXi -> 22  kOhm (code: 223)
 *
 *
 *
\******************************************************************************/

#include "TMath.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TF1.h"
#include "TH1D.h"
#include "TMultiGraph.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TLegend.h"

#define n_measured 15
#define n_used 6

#define h 600
#define w 1000


void fit_linear(TGraphErrors *g, TCanvas *c1);
void fit_linear_range(TGraphErrors *g, TCanvas *c1, double low, double high);
void Misure_PCB_DRAGON_20180711();

double m_V_out, errm_V_out, q_V_out, errq_V_out;
double m_V_out_V_vmon, errm_V_out_V_vmon, q_V_out_V_vmon, errq_V_out_V_vmon;

bool fit_V_out_bool = false;
bool fit_V_out_V_vmon_bool = false;

void Misure_PCB_DRAGON_20180712(){
    // HV_set (NOT TRIMMER as before); the voltage of the pin 38 is about 30 - 35 mV more than HV_set[i]
    // from Terminal: $ ApplyHV 192.168.1.5 HV_set[i]
    double HV_set[] =                {100,    200,     300,     400,     500,     600,    700,    800,     900,    1000,   1100,     1200,   1300,   1400,   1500};

    // LOAD: INF
    double V_out_SuperMINi[] =       {2.615,  5.241,   7.88,    10.48,   13.10,   15.72,  18.35,  20.98,   23.60,  26.23,  28.86,    31.49,  34.11,  36.74,  39.37};
    double V_vmon_SuperMINi[] =      {-0.063, -0.159,  -0.256,  -0.350,  -0.446,  -0.543, -0.639, -0.734,  -0.831, -0.927, -1.023,   -1.120, -1.215, -1.312, -1.408};
    double V_imon_SuperMINi[] =      {-0.00,  -0.00,    -0.00,  -0.03,   -0.03,   -0.03,  -0.04,   -0.04,   -0.04,  -0.05,  -0.05,    -0.05,  -0.05,  -0.06,  -0.06};

    // LOAD: 224 (220 kOhm)
    double V_out_MINi[] =            {2.615,  5.241,   7.88,    10.48,   13.10,   15.72,  18.35,  20.98,   23.60,  26.23,  28.86,    31.48,  34.11,  36.74,  39.36};
    double V_vmon_MINi[] =           {-0.063, -0.159,  -0.255,  -0.350,  -0.446,  -0.542, -0.638, -0.734,  -0.830, -0.926, -1.023,   -1.119, -1.215, -1.311, -1.407};
    double V_imon_MINi[] =           {-0.00,  -0.01,   -0.02,   -0.05,   -0.06,   -0.07,  -0.08,  -0.09,   -0.10,  -0.11,  -0.12,    -0.13,  -0.13,  -0.14,  -0.15};

    // LOAD: 39 kOhm
    double V_out_MAXi[] =            {2.615,  5.240,  7.88,     10.47,   13.09,   15.72,  18.34,  20.97,    23.60,  26.23,  28.85,   31.48,  34.11,  36.74,  39.36};
    double V_vmon_MAXi[] =           {-0.063, -0.159, -0.255,   -0.350,  -0.446,  -0.542, -0.638, -0.734,   -0.830, -0.925, -1.022,  -1.118, -1.214, -1.311, -1.406};
    double V_imon_MAXi[] =           {-0.03,  -0.06,  -0.1,     -0.15,   -0.19,   -0.22,  -0.25,  -0.29,    -0.32,  -0.36,  -0.39,   -0.42,  -0.45,  -0.49,  -0.53};

    // LOAD: 223 (22 kOhm)
    double V_out_SuperMAXi[] =       {2.614,  5.239,   7.88,    10.47,   13.09,   15.72,  18.34,   20.97,   23.60,  26.22,  28.85,   31.47,  34.10,  36.72,  39.35};
    double V_vmon_SuperMAXi[] =      {-0.063, -0.159,  -0.255,  -0.349,  -0.445,  -0.541, -0.637,  -0.733,  -0.829, -0.925, -1.021,  -1.117, -1.213, -1.310, -1.405};
    double V_imon_SuperMAXi[] =      {-0.06,  -0.11,   -0.17,   -0.25,   -0.31,   -0.37,  -0.43,   -0.49,   -0.54,  -0.60,  -0.66,   -0.72,  -0.78,  -0.83,  -0.89};

    // errors: always the same:
    double errHV_set[n_measured], errV_out_MAXi[n_measured], errV_vmon_MAXi[n_measured], errV_imon_MAXi[n_measured], errV_out_MINi[n_measured], errV_vmon_MINi[n_measured], errV_imon_MINi[n_measured], errV_imon_MINi_h[n_measured], errV_out_SuperMAXi[n_measured], errV_vmon_SuperMAXi[n_measured], errV_imon_SuperMAXi[n_measured], errV_out_SuperMINi[n_measured], errV_vmon_SuperMINi[n_measured], errV_imon_SuperMINi[n_measured], errV_imon_SuperMINi_h[n_measured];
    for(int i=0; i<n_measured; i++){
      errHV_set[i] = 1;
      if(HV_set[i] <= 200) errV_out_MAXi[i] = errV_out_MINi[i] = errV_out_SuperMAXi[i] = errV_out_SuperMINi[i] = 0.001;
      else                 errV_out_MAXi[i] = errV_out_MINi[i] = errV_out_SuperMAXi[i] = errV_out_SuperMINi[i] = 0.01;
      errV_vmon_MAXi[i] = errV_vmon_MINi[i] = errV_vmon_SuperMAXi[i] = errV_vmon_SuperMINi[i] = 0.001;
      errV_imon_MAXi[i] = errV_imon_MINi[i] = errV_imon_SuperMAXi[i] = errV_imon_SuperMINi[i] = 0.01;
      if(HV_set[i] <= 300) errV_imon_MINi_h[i] = errV_imon_SuperMINi_h[i] = 0.00;
      else                 errV_imon_MINi_h[i] = errV_imon_SuperMINi_h[i] = 0.01;
    }

    //------------------------------

    TGraphErrors *gV_out_MAXi  = new TGraphErrors(n_measured, HV_set, V_out_MAXi, errHV_set, errV_out_MAXi );
    TGraphErrors *gV_vmon_MAXi = new TGraphErrors(n_measured, HV_set, V_vmon_MAXi,errHV_set, errV_vmon_MAXi);
    TGraphErrors *gV_imon_MAXi = new TGraphErrors(n_measured, HV_set, V_imon_MAXi,errHV_set, errV_imon_MAXi);

    TGraphErrors *gV_out_MINi  = new TGraphErrors(n_measured, HV_set, V_out_MINi, errHV_set, errV_out_MINi );
    TGraphErrors *gV_vmon_MINi = new TGraphErrors(n_measured, HV_set, V_vmon_MINi,errHV_set, errV_vmon_MINi);
    TGraphErrors *gV_imon_MINi = new TGraphErrors(n_measured, HV_set, V_imon_MINi,errHV_set, errV_imon_MINi);
    // TGraphAsymmErrors *gV_imon_MINi = new TGraphAsymmErrors(n_measured, HV_set, V_imon_MINi,errHV_set,errHV_set,errV_imon_MINi,errV_imon_MINi_h);

    TGraphErrors *gV_out_SuperMINi  = new TGraphErrors(n_measured, HV_set, V_out_SuperMINi, errHV_set, errV_out_SuperMINi );
    TGraphErrors *gV_vmon_SuperMINi = new TGraphErrors(n_measured, HV_set, V_vmon_SuperMINi,errHV_set, errV_vmon_SuperMINi);
    TGraphErrors *gV_imon_SuperMINi = new TGraphErrors(n_measured, HV_set, V_imon_SuperMINi,errHV_set, errV_imon_SuperMINi);
    // TGraphAsymmErrors *gV_imon_SuperMINi = new TGraphAsymmErrors(n_measured, HV_set, V_imon_SuperMINi,errHV_set,errHV_set,errV_imon_SuperMINi,errV_imon_SuperMINi_h);

    TGraphErrors *gV_out_SuperMAXi  = new TGraphErrors(n_measured, HV_set, V_out_SuperMAXi, errHV_set, errV_out_SuperMAXi );
    TGraphErrors *gV_vmon_SuperMAXi = new TGraphErrors(n_measured, HV_set, V_vmon_SuperMAXi,errHV_set, errV_vmon_SuperMAXi);
    TGraphErrors *gV_imon_SuperMAXi = new TGraphErrors(n_measured, HV_set, V_imon_SuperMAXi,errHV_set, errV_imon_SuperMAXi);

    //------------------------------

    gV_out_MAXi->SetMarkerStyle(20);
    gV_vmon_MAXi->SetMarkerStyle(20);
    gV_imon_MAXi->SetMarkerStyle(20);

    gV_out_MINi->SetMarkerStyle(20);
    gV_vmon_MINi->SetMarkerStyle(20);
    gV_imon_MINi->SetMarkerStyle(20);

    gV_out_SuperMAXi->SetMarkerStyle(20);
    gV_vmon_SuperMAXi->SetMarkerStyle(20);
    gV_imon_SuperMAXi->SetMarkerStyle(20);

    gV_out_SuperMINi->SetMarkerStyle(20);
    gV_vmon_SuperMINi->SetMarkerStyle(20);
    gV_imon_SuperMINi->SetMarkerStyle(20);

    // gV_out_MAXi->SetMarkerSize(3);
    // gV_vmon_MAXi->SetMarkerSize(3);
    // gV_imon_MAXi->SetMarkerSize(3);
    //
    // gV_out_MINi->SetMarkerSize(3);
    // gV_vmon_MINi->SetMarkerSize(3);
    // gV_imon_MINi->SetMarkerSize(3);

    int color_SuperMINi = kGreen+1;
    int color_MINi      = kBlue;
    int color_MAXi      = kOrange;
    int color_SuperMAXi = kRed;

    gV_out_SuperMINi->SetMarkerColor(color_SuperMINi);
    gV_vmon_SuperMINi->SetMarkerColor(color_SuperMINi);
    gV_imon_SuperMINi->SetMarkerColor(color_SuperMINi);

    gV_out_MINi->SetMarkerColor(color_MINi);
    gV_vmon_MINi->SetMarkerColor(color_MINi);
    gV_imon_MINi->SetMarkerColor(color_MINi);

    gV_out_MAXi->SetMarkerColor(color_MAXi);
    gV_vmon_MAXi->SetMarkerColor(color_MAXi);
    gV_imon_MAXi->SetMarkerColor(color_MAXi);


    gV_out_SuperMAXi->SetMarkerColor(color_SuperMAXi);
    gV_vmon_SuperMAXi->SetMarkerColor(color_SuperMAXi);
    gV_imon_SuperMAXi->SetMarkerColor(color_SuperMAXi);

    gV_out_MAXi->SetTitle();
    gV_vmon_MAXi->SetTitle();
    gV_imon_MAXi->SetTitle();

    gV_out_MINi->SetTitle();
    gV_vmon_MINi->SetTitle();
    gV_imon_MINi->SetTitle();

    gV_out_SuperMAXi->SetTitle();
    gV_vmon_SuperMAXi->SetTitle();
    gV_imon_SuperMAXi->SetTitle();

    gV_out_SuperMINi->SetTitle();
    gV_vmon_SuperMINi->SetTitle();
    gV_imon_SuperMINi->SetTitle();


    //------------------------------
    //----------- GRAPHS -----------
    //------------------------------
    //
    // TCanvas *cV_out_MAXi = new TCanvas("cV_out_MAXi", "cV_out_MAXi");
    // cV_out_MAXi->SetGrid();
    // gV_out_MAXi->Draw("AP");
    //
    // TCanvas *cV_vmon_MAXi = new TCanvas("cV_vmon_MAXi", "cV_vmon_MAXi");
    // cV_vmon_MAXi->SetGrid();
    // gV_vmon_MAXi->Draw("AP");
    //
    // TCanvas *cV_imon_MAXi = new TCanvas("cV_imon_MAXi", "cV_imon_MAXi");
    // cV_imon_MAXi->SetGrid();
    // gV_imon_MAXi->Draw("AP");
    //
    //
    // TCanvas *cV_out_MINi = new TCanvas("cV_out_MINi", "cV_out_MINi");
    // cV_out_MINi->SetGrid();
    // gV_out_MINi->Draw("AP");
    //
    // TCanvas *cV_vmon_MINi = new TCanvas("cV_vmon_MINi", "cV_vmon_MINi");
    // cV_vmon_MINi->SetGrid();
    // gV_vmon_MINi->Draw("AP");
    //
    // TCanvas *cV_imon_MINi = new TCanvas("cV_imon_MINi", "cV_imon_MINi");
    // cV_imon_MINi->SetGrid();
    // gV_imon_MINi->Draw("AP");
    //
    //
    // TCanvas *cV_out_SuperMINi = new TCanvas("cV_out_SuperMINi", "cV_out_SuperMINi");
    // cV_out_SuperMINi->SetGrid();
    // gV_out_SuperMINi->Draw("AP");
    //
    // TCanvas *cV_vmon_SuperMINi = new TCanvas("cV_vmon_SuperMINi", "cV_vmon_SuperMINi");
    // cV_vmon_SuperMINi->SetGrid();
    // gV_vmon_SuperMINi->Draw("AP");
    //
    // TCanvas *cV_imon_SuperMINi = new TCanvas("cV_imon_SuperMINi", "cV_imon_SuperMINi");
    // cV_imon_SuperMINi->SetGrid();
    // gV_imon_SuperMINi->Draw("AP");
    //
    //
    // TCanvas *cV_out_SuperMAXi = new TCanvas("cV_out_SuperMAXi", "cV_out_SuperMAXi");
    // cV_out_SuperMAXi->SetGrid();
    // gV_out_SuperMAXi->Draw("AP");
    //
    // TCanvas *cV_vmon_SuperMAXi = new TCanvas("cV_vmon_SuperMAXi", "cV_vmon_SuperMAXi");
    // cV_vmon_SuperMAXi->SetGrid();
    // gV_vmon_SuperMAXi->Draw("AP");
    //
    // TCanvas *cV_imon_SuperMAXi = new TCanvas("cV_imon_SuperMAXi", "cV_imon_SuperMAXi");
    // cV_imon_SuperMAXi->SetGrid();
    // gV_imon_SuperMAXi->Draw("AP");


    //------------------------------
    //-------- MULTIGRAPHs ---------
    //------------------------------

    TMultiGraph *V_out = new TMultiGraph("V_out", ";DAC_{ApplyHV} (adim); V_{HV} (V)");
    V_out->Add(gV_out_MINi);
    V_out->Add(gV_out_MAXi);
    V_out->Add(gV_out_SuperMAXi);
    V_out->Add(gV_out_SuperMINi);


    TMultiGraph *V_vmon = new TMultiGraph("V_vmon", ";DAC_{ApplyHV} (adim); V_{MON} (V)");
    V_vmon->Add(gV_vmon_MINi);
    V_vmon->Add(gV_vmon_MAXi);
    V_vmon->Add(gV_vmon_SuperMAXi);
    V_vmon->Add(gV_vmon_SuperMINi);



    TMultiGraph *V_imon = new TMultiGraph("V_imon", ";DAC_{ApplyHV} (adim); V_{I_{MON}} (V)");
    V_imon->Add(gV_imon_MINi);
    V_imon->Add(gV_imon_MAXi);
    V_imon->Add(gV_imon_SuperMAXi);
    V_imon->Add(gV_imon_SuperMINi);

    //------------------------------
    // Since I always use the same colors:
    auto legend_TL = new TLegend(0.15,0.725,0.325,0.9);
    legend_TL->AddEntry(gV_out_SuperMINi,"R_{LOAD} \\rightarrow inf", "p");
    legend_TL->AddEntry(gV_out_MINi,"R_{LOAD} \\simeq 220k\\Omega", "p");
    legend_TL->AddEntry(gV_out_MAXi,"R_{LOAD} \\simeq 39k\\Omega","p");
    legend_TL->AddEntry(gV_out_SuperMAXi,"R_{LOAD} \\simeq 22k\\Omega", "p");

    auto legend_TR = new TLegend(0.725,0.725,0.9,0.9);
    legend_TR->AddEntry(gV_out_SuperMINi,"R_{LOAD} \\rightarrow inf", "p");
    legend_TR->AddEntry(gV_out_MINi,"R_{LOAD} \\simeq 220k\\Omega", "p");
    legend_TR->AddEntry(gV_out_MAXi,"R_{LOAD} \\simeq 39k\\Omega","p");
    legend_TR->AddEntry(gV_out_SuperMAXi,"R_{LOAD} \\simeq 22k\\Omega", "p");

    auto legend_BL = new TLegend(0.15,0.15,0.325,0.325);
    legend_BL->AddEntry(gV_out_SuperMINi,"R_{LOAD} \\rightarrow inf", "p");
    legend_BL->AddEntry(gV_out_MINi,"R_{LOAD} \\simeq 220k\\Omega", "p");
    legend_BL->AddEntry(gV_out_MAXi,"R_{LOAD} \\simeq 39k\\Omega","p");
    legend_BL->AddEntry(gV_out_SuperMAXi,"R_{LOAD} \\simeq 22k\\Omega", "p");




    //------------------------------

    TCanvas *cV_out = new TCanvas("cV_out", "cV_out", w, h);
    cV_out->SetGrid();
    V_out->Draw("AP");
    auto legend_V_out = new TLegend(0.15,0.7,0.35,0.9);
    legend_V_out->AddEntry(gV_out_SuperMINi,"V_{HV},R_{LOAD} \\rightarrow inf", "p");
    legend_V_out->AddEntry(gV_out_MINi,"V_{HV},R_{LOAD} \\simeq 220k\\Omega", "p");
    legend_V_out->AddEntry(gV_out_MAXi,"V_{HV},R_{LOAD} \\simeq 39k\\Omega","p");
    legend_V_out->AddEntry(gV_out_SuperMAXi,"V_{HV},R_{LOAD} \\simeq 22k\\Omega", "p");
    // legend_V_out->Draw();
    legend_TL->Draw();

    TCanvas *cV_vmon = new TCanvas("cV_vmon", "cV_vmon", w, h);
    cV_vmon->SetGrid();
    V_vmon ->Draw("AP");
    auto legend_V_vmon = new TLegend(0.7,0.7,0.9,0.9);
    legend_V_vmon->AddEntry(gV_vmon_SuperMINi,"V_{MON},R_{LOAD} \\rightarrow inf","p");
    legend_V_vmon->AddEntry(gV_vmon_MINi,"V_{MON},R_{LOAD} \\simeq 220k\\Omega", "p");
    legend_V_vmon->AddEntry(gV_vmon_MAXi,"V_{MON},R_{LOAD} \\simeq 39k\\Omega","p");
    legend_V_vmon->AddEntry(gV_vmon_SuperMAXi,"V_{MON},R_{LOAD} \\simeq 22k\\Omega","p");
    // legend_V_vmon->Draw();
    legend_TR->Draw();

    TCanvas *cV_imon = new TCanvas("cV_imon", "cV_imon", w, h);
    cV_imon->SetGrid();
    V_imon ->Draw("AP");
    auto legend_V_imon = new TLegend(0.15,0.15,0.35,0.35);
    legend_V_imon->AddEntry(gV_imon_SuperMINi,"V_{I_{MON}},R_{LOAD} \\rightarrow inf","p");
    legend_V_imon->AddEntry(gV_imon_MINi,"V_{I_{MON}},R_{LOAD} \\simeq 220k\\Omega", "p");
    legend_V_imon->AddEntry(gV_imon_MAXi,"V_{I_{MON}},R_{LOAD} \\simeq 39k\\Omega","p");
    legend_V_imon->AddEntry(gV_imon_SuperMAXi,"V_{I_{MON}},R_{LOAD} \\simeq 22k\\Omega","p");
    // legend_V_imon->Draw();
    legend_BL->Draw();


    //------------------------------
    //------ TEST_Z and cout -------
    //------------------------------
    /*
    double TEST_Z;
    bool OK;

    // V_out
    cout<<endl;
    cout<<"/////// ONLY FOR MINi and MAXi ///////"<<endl;
    cout<<"HV_set\t\t | V_out, I_min\t   | V_out, I_max     | TEST_Z"<<endl;
    cout<<"\t\t |\t\t   |\t\t      |"<<endl;
    for(int i=0; i<n; i++){
        TEST_Z = TMath::Abs(V_out_MINi[i] - V_out_MAXi[i]) / TMath::Sqrt( errV_out_MINi[i]*errV_out_MINi[i] + errV_out_MAXi[i]*errV_out_MAXi[i]);
        OK = TEST_Z<1.96;
        printf("%.3lf +- %.3lf\t | %.2lf +- %.2lf   | %.2lf +- %.2lf    | %.3lf\t", HV_set[i], errHV_set[i], V_out_MINi[i], errV_out_MINi[i], V_out_MAXi[i], errV_out_MAXi[i], TEST_Z);
        if(OK==true) cout<<"true"<<endl;
        else         cout<<"false"<<endl;
    }

    // V_vmon
    cout<<endl;
    cout<<"HV_set\t\t | V_vmon, I_min   | V_vmon, I_max    | TEST_Z"<<endl;
    cout<<"\t\t |\t\t   |\t\t      |"<<endl;
    for(int i=0; i<n; i++){
        TEST_Z = TMath::Abs(V_vmon_MINi[i] - V_vmon_MAXi[i]) / TMath::Sqrt( errV_vmon_MINi[i]*errV_vmon_MINi[i] + errV_vmon_MAXi[i]*errV_vmon_MAXi[i]);
        OK = TEST_Z<1.96;
        printf("%.3lf +- %.3lf\t | %.3lf +- %.3lf | %.3lf +- %.3lf  | %.3lf\t", HV_set[i], errHV_set[i], V_vmon_MINi[i], errV_vmon_MINi[i], V_vmon_MAXi[i], errV_vmon_MAXi[i], TEST_Z);
        if(OK==true) cout<<"true"<<endl;
        else         cout<<"false"<<endl;
    }

    // V_imon
    cout<<endl;
    cout<<"HV_set\t\t | V_imon, I_min      | V_imon, I_max"<<endl;
    cout<<"\t\t |\t\t   |\t\t"<<endl;
    for(int i=0; i<n; i++){
        printf("%.3lf +- %.3lf\t | %.2lf +- %.2lf   | %.2lf +- %.2lf\n", HV_set[i], errHV_set[i], V_imon_MINi[i], errV_imon_MINi[i], V_imon_MAXi[i], errV_imon_MAXi[i]);
    }
    */

    //------------------------------
    //------------ FIT -------------
    //------------------------------

    //----------- V_out ------------
    m_V_out = errm_V_out = q_V_out = errq_V_out = 0.;
    int n_mean;

    cout<<endl;
    cout<<"///// V_out /////"<<endl;
    n_mean = 4; // fit_linear is called n_mean times
    fit_V_out_bool = true;
    fit_V_out_V_vmon_bool = false;
    fit_linear(gV_out_SuperMINi, cV_out);
    fit_linear(gV_out_MINi, cV_out);
    fit_linear(gV_out_MAXi, cV_out);
    fit_linear(gV_out_SuperMAXi, cV_out);

    // mean values
    q_V_out /= n_mean;
    m_V_out /= n_mean;
    errq_V_out = TMath::Sqrt(errq_V_out) / n_mean;
    errm_V_out = TMath::Sqrt(errm_V_out) / n_mean;
    cout<<"Mean"<<endl;
    printf("V_out = (%lf +- %lf) * HV_set + (%lf +- %lf)\n",  m_V_out, errm_V_out, q_V_out, errq_V_out);

    // err on V_out
    // Less important, I can use the V_vmon to adjust the voltage
    int n_est = 1500;
    double V_out_est[n_est], errV_out_est[n_est], HV_set_est[n_est], errHV_set_est[n_est];
    double HV_set_est_min, HV_set_est_max, binw;
    HV_set_est_min = 0;
    HV_set_est_max = 1500;
    binw = (HV_set_est_max-HV_set_est_min)/n_est;
    for(int i=0; i<n_est; i++){
      HV_set_est[i] = HV_set_est_min + i * binw;
      errHV_set_est[i] = 1;
    }
    for(int i=0; i<n_est; i++){
      V_out_est[i] = m_V_out * HV_set_est[i] + q_V_out;
      errV_out_est[i] = TMath::Sqrt( HV_set_est[i]*HV_set_est[i]*errm_V_out*errm_V_out + m_V_out*m_V_out*errHV_set_est[i]*errHV_set_est[i] + errq_V_out*errq_V_out );
      errV_out_est[i] *= 1000; // in mV
    }
    // TCanvas *cerrV_out_est = new TCanvas("cerrV_out_est", "cerrV_out_est", w, h);
    TGraph* gerrV_out_est = new TGraph(n_est,V_out_est,errV_out_est);
    gerrV_out_est->SetTitle();
    gerrV_out_est->GetXaxis()->SetTitle("V_out [V]");
    gerrV_out_est->GetYaxis()->SetTitle("errV_out [mV]");
    // gerrV_out_est->Draw("AL");



    //------------------------------
    //---------- CONSTANT ----------
    //------------------------------
    double min_K = -40;
    double max_K = -20;
    double binw_K = 0.1;
    int nbins_K;

    nbins_K = (int)((max_K-min_K)/binw_K);

    TH1D* ptrhistK           = new TH1D("histK", "", nbins_K, min_K, max_K);
    TH1D* ptrhistK_SuperMINi = new TH1D("histK_SuperMINi", "", nbins_K, min_K, max_K);
    TH1D* ptrhistK_MINi      = new TH1D("histK_MINi", "", nbins_K, min_K, max_K);
    TH1D* ptrhistK_MAXi      = new TH1D("histK_MAXi", "", nbins_K, min_K, max_K);
    TH1D* ptrhistK_SuperMAXi = new TH1D("histK_SuperMAXi", "", nbins_K, min_K, max_K);


    //-----------------------------------
    //---------- FROM 26 TO 40 ----------
    //-----------------------------------
    double K_SuperMINi[n_used], K_MINi[n_used], K_MAXi[n_used], K_SuperMAXi[n_used];
    double V_out_SuperMINi_used[n_used], V_out_MINi_used[n_used], V_out_MAXi_used[n_used], V_out_SuperMAXi_used[n_used];
    double errV_out_SuperMINi_used[n_used], errV_out_MINi_used[n_used], errV_out_MAXi_used[n_used], errV_out_SuperMAXi_used[n_used];

    int index;
    int start = n_measured-n_used;

    for(int i=start; i<n_measured; i++){

      index=i-start;

      K_SuperMINi[index] = V_out_SuperMINi[i]/V_vmon_SuperMINi[i];
      K_MINi[index]      = V_out_MINi[i]/V_vmon_MINi[i];
      K_MAXi[index]      = V_out_MAXi[i]/V_vmon_MAXi[i];
      K_SuperMAXi[index] = V_out_SuperMAXi[i]*TMath::Power(V_vmon_SuperMAXi[i],-1);

      V_out_SuperMINi_used[index] = V_out_SuperMINi[i];
      V_out_MINi_used[index] = V_out_MINi[i];
      V_out_MAXi_used[index] = V_out_MAXi[i];
      V_out_SuperMAXi_used[index] = V_out_SuperMAXi[i];

      errV_out_SuperMINi_used[index] = errV_out_SuperMINi[i];
      errV_out_MINi_used[index] = errV_out_MINi[i];
      errV_out_MAXi_used[index] = errV_out_MAXi[i];
      errV_out_SuperMAXi_used[index] = errV_out_SuperMAXi[i];



      cout<<index<<"\t"<<K_SuperMINi[index]<<" \t"<<K_MINi[index]<<" \t"<<K_MAXi[index]<<" \t"<<K_SuperMAXi[index]<<endl;

      ptrhistK->Fill(K_SuperMINi[index]);
      ptrhistK->Fill(K_MINi[index]);
      ptrhistK->Fill(K_MAXi[index]);
      ptrhistK->Fill(K_SuperMAXi[index]);

      ptrhistK_SuperMINi->Fill(K_SuperMINi[index]);
      ptrhistK_MINi->Fill(K_MINi[index]);
      ptrhistK_MAXi->Fill(K_MAXi[index]);
      ptrhistK_SuperMAXi->Fill(K_SuperMAXi[index]);

    }


    TCanvas *chistK = new TCanvas("chistK", "chistK", w, h);

    ptrhistK->SetLineColor(kBlack);
    ptrhistK_SuperMINi->SetLineColor(kGreen);
    ptrhistK_MINi->SetLineColor(kBlue);
    ptrhistK_MAXi->SetLineColor(kOrange);
    ptrhistK_SuperMAXi->SetLineColor(kRed);

    ptrhistK_SuperMINi->SetFillColorAlpha(kGreen, 0.15);
    ptrhistK_MINi->SetFillColorAlpha(kBlue, 0.15);
    ptrhistK_MAXi->SetFillColorAlpha(kOrange, 0.15);
    ptrhistK_SuperMAXi->SetFillColorAlpha(kRed, 0.15);

    ptrhistK->Draw();
    ptrhistK_SuperMINi->Draw("same");
    ptrhistK_MINi->Draw("same");
    ptrhistK_MAXi->Draw("same");
    ptrhistK_SuperMAXi->Draw("same");

    TCanvas *cK = new TCanvas("cK", "cK", w, h);

    TGraphErrors *gK_SuperMINi  = new TGraphErrors(n_used, V_out_SuperMINi_used, K_SuperMINi, errV_out_SuperMINi_used, 0);
    TGraphErrors *gK_MINi  = new TGraphErrors(n_used, V_out_MINi_used, K_MINi, errV_out_MINi_used, 0);
    TGraphErrors *gK_MAXi  = new TGraphErrors(n_used, V_out_MAXi_used, K_MAXi, errV_out_MAXi_used, 0);
    TGraphErrors *gK_SuperMAXi  = new TGraphErrors(n_used, V_out_SuperMAXi_used, K_SuperMAXi, errV_out_SuperMAXi_used, 0);


    gK_SuperMINi->SetMarkerSize(3);
    gK_MINi->SetMarkerSize(3);
    gK_MAXi->SetMarkerSize(3);
    gK_SuperMAXi->SetMarkerSize(3);

    gK_SuperMINi->SetMarkerColor(kOrange-1);
    gK_MINi->SetMarkerColor(kOrange);
    gK_MAXi->SetMarkerColor(kOrange+1);
    gK_SuperMAXi->SetMarkerColor(kOrange+2);

    gK_SuperMINi->SetLineColor(kOrange-1);
    gK_MINi->SetLineColor(kOrange);
    gK_MAXi->SetLineColor(kOrange+1);
    gK_SuperMAXi->SetLineColor(kOrange+2);


    TMultiGraph *K_graph = new TMultiGraph("V_K", ";Vout (V); Vout/Vmon (adim)");
    K_graph->Add(gK_SuperMINi);
    K_graph->Add(gK_MINi);
    K_graph->Add(gK_MAXi);
    K_graph->Add(gK_SuperMAXi);

    K_graph->Draw("APL");




    //-------------------------
    //---------- ALL ----------
    //-------------------------


    double K_SuperMINi_all[n_measured], K_MINi_all[n_measured ], K_MAXi_all[n_measured ], K_SuperMAXi_all[n_measured];

    TH1D* ptrhistK_all           = new TH1D("histK_all", "", nbins_K, min_K, max_K);
    TH1D* ptrhistK_SuperMINi_all = new TH1D("histK_SuperMINi_all", "", nbins_K, min_K, max_K);
    TH1D* ptrhistK_MINi_all      = new TH1D("histK_MINi_all", "", nbins_K, min_K, max_K);
    TH1D* ptrhistK_MAXi_all      = new TH1D("histK_MAXi_all", "", nbins_K, min_K, max_K);
    TH1D* ptrhistK_SuperMAXi_all = new TH1D("histK_SuperMAXi_all", "", nbins_K, min_K, max_K);

    for(int i=0; i<n_measured; i++){


      K_SuperMINi_all[i] = V_out_SuperMINi[i]/V_vmon_SuperMINi[i];
      K_MINi_all[i]      = V_out_MINi[i]/V_vmon_MINi[i];
      K_MAXi_all[i]      = V_out_MAXi[i]/V_vmon_MAXi[i];
      K_SuperMAXi_all[i] = V_out_SuperMAXi[i]*TMath::Power(V_vmon_SuperMAXi[i],-1);

      cout<<i<<"\t"<<K_SuperMINi_all[i]<<" \t"<<K_MINi_all[i]<<" \t"<<K_MAXi_all[i]<<" \t"<<K_SuperMAXi_all[i]<<endl;

      ptrhistK_all->Fill(K_SuperMINi_all[i]);
      ptrhistK_all->Fill(K_MINi_all[i]);
      ptrhistK_all->Fill(K_MAXi_all[i]);
      ptrhistK_all->Fill(K_SuperMAXi_all[i]);

      ptrhistK_SuperMINi_all->Fill(K_SuperMINi_all[i]);
      ptrhistK_MINi_all->Fill(K_MINi_all[i]);
      ptrhistK_MAXi_all->Fill(K_MAXi_all[i]);
      ptrhistK_SuperMAXi_all->Fill(K_SuperMAXi_all[i]);

    }


    TCanvas *chistK_all = new TCanvas("chistK_all", "chistK_all", w, h);

    ptrhistK->SetLineColor(kBlack);
    ptrhistK_SuperMINi_all->SetLineColor(kGreen);
    ptrhistK_MINi_all->SetLineColor(kBlue);
    ptrhistK_MAXi_all->SetLineColor(kOrange);
    ptrhistK_SuperMAXi_all->SetLineColor(kRed);

    ptrhistK_SuperMINi_all->SetFillColorAlpha(kGreen, 0.15);
    ptrhistK_MINi_all->SetFillColorAlpha(kBlue, 0.15);
    ptrhistK_MAXi_all->SetFillColorAlpha(kOrange, 0.15);
    ptrhistK_SuperMAXi_all->SetFillColorAlpha(kRed, 0.15);

    ptrhistK->Draw();
    ptrhistK_SuperMINi_all->Draw("same");
    ptrhistK_MINi_all->Draw("same");
    ptrhistK_MAXi_all->Draw("same");
    ptrhistK_SuperMAXi_all->Draw("same");

    TCanvas *cK_all = new TCanvas("cK_all", "cK_all", w, h);

    TGraphErrors *gK_SuperMINi_all  = new TGraphErrors(n_measured, V_out_SuperMINi, K_SuperMINi_all, errV_out_SuperMINi, 0);
    TGraphErrors *gK_MINi_all  = new TGraphErrors(n_measured, V_out_MINi, K_MINi_all, errV_out_MINi, 0);
    TGraphErrors *gK_MAXi_all  = new TGraphErrors(n_measured, V_out_MAXi, K_MAXi_all, errV_out_MAXi, 0);
    TGraphErrors *gK_SuperMAXi_all  = new TGraphErrors(n_measured, V_out_SuperMAXi, K_SuperMAXi_all, errV_out_SuperMAXi, 0);


    gK_SuperMINi_all->SetMarkerSize(3);
    gK_MINi_all->SetMarkerSize(3);
    gK_MAXi_all->SetMarkerSize(3);
    gK_SuperMAXi_all->SetMarkerSize(3);

    gK_SuperMINi_all->SetMarkerColor(kOrange-1);
    gK_MINi_all->SetMarkerColor(kOrange);
    gK_MAXi_all->SetMarkerColor(kOrange+1);
    gK_SuperMAXi_all->SetMarkerColor(kOrange+2);

    gK_SuperMINi_all->SetLineColor(kOrange-1);
    gK_MINi_all->SetLineColor(kOrange);
    gK_MAXi_all->SetLineColor(kOrange+1);
    gK_SuperMAXi_all->SetLineColor(kOrange+2);


    TMultiGraph *K_graph_all = new TMultiGraph("V_K_all", ";Vout (V); Vout/Vmon (adim)");
    K_graph_all->Add(gK_SuperMINi_all);
    K_graph_all->Add(gK_MINi_all);
    K_graph_all->Add(gK_MAXi_all);
    K_graph_all->Add(gK_SuperMAXi_all);

    K_graph_all->Draw("APL");



    //-----------------------
    //--- V_out vs V_vmon ---
    //-----------------------


    TGraphErrors *gV_out_V_vmon_SuperMINi  = new TGraphErrors(n_measured, V_vmon_SuperMINi, V_out_SuperMINi, errV_vmon_SuperMINi, errV_out_SuperMINi );
    TGraphErrors *gV_out_V_vmon_MINi       = new TGraphErrors(n_measured, V_vmon_MINi, V_out_MINi, errV_vmon_MINi, errV_out_MINi );
    TGraphErrors *gV_out_V_vmon_MAXi       = new TGraphErrors(n_measured, V_vmon_MAXi, V_out_MAXi, errV_vmon_MAXi, errV_out_MAXi );
    TGraphErrors *gV_out_V_vmon_SuperMAXi  = new TGraphErrors(n_measured, V_vmon_SuperMAXi, V_out_SuperMAXi, errV_vmon_SuperMAXi, errV_out_SuperMAXi );

    gV_out_V_vmon_SuperMINi->SetMarkerStyle(20);
    gV_out_V_vmon_SuperMINi->SetMarkerColor(color_SuperMINi);
    gV_out_V_vmon_SuperMINi->SetTitle();

    gV_out_V_vmon_MINi->SetMarkerStyle(20);
    gV_out_V_vmon_MINi->SetMarkerColor(color_MINi);
    gV_out_V_vmon_MINi->SetTitle();

    gV_out_V_vmon_MAXi->SetMarkerStyle(20);
    gV_out_V_vmon_MAXi->SetMarkerColor(color_MAXi);
    gV_out_V_vmon_MAXi->SetTitle();

    gV_out_V_vmon_SuperMAXi->SetMarkerStyle(20);
    gV_out_V_vmon_SuperMAXi->SetMarkerColor(color_SuperMAXi);
    gV_out_V_vmon_SuperMAXi->SetTitle();

    TCanvas *cV_out_V_vmon = new TCanvas("cV_out_V_vmon", "cV_out_V_vmon", w, h);
    cV_out_V_vmon->SetGrid();

    TMultiGraph *mgV_out_V_vmon = new TMultiGraph("mgV_out_V_vmon", ";V_{MON} (V); V_{HV} (V)");
    mgV_out_V_vmon->Add(gV_out_V_vmon_SuperMAXi);
    mgV_out_V_vmon->Add(gV_out_V_vmon_MAXi);
    mgV_out_V_vmon->Add(gV_out_V_vmon_MINi);
    mgV_out_V_vmon->Add(gV_out_V_vmon_SuperMINi);

    mgV_out_V_vmon->Draw("AP");
    legend_TR->Draw();

    ////////
    // FIT
    ////////
    cout<<endl<<endl;
    q_V_out_V_vmon = 0;
    m_V_out_V_vmon = 0;
    errq_V_out_V_vmon = 0;
    errm_V_out_V_vmon = 0;
    n_mean = 4; // fit_linear is called n_mean times
    fit_V_out_V_vmon_bool = true;
    fit_V_out_bool = false;
    double low = -1.5;
    double high = -0.8;

    cout<<"R_LOAD = 22 k:"<<endl;
    fit_linear(gV_out_V_vmon_SuperMAXi,cV_out_V_vmon);
    cout<<"R_LOAD = 39 k:"<<endl;
    fit_linear(gV_out_V_vmon_MAXi,cV_out_V_vmon);
    cout<<"R_LOAD = 220 k:"<<endl;
    fit_linear(gV_out_V_vmon_MINi,cV_out_V_vmon);
    cout<<"R_LOAD = inf:"<<endl;
    fit_linear(gV_out_V_vmon_SuperMINi,cV_out_V_vmon);

    // cout<<"R_LOAD = 22 k:"<<endl;
    // fit_linear_range(gV_out_V_vmon_SuperMAXi,cV_out_V_vmon, low, high);
    // cout<<"R_LOAD = 39 k:"<<endl;
    // fit_linear_range(gV_out_V_vmon_MAXi,cV_out_V_vmon, low, high);
    // cout<<"R_LOAD = 220 k:"<<endl;
    // fit_linear_range(gV_out_V_vmon_MINi,cV_out_V_vmon, low, high);
    // cout<<"R_LOAD = inf:"<<endl;
    // fit_linear_range(gV_out_V_vmon_SuperMINi,cV_out_V_vmon, low, high);



    // mean values
    q_V_out_V_vmon /= n_mean;
    m_V_out_V_vmon /= n_mean;
    errq_V_out_V_vmon = TMath::Sqrt(errq_V_out_V_vmon) / n_mean;
    errm_V_out_V_vmon = TMath::Sqrt(errm_V_out_V_vmon) / n_mean;
    cout<<"Mean"<<endl;
    printf("V_out = (%lf +- %lf) * V_vmon + (%lf +- %lf)\n",  m_V_out_V_vmon, errm_V_out_V_vmon, q_V_out_V_vmon, errq_V_out_V_vmon);

    return;


}

void fit_linear(TGraphErrors *g, TCanvas *c1){
  // Fit
  TF1 *line = new TF1("line","[0]*x+[1]");
  // TF1 *line = new TF1("line","[0]*x");


  g->Fit("line", "q0");
  c1->cd();
  // g->Draw("same");
  // c2->cd();

  //g->Draw("same");

  // Get Parameters
  double m, q, errm, errq; // y = m x + q
  q = line->GetParameter(1);
  errq = line->GetParError(1);
  // q=0; errq=0;
  m = line->GetParameter(0);
  errm = line->GetParError(0);


  if(fit_V_out_bool){
      printf("V_out = (%lf +- %lf) * HV_set + (%lf +- %lf)\n",  m, errm, q, errq);
      // for the mean:
      q_V_out += q;
      m_V_out += m;

      errq_V_out += errq*errq;
      errm_V_out += errm*errm;
  }else{
      if(fit_V_out_V_vmon_bool){
          printf("V_out = (%.2lf +- %.2lf) * V_vmon + (%.2lf +- %.2lf)\n",  m, errm, q, errq);
          // for the mean:
          q_V_out_V_vmon += q;
          m_V_out_V_vmon += m;

          errq_V_out_V_vmon += errq*errq;
          errm_V_out_V_vmon += errm*errm;
      }
 }

}


void fit_linear_range(TGraphErrors *g, TCanvas *c1, double low, double high){
   // Fit
   TF1 *line = new TF1("line","[0]*x+[1]", low, high);
   // TF1 *line = new TF1("line","[0]*x");


   g->Fit("line", "qr");
   c1->cd();
   g->Draw("same");

   // Get Parameters
   double m, q, errm, errq; // y = m x + q
   q = line->GetParameter(1);
   errq = line->GetParError(1);
   // q=0; errq=0;
   m = line->GetParameter(0);
   errm = line->GetParError(0);


   if(fit_V_out_bool){
       printf("V_out = (%lf +- %lf) * HV_set + (%lf +- %lf)\n",  m, errm, q, errq);
       // for the mean:
       q_V_out += q;
       m_V_out += m;

       errq_V_out += errq*errq;
       errm_V_out += errm*errm;
   }else{
       if(fit_V_out_V_vmon_bool){
           printf("V_out = (%.2lf +- %.2lf) * V_vmon + (%.2lf +- %.2lf)\n",  m, errm, q, errq);
           // for the mean:
           q_V_out_V_vmon += q;
           m_V_out_V_vmon += m;

           errq_V_out_V_vmon += errq*errq;
           errm_V_out_V_vmon += errm*errm;
       }
  }




}
