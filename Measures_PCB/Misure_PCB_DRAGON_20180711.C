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
\******************************************************************************/

#include "TMath.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TF1.h"
#include "TMultiGraph.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TLegend.h"

#define n_meas 6


void fit_linear(TGraphErrors *g, TCanvas *c1, TCanvas *c2);
void Misure_PCB_DRAGON_20180711();

double m_V_out, errm_V_out, q_V_out, errq_V_out;

void Misure_PCB_DRAGON_20180711(){
    // HV_set (NOT TRIMMER as before); the voltage of the pin 38 is about 30 - 35 mV more than HV_set[i]
    // from Terminal: $ ApplyHV 192.168.1.5 HV_set[i]
    double HV_set[] =                {1000,   1100,     1200,   1300,   1400,   1500};
    // LOAD: 39 kOhm
    double V_out_MAXi[] =            {26.23,  28.85,   31.48,  34.11,  36.74,  39.36};
    double V_vmon_MAXi[] =           {-0.927, -1.023,  -1.119, -1.215, -1.311, -1.407};
    double V_imon_MAXi[] =           {-0.36,  -0.39,   -0.42,  -0.45,  -0.49,  -0.53};

    // LOAD: 224 (220 kOhm)
    double V_out_MINi[] =            {26.23,  28.86,   31.49,  34.11,  36.74,  39.37};
    double V_vmon_MINi[] =           {-0.927, -1.024,  -1.120, -1.215, -1.312, -1.408};
    double V_imon_MINi[] =           {-0.11,  -0.12,   -0.13,  -0.13,  -0.14,  -0.15};

    // LOAD: 223 (22 kOhm)
    double V_out_SuperMAXi[] =       {26.22,  28.85,   31.47,  34.10,  36.72,  39.35};
    double V_vmon_SuperMAXi[] =      {-0.925, -1.021,  -1.116, -1.212, -1.309, -1.405};
    double V_imon_SuperMAXi[] =      {-0.60,  -0.66,   -0.72,  -0.78,  -0.83,  -0.89};

    // errors: always the same:
    double errHV_set[n_meas], errV_out_MAXi[n_meas], errV_vmon_MAXi[n_meas], errV_imon_MAXi[n_meas], errV_out_MINi[n_meas], errV_vmon_MINi[n_meas], errV_imon_MINi[n_meas], errV_out_SuperMAXi[n_meas], errV_vmon_SuperMAXi[n_meas], errV_imon_SuperMAXi[n_meas];
    for(int i=0; i<n_meas; i++){
      errHV_set[i] = 1;
      errV_out_MAXi[i] = errV_out_MINi[i] = errV_out_SuperMAXi[i] = 0.01;
      errV_vmon_MAXi[i] = errV_vmon_MINi[i] = errV_vmon_SuperMAXi[i] = 0.001;
      errV_imon_MAXi[i] = errV_imon_MINi[i] = errV_imon_SuperMAXi[i] = 0.01;
    }

    //------------------------------

    TGraphErrors *gV_out_MAXi  = new TGraphErrors(n_meas, HV_set, V_out_MAXi, errHV_set, errV_out_MAXi );
    TGraphErrors *gV_vmon_MAXi = new TGraphErrors(n_meas, HV_set, V_vmon_MAXi,errHV_set, errV_vmon_MAXi);
    TGraphErrors *gV_imon_MAXi = new TGraphErrors(n_meas, HV_set, V_imon_MAXi,errHV_set, errV_imon_MAXi);

    TGraphErrors *gV_out_MINi  = new TGraphErrors(n_meas, HV_set, V_out_MINi, errHV_set, errV_out_MINi );
    TGraphErrors *gV_vmon_MINi = new TGraphErrors(n_meas, HV_set, V_vmon_MINi,errHV_set, errV_vmon_MINi);
    TGraphErrors *gV_imon_MINi = new TGraphErrors(n_meas, HV_set, V_imon_MINi,errHV_set, errV_imon_MINi);

    TGraphErrors *gV_out_SuperMAXi  = new TGraphErrors(n_meas, HV_set, V_out_SuperMAXi, errHV_set, errV_out_SuperMAXi );
    TGraphErrors *gV_vmon_SuperMAXi = new TGraphErrors(n_meas, HV_set, V_vmon_SuperMAXi,errHV_set, errV_vmon_SuperMAXi);
    TGraphErrors *gV_imon_SuperMAXi = new TGraphErrors(n_meas, HV_set, V_imon_SuperMAXi,errHV_set, errV_imon_SuperMAXi);

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

    // gV_out_MAXi->SetMarkerSize(3);
    // gV_vmon_MAXi->SetMarkerSize(3);
    // gV_imon_MAXi->SetMarkerSize(3);
    //
    // gV_out_MINi->SetMarkerSize(3);
    // gV_vmon_MINi->SetMarkerSize(3);
    // gV_imon_MINi->SetMarkerSize(3);

    gV_out_MAXi->SetMarkerColor(kOrange+2);
    gV_vmon_MAXi->SetMarkerColor(kRed+2);
    gV_imon_MAXi->SetMarkerColor(kMagenta+2);

    gV_out_MINi->SetMarkerColor(kOrange);
    gV_vmon_MINi->SetMarkerColor(kRed);
    gV_imon_MINi->SetMarkerColor(kMagenta);

    gV_out_SuperMAXi->SetMarkerColor(kOrange+3);
    gV_vmon_SuperMAXi->SetMarkerColor(kRed+3);
    gV_imon_SuperMAXi->SetMarkerColor(kMagenta+3);

    gV_out_MAXi->SetTitle();
    gV_vmon_MAXi->SetTitle();
    gV_imon_MAXi->SetTitle();


    gV_out_MINi->SetTitle();
    gV_vmon_MINi->SetTitle();
    gV_imon_MINi->SetTitle();

    gV_out_SuperMAXi->SetTitle();
    gV_vmon_SuperMAXi->SetTitle();
    gV_imon_SuperMAXi->SetTitle();


    //------------------------------
    //----------- GRAPHS -----------
    //------------------------------

    TCanvas *cV_out_MAXi = new TCanvas("cV_out_MAXi", "cV_out_MAXi");
    cV_out_MAXi->SetGrid();
    gV_out_MAXi->Draw("AP");

    TCanvas *cV_vmon_MAXi = new TCanvas("cV_vmon_MAXi", "cV_vmon_MAXi");
    cV_vmon_MAXi->SetGrid();
    gV_vmon_MAXi->Draw("AP");

    TCanvas *cV_imon_MAXi = new TCanvas("cV_imon_MAXi", "cV_imon_MAXi");
    cV_imon_MAXi->SetGrid();
    gV_imon_MAXi->Draw("AP");


    TCanvas *cV_out_MINi = new TCanvas("cV_out_MINi", "cV_out_MINi");
    cV_out_MINi->SetGrid();
    gV_out_MINi->Draw("AP");

    TCanvas *cV_vmon_MINi = new TCanvas("cV_vmon_MINi", "cV_vmon_MINi");
    cV_vmon_MINi->SetGrid();
    gV_vmon_MINi->Draw("AP");

    TCanvas *cV_imon_MINi = new TCanvas("cV_imon_MINi", "cV_imon_MINi");
    cV_imon_MINi->SetGrid();
    gV_imon_MINi->Draw("AP");


    TCanvas *cV_out_SuperMAXi = new TCanvas("cV_out_SuperMAXi", "cV_out_SuperMAXi");
    cV_out_SuperMAXi->SetGrid();
    gV_out_SuperMAXi->Draw("AP");

    TCanvas *cV_vmon_SuperMAXi = new TCanvas("cV_vmon_SuperMAXi", "cV_vmon_SuperMAXi");
    cV_vmon_SuperMAXi->SetGrid();
    gV_vmon_SuperMAXi->Draw("AP");

    TCanvas *cV_imon_SuperMAXi = new TCanvas("cV_imon_SuperMAXi", "cV_imon_SuperMAXi");
    cV_imon_SuperMAXi->SetGrid();
    gV_imon_SuperMAXi->Draw("AP");


    //------------------------------
    //-------- MULTIGRAPHs ---------
    //------------------------------

    TMultiGraph *V_out = new TMultiGraph("V_out", ";HV_set (V); V_out (V)");
    V_out->Add(gV_out_MINi);
    V_out->Add(gV_out_MAXi);
    V_out->Add(gV_out_SuperMAXi);


    TMultiGraph *V_vmon = new TMultiGraph("V_vmon", ";HV_set (V); V_vmon (V)");
    V_vmon->Add(gV_vmon_MINi);
    V_vmon->Add(gV_vmon_MAXi);
    V_vmon->Add(gV_vmon_SuperMAXi);


    TMultiGraph *V_imon = new TMultiGraph("V_imon", ";HV_set (V); V_imon (V)");
    V_imon->Add(gV_imon_MINi);
    V_imon->Add(gV_imon_MAXi);
    V_imon->Add(gV_imon_SuperMAXi);


    //------------------------------

    TCanvas *cV_out = new TCanvas("cV_out", "cV_out");
    cV_out->SetGrid();
    V_out->Draw("AP");
    auto legend_V_out = new TLegend(0.15,0.75,0.35,0.9);
    legend_V_out->AddEntry(gV_out_MINi,"V_out_MINi", "p");
    legend_V_out->AddEntry(gV_out_MAXi,"V_out_MAXi","p");
    legend_V_out->AddEntry(gV_out_SuperMAXi,"V_out_SuperMAXi", "p");
    legend_V_out->Draw();

    TCanvas *cV_vmon = new TCanvas("cV_vmon", "cV_vmon");
    cV_vmon->SetGrid();
    V_vmon ->Draw("AP");
    auto legend_V_vmon = new TLegend(0.15,0.75,0.35,0.9);
    legend_V_vmon->AddEntry(gV_vmon_MINi,"V_vmon_MINi", "p");
    legend_V_vmon->AddEntry(gV_vmon_MAXi,"V_vmon_MAXi","p");
    legend_V_vmon->AddEntry(gV_vmon_SuperMAXi,"V_vmon_SuperMAXi","p");
    legend_V_vmon->Draw();

    TCanvas *cV_imon = new TCanvas("cV_imon", "cV_imon");
    cV_imon->SetGrid();
    V_imon ->Draw("AP");
    auto legend_V_imon = new TLegend(0.15,0.75,0.35,0.9);
    legend_V_imon->AddEntry(gV_imon_MINi,"V_imon_MINi", "p");
    legend_V_imon->AddEntry(gV_imon_MAXi,"V_imon_MAXi","p");
    legend_V_imon->AddEntry(gV_imon_SuperMAXi,"V_imon_SuperMAXi","p");
    legend_V_imon->Draw();


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
/*
    //----------- V_out ------------
    m_V_out = errm_V_out = q_V_out = errq_V_out = 0.;
    int n_mean;

    cout<<endl;
    cout<<"///// V_out /////"<<endl;
    n_mean = 3;
    fit_linear(gV_out_MINi, cV_out, cV_out_MINi);
    fit_linear(gV_out_MAXi, cV_out, cV_out_MAXi);
    fit_linear(gV_out_SuperMAXi, cV_out, cV_out_SuperMAXi);

    // mean values
    q_V_out /= n_mean;
    m_V_out /= n_mean;
    errq_V_out = TMath::Sqrt(errq_V_out) / n_mean;
    errm_V_out = TMath::Sqrt(errm_V_out) / n_mean;
    printf("V_out = (%lf +- %lf) * HV_set + (%lf +- %lf)\n",  m_V_out, errm_V_out, q_V_out, errq_V_out);

    // err on V_out
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
    TCanvas *cerrV_out_est = new TCanvas("cerrV_out_est", "cerrV_out_est");
    TGraph* gerrV_out_est = new TGraph(n_est,V_out_est,errV_out_est);
    gerrV_out_est->SetTitle();
    gerrV_out_est->GetXaxis()->SetTitle("V_out [V]");
    gerrV_out_est->GetYaxis()->SetTitle("errV_out [mV]");
    gerrV_out_est->Draw("AL");

*/

    //------------------------------
    //---------- CONSTANT ----------
    //------------------------------
    double min_K = -29;
    double max_K = -27;
    double binw_K = 0.1;
    int nbins_K;

    nbins_K = (int)((max_K-min_K)/binw_K);

    TH1* ptrhistK = new TH1I("histK", "", nbins_K, min_K, max_K);

    for (int i=0; i<n_meas; i++){
      ptrhistK->Fill(V_out_MINi[i]/V_vmon_MINi[i]);
      ptrhistK->Fill(V_out_MAXi[i]/V_vmon_MAXi[i]);
      ptrhistK->Fill(V_out_SuperMAXi[i]/V_vmon_SuperMAXi[i]);
    }

    TCanvas *chistK = new TCanvas("chistK", "chistK");
    ptrhistK->Draw();


    return;


}

void fit_linear(TGraphErrors *g, TCanvas *c1, TCanvas *c2){
  // Fit
  TF1 *line = new TF1("line","[0]*x+[1]");
  // TF1 *line = new TF1("line","[0]*x");


  g->Fit("line", "q");
  c1->cd();
  g->Draw("same");
  c2->cd();
  g->Draw("same");

  // Get Parameters
  double m, q, errm, errq; // y = m x + q
  q = line->GetParameter(1);
  errq = line->GetParError(1);
  // q=0; errq=0;
  m = line->GetParameter(0);
  errm = line->GetParError(0);

  printf("V_out = (%lf +- %lf) * HV_set + (%lf +- %lf)\n",  m, errm, q, errq);

  // for the mean
  q_V_out += q;
  m_V_out += m;

  errq_V_out += errq*errq;
  errm_V_out += errm*errm;


}
