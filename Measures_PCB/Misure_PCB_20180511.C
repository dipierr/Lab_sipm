/******************************************************************************\
 *  MEASURES FROM PCB OF INTERFACE BOARD, MADE BY SILVANO GALLIAN
 *  FILTER 3 STAGES (92dB @ 100kHz)
 *  
 *  MODEL based on
 *              > 20180510_SiPM_BOARD_09.pdf
 *  
 *  POWER SUPPLY:
 *              > Kenwood Regulated DC Power Supply PW18-1.8Q
 *  
 *  MEASURES TAKEN WITH 
 *              > Fluke 179 True RMS Multimeter 
 * 
 *  LOADs USED:
 *              > MAXi -> 39  kOhm
 *              > MINi -> 220 kOhm (code: 224)
 *  
 *  SCHEMATIC:
 *              > 20180426_PCB_TEST_KENWOOD_01.jpg
 * 
\******************************************************************************/

#define n 5

void Misure_PCB_20180511(){
    // TRIMMER
    double DAC[] =                {0.954, 1.049, 1.145, 1.240, 1.336};
    double errDAC[] =             {0.001,0.001,0.001,0.001,0.001};
    
    // LOAD: 39 kOhm
    double V_out_MAXi[] =         {25.00, 27.49, 30.00,  32.49  , 34.99};
    double errV_out_MAXi[] =      {0.01,0.01,0.01,0.01,0.01};
    double V_vmon_MAXi[] =        {-0.917, -1.008, -1.100, -1.191 ,     -1.283};
    double errV_vmon_MAXi[] =     {0.001,0.001,0.001,0.001,0.001};
    double V_imon_MAXi[] =        {-0.38, -0.41, -0.44, -0.48,  -0.50};
    double errV_imon_MAXi[] =     {0.01,0.01,0.01,0.01,0.01};
    
    // LOAD: 224 (220 kOhm)
    double V_out_MINi[] =         {25.00, 27.50,  30.00,  32.50 ,   35.00};
    double errV_out_MINi[] =      {0.01,0.01,0.01,0.01,0.01};
    double V_vmon_MINi[] =        {-0.917, -1.008, -1.100, -1.191,  -1.283 };
    double errV_vmon_MINi[] =     {0.001,0.001,0.001,0.001,0.001 };
    double V_imon_MINi[] =        {-0.14,  -0.15, -0.16,  -0.17,  -0.18 };
    double errV_imon_MINi[] =     {0.01,0.01,0.01,0.01,0.01};

    //------------------------------
    
    TGraphErrors *gV_out_MAXi  = new TGraphErrors(n, DAC, V_out_MAXi, errDAC, errV_out_MAXi );
    TGraphErrors *gV_vmon_MAXi = new TGraphErrors(n, DAC, V_vmon_MAXi,errDAC, errV_vmon_MAXi);
    TGraphErrors *gV_imon_MAXi = new TGraphErrors(n, DAC, V_imon_MAXi,errDAC, errV_imon_MAXi);
    
    TGraphErrors *gV_out_MINi  = new TGraphErrors(n, DAC, V_out_MINi, errDAC, errV_out_MINi );
    TGraphErrors *gV_vmon_MINi = new TGraphErrors(n, DAC, V_vmon_MINi,errDAC, errV_vmon_MINi);
    TGraphErrors *gV_imon_MINi = new TGraphErrors(n, DAC, V_imon_MINi,errDAC, errV_imon_MINi);
    
    //------------------------------
    
    gV_out_MAXi->SetMarkerStyle(20);
    gV_vmon_MAXi->SetMarkerStyle(20);
    gV_imon_MAXi->SetMarkerStyle(20);
    
    gV_out_MINi->SetMarkerStyle(20);
    gV_vmon_MINi->SetMarkerStyle(20);
    gV_imon_MINi->SetMarkerStyle(20);
    
    gV_out_MAXi->SetMarkerColor(kOrange+2);
    gV_vmon_MAXi->SetMarkerColor(kRed+2);
    gV_imon_MAXi->SetMarkerColor(kMagenta+2);
    
    gV_out_MINi->SetMarkerColor(kOrange);
    gV_vmon_MINi->SetMarkerColor(kRed);
    gV_imon_MINi->SetMarkerColor(kMagenta);
    
    gV_out_MAXi->SetTitle();
    gV_vmon_MAXi->SetTitle();
    gV_imon_MAXi->SetTitle();
    
    gV_out_MINi->SetTitle();
    gV_vmon_MINi->SetTitle();
    gV_imon_MINi->SetTitle();
    
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
    
    
    //------------------------------
    //-------- MULTIGRAPHs ---------
    //------------------------------
    
    TMultiGraph *V_out = new TMultiGraph("V_out", ";DAC (V); V_out (V)");
    V_out->Add(gV_out_MAXi);
    V_out->Add(gV_out_MINi);
    
    TMultiGraph *V_vmon = new TMultiGraph("V_vmon", ";DAC (V); V_vmon (V)");
    V_vmon->Add(gV_vmon_MAXi);
    V_vmon->Add(gV_vmon_MINi);
    
    TMultiGraph *V_imon = new TMultiGraph("V_imon", ";DAC (V); V_imon (V)");
    V_imon->Add(gV_imon_MAXi);
    V_imon->Add(gV_imon_MINi);

    //------------------------------

    TCanvas *cV_out = new TCanvas("cV_out", "cV_out");
    cV_out->SetGrid();
    V_out->Draw("AP");
    auto legend_V_out = new TLegend(0.15,0.75,0.35,0.9);
    legend_V_out->AddEntry(gV_out_MAXi,"V_out_MAXi","p");
    legend_V_out->AddEntry(gV_out_MINi,"V_out_MINi", "p");
    legend_V_out->Draw();
    
    TCanvas *cV_vmon = new TCanvas("cV_vmon", "cV_vmon");
    cV_vmon->SetGrid();
    V_vmon ->Draw("AP");
    auto legend_V_vmon = new TLegend(0.15,0.75,0.35,0.9);
    legend_V_vmon->AddEntry(gV_vmon_MAXi,"V_vmon_MAXi","p");
    legend_V_vmon->AddEntry(gV_vmon_MINi,"V_vmon_MINi", "p");
    legend_V_vmon->Draw();
    
    TCanvas *cV_imon = new TCanvas("cV_imon", "cV_imon");
    cV_imon->SetGrid();
    V_imon ->Draw("AP");
    auto legend_V_imon = new TLegend(0.15,0.75,0.35,0.9);
    legend_V_imon->AddEntry(gV_imon_MAXi,"V_imon_MAXi","p");
    legend_V_imon->AddEntry(gV_imon_MINi,"V_imon_MINi", "p");
    legend_V_imon->Draw();
    
    
    //------------------------------
    //------ TEST_Z and cout -------
    //------------------------------
    double TEST_Z;
    bool OK;
    
    // V_out
    cout<<endl;
    cout<<"DAC\t\t | V_out, I_min\t   | V_out, I_max     | TEST_Z"<<endl;
    cout<<"\t\t |\t\t   |\t\t      |"<<endl;
    for(int i=0; i<n; i++){
        TEST_Z = TMath::Abs(V_out_MINi[i] - V_out_MAXi[i]) / TMath::Sqrt( errV_out_MINi[i]*errV_out_MINi[i] + errV_out_MAXi[i]*errV_out_MAXi[i]);
        OK = TEST_Z<1.96;
        printf("%.3lf +- %.3lf\t | %.2lf +- %.2lf   | %.2lf +- %.2lf    | %.3lf\t", DAC[i], errDAC[i], V_out_MINi[i], errV_out_MINi[i], V_out_MAXi[i], errV_out_MAXi[i], TEST_Z);
        if(OK==true) cout<<"true"<<endl;
        else         cout<<"false"<<endl;
    }
    
    // V_vmon
    cout<<endl;
    cout<<"DAC\t\t | V_vmon, I_min   | V_vmon, I_max    | TEST_Z"<<endl;
    cout<<"\t\t |\t\t   |\t\t      |"<<endl;
    for(int i=0; i<n; i++){
        TEST_Z = TMath::Abs(V_vmon_MINi[i] - V_vmon_MAXi[i]) / TMath::Sqrt( errV_vmon_MINi[i]*errV_vmon_MINi[i] + errV_vmon_MAXi[i]*errV_vmon_MAXi[i]);
        OK = TEST_Z<1.96;
        printf("%.3lf +- %.3lf\t | %.3lf +- %.3lf | %.3lf +- %.3lf  | %.3lf\t", DAC[i], errDAC[i], V_vmon_MINi[i], errV_vmon_MINi[i], V_vmon_MAXi[i], errV_vmon_MAXi[i], TEST_Z);
        if(OK==true) cout<<"true"<<endl;
        else         cout<<"false"<<endl;
    }
    
    // V_imon
    cout<<endl;
    cout<<"DAC\t\t | V_imon, I_min      | V_imon, I_max"<<endl;
    cout<<"\t\t |\t\t   |\t\t"<<endl;
    for(int i=0; i<n; i++){
        printf("%.3lf +- %.3lf\t | %.2lf +- %.2lf   | %.2lf +- %.2lf\n", DAC[i], errDAC[i], V_imon_MINi[i], errV_imon_MINi[i], V_imon_MAXi[i], errV_imon_MAXi[i]);
    }
    
}
