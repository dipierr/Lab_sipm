/******************************************************************************\
 * OPERATION POINT
 * 
 * Setup LED
 * 
 * 20180522_HD3-2_1_LED_5583_AGILENT_36_AS_2_50000_04.dat
 * 20180523_HD3-2_1_LED_5583_AGILENT_35.5_AS_2_50000_01.dat
 * 20180511_HD3-2_1_LED_5583_AGILENT_35_AS_2_01-50000eve.dat
 * 20180523_HD3-2_1_LED_5583_AGILENT_34.5_AS_2_50000_01.dat
 * 20180522_HD3-2_1_LED_5583_AGILENT_34_AS_2_50000_01.dat
 * 20180522_HD3-2_1_LED_5583_AGILENT_33_AS_2_50000_02.dat
 * 20180522_HD3-2_1_LED_5583_AGILENT_32_AS_2_50000_03.dat
 * 
\******************************************************************************/

#define n 7

void OperationPoint(){
    // HV
    double HV[] =                {32.00,33.00, 34.00,34.50,35.00,35.50,36.00};
    double errHV[]=              {0.01,0.01,0.01,0.01, 0.01,0.01,0.01};
    
    // Sigma_peak_0
    double Sigma_peak_0[] =      {2.19778e+00,2.18775e+00,2.35603,2.26689e+00,2.70367e+00,2.52442e+00,3.15461e+00};
    double errSigma_peak_0[] =   {5.09429e-02,5.55387e-02,7.16945e-02,6.32656e-02,7.25822e-02,7.36245e-02,1.07703e-01};
    
    // Sigma_peak_1
    double Sigma_peak_1[] =      {5.75964e+00,5.83062e+00,5.60473e+00,5.13521e+00,5.81193e+00,5.95813e+00,8.44510e+00};
    double errSigma_peak_1[] =   {2.49730e-01,1.53347e-01,1.22846e-01,8.78786e-02,1.77725e-01,1.44267e-01,3.08829e-01};
    
    // Mean_peak_1
    double Mean_peak_1[] =      {1.53951e+01,1.75317e+01,1.88485e+01,2.01527e+01,2.04909e+01,2.16402e+01,2.29223e+01};
    double errMean_peak_1[] =   {2.48184e-01,1.22079e-01,1.09392e-01,8.35216e-02,1.58792e-01,1.43441e-01,3.58637e-01};
      
    // Mean_peak_2
    double Mean_peak_2[] =      {2.63210e+01,3.11680e+01,3.60932e+01,3.88570e+01,4.20238e+01,4.34901e+01,4.48226e+01};
    double errMean_peak_2[] =   {1.97001e-01,2.54707e-01,1.64642e-01,1.62696e-01,2.18771e-01,3.00476e-01,6.91813e-01};
    
    // Mean_hist_global
    double Mean_hist_global[] =      {38.89,51.58,65.41,73.94,80.43,89.72,98.15};
    double errMean_hist_global[] =   {0.1069,0.1441,0.1832,0.2069,0.229,0.254,0.2756};
    
    // Standard_dev_hist_global
    double St_dev_hist_global[] =      {23.9,32.22,40.97,46.26,51.21,56.8,69.37};
    double errSt_dev_hist_global[] =   {0.07559,0.1019,0.1295,0.1463,0.1619,0.1796,0.1949};
    
   
    double GAIN[n];    
    double errGAIN[n];
    double GainWeighted[n];
    double errGainWeighted[n];
    double Mean_hist_St_dev[n];
    double errMean_hist_St_dev[n];
    
   
    //------------------------------
    
    for (int i=0; i<n; i++){
        
         // GAIN
        GAIN[i] = Mean_peak_2[i]-Mean_peak_1[i];
        errGAIN[i] = TMath::Sqrt((errMean_peak_2[i]/Mean_peak_2[i])*(errMean_peak_2[i]/Mean_peak_2[i])+(errMean_peak_1[i]/Mean_peak_1[i])*(errMean_peak_1[i]/Mean_peak_1[i]));
        
        // GainWeighted
        GainWeighted[i] = GAIN[i]/TMath::Sqrt(Sigma_peak_1[i]*Sigma_peak_1[i]-Sigma_peak_0[i]*Sigma_peak_0[i]);
        errGainWeighted[i] = errGAIN[i]*errGAIN[i]/(Sigma_peak_1[i]*Sigma_peak_1[i]-Sigma_peak_0[i]*Sigma_peak_0[i]);
        errGainWeighted[i] += errSigma_peak_1[i]*errSigma_peak_1[i]*Sigma_peak_1[i]*Sigma_peak_1[i]*GAIN[i]*GAIN[i]/TMath::Power(Sigma_peak_1[i]*Sigma_peak_1[i]-Sigma_peak_0[i]*Sigma_peak_0[i],3);
        errGainWeighted[i] += errSigma_peak_0[i]*errSigma_peak_0[i]*Sigma_peak_0[i]*Sigma_peak_0[i]*GAIN[i]*GAIN[i]/TMath::Power(Sigma_peak_1[i]*Sigma_peak_1[i]-Sigma_peak_0[i]*Sigma_peak_0[i],3);
        errGainWeighted[i] = TMath::Sqrt(errGainWeighted[i]);
        
        // Mean_hist_global/St_dev
        Mean_hist_St_dev[i]= Mean_hist_global[i]/St_dev_hist_global[i];
        errMean_hist_St_dev[i] = TMath::Power(St_dev_hist_global[i],-1)*TMath::Power(errMean_hist_global[i]*errMean_hist_global[i]+Mean_hist_global[i]*Mean_hist_global[i]*errSt_dev_hist_global[i]*errSt_dev_hist_global[i]/(St_dev_hist_global[i]*St_dev_hist_global[i]),0.5);
        
    }
    
    //------------------------------
    // Gain_Weighted
    //------------------------------
    
    TGraphErrors *gV_GW  = new TGraphErrors(n, HV, GainWeighted, errHV, errGainWeighted);
    
    
    //------------------------------
    
    gV_GW->SetMarkerStyle(20);
    gV_GW->SetMarkerColor(kOrange+2);
    gV_GW->SetTitle();
    gV_GW->GetXaxis()->SetTitle("Voltage (V)");
    gV_GW->GetYaxis()->SetTitle("Weighted Gain");
    
    //------------------------------
    
    TCanvas *cV_GW = new TCanvas("cV_GW", "cV_GW");
    cV_GW->SetGrid();
    gV_GW->Draw("AP");
    
     //------------------------------
     // Mean/Dev_st
     //------------------------------
    
    TGraphErrors *gV_MS  = new TGraphErrors(n, HV, Mean_hist_St_dev, errHV, errMean_hist_St_dev);
    
    
    //------------------------------
    
    gV_MS->SetMarkerStyle(20);
    gV_MS->SetMarkerColor(kOrange+2);
    gV_MS->SetTitle();
    gV_MS->GetXaxis()->SetTitle("Voltage (V)");
    gV_MS->GetYaxis()->SetTitle("Mean/St_Dev ()");
    
    //------------------------------
    
    TCanvas *cV_MS = new TCanvas("cV_MS", "cV_MS");
    cV_MS->SetGrid();
    gV_MS->Draw("AP");


    
    
    
}
