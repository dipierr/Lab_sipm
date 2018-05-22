/******************************************************************************\
 * OPERATION POINT
 * 
 * Setup LED
 * 
 * 20180511_HD3-2_1_LED_5583_AGILENT_35_AS_2_01-50000eve.dat
 * 20180522_HD3-2_1_LED_5583_AGILENT_34_AS_2_50000_01.dat
 * 20180522_HD3-2_1_LED_5583_AGILENT_33_AS_2_50000_02.dat
 * 
\******************************************************************************/

#define n 3

void OperationPoint(){
    // HV
    double HV[] =                {33.00, 34.00, 35.00};
    double errHV[]=              {0.01, 0.01, 0.01};
    
    // Sigma_peak_0
    double Sigma_peak_0[] =      {2.18775e+00,2.35603,2.42307e+00};
    double errSigma_peak_0[] =   {5.55387e-02,7.16945e-02,5.57088e-02};
    
    // Sigma_peak_1
    double Sigma_peak_1[] =      {5.83062e+00,5.60473e+00,5.89185e+00};
    double errSigma_peak_1[] =   {1.53347e-01,1.22846e-01,1.13433e-01};
    
    // Mean_peak_1
    double Mean_peak_1[] =      {1.75317e+01,1.88485e+01,2.02275e+01};
    double errMean_peak_1[] =   {1.22079e-01,1.09392e-01,1.02673e-01};
      
    // Mean_peak_2
    double Mean_peak_2[] =      {3.11680e+01,3.60932e+01,4.02396e+01};
    double errMean_peak_2[] =   {2.54707e-01,1.64642e-01,1.65936e-01};
    
   
    double GAIN[n];    
    double errGAIN[n];
    double GainWeighted[n];
    double errGainWeighted[n];
    
   
    //------------------------------
    
    for (int i=0; i<n; i++){
        
         // GAIN
        GAIN[i]=Mean_peak_2[i]-Mean_peak_1[i];
        errGAIN[i]=TMath::Sqrt((errMean_peak_2[i]/Mean_peak_2[i])*(errMean_peak_2[i]/Mean_peak_2[i])+(errMean_peak_1[i]/Mean_peak_1[i])*(errMean_peak_1[i]/Mean_peak_1[i]));
        
        // GainWeighted
        GainWeighted[i]=GAIN[i]/TMath::Sqrt(Sigma_peak_1[i]*Sigma_peak_1[i]-Sigma_peak_0[i]*Sigma_peak_0[i]);
        errGainWeighted[i]= errGAIN[i]*errGAIN[i]/(Sigma_peak_1[i]*Sigma_peak_1[i]-Sigma_peak_0[i]*Sigma_peak_0[i]);
        errGainWeighted[i]+= errSigma_peak_1[i]*errSigma_peak_1[i]*Sigma_peak_1[i]*Sigma_peak_1[i]*GAIN[i]*GAIN[i]/TMath::Power(Sigma_peak_1[i]*Sigma_peak_1[i]-Sigma_peak_0[i]*Sigma_peak_0[i],3);
        errGainWeighted[i]+= errSigma_peak_0[i]*errSigma_peak_0[i]*Sigma_peak_0[i]*Sigma_peak_0[i]*GAIN[i]*GAIN[i]/TMath::Power(Sigma_peak_1[i]*Sigma_peak_1[i]-Sigma_peak_0[i]*Sigma_peak_0[i],3);
        errGainWeighted[i]=TMath::Sqrt(errGainWeighted[i]);
        
        
    }
    
    //------------------------------
    
    TGraphErrors *gV_GW  = new TGraphErrors(n, HV, GainWeighted, errHV, errGainWeighted);
    
    
    //------------------------------
    
    gV_GW->SetMarkerStyle(20);
        
    gV_GW->SetMarkerColor(kOrange+2);
    
    
    gV_GW->SetTitle();
    
    //------------------------------
    
    TCanvas *cV_GW = new TCanvas("cV_GW", "cV_GW");
    cV_GW->SetGrid();
    gV_GW->Draw("AP");

    
    
    
}
