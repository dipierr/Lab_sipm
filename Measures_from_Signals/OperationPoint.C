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
    double Sigma_peak_0[] =      {};
    double errSigma_peak_0[] =   {};
    
    // Sigma_peak_1
    double Sigma_peak_1[] =      {};
    double errSigma_peak_1[] =   {};
    
    // GAIN
    double GAIN[] =              {};
    double errGAIN[] =           {};
    
   //------------------------------
    
    TGraphErrors *gV_out_MAXi  = new TGraphErrors(n, DAC, V_out_MAXi, errDAC, errV_out_MAXi );
    
    
    //------------------------------
    
    gV_out_MAXi->SetMarkerStyle(20);
        
    gV_out_MAXi->SetMarkerColor(kOrange+2);
    
    
    gV_out_MAXi->SetTitle();
    
    //------------------------------
    
    TCanvas *cV_out_MAXi = new TCanvas("cV_out_MAXi", "cV_out_MAXi");
    cV_out_MAXi->SetGrid();
    gV_out_MAXi->Draw("AP");

    
    
    
}
