// ****************************************************************************
// *********************************   OLD   **********************************
// ****************************************************************************

/******************************************************************************\
 * OPERATION POINT and MORE
 *
 * Setup LED:
 *    > HV:         AGILENT E3641A
 *    > SUPPLY:     Kenwood Regulated DC Power Supply PW18-1.8Q
 *    > AMPLIFIER:  ADVANSID OUT 2
 *    > DIGITIZER:  DRS4 Evaluation Board
 *    > LASER:      Advanced Laser Diode System, Picosecond Laser System,
 *                  Controller EIG2000DX and a 406 nm laser head
 *
 * > File obtained using Ana_LED(...) function in Ana_Traces_SiPM.cxx using
 *   fit_hist_peaks_gaus_sum_012(...):
 *      [0]*TMath::Exp( - (x-[6])*(x-[6])/( 2*[4]*[4] ) ) +
 *      + [1]*TMath::Exp(-(x-[6]-[3])*(x-[6]-[3])/( 2*([4]*[4] + [5]*[5] )) ) +
 *      + [2]*TMath::Exp(-(x-[6]-2*[3])*(x-[6]-2*[3])/(2*([4]*[4] + 4*[5]*[5])))
 *
 * > Ana_Traces_SiPM.cxx github: 09/08/2018 # 1
 *
 * > Fit range and the initial values of the parameters are set manually for
 *   each file
 *
 * File used:
 *    20180614_HD3-2_1_LASER_PLS_81_PAPER_AGILENT_29_AS_2_50000_01.dat
 *    20180614_HD3-2_1_LASER_PLS_81_PAPER_AGILENT_30_AS_2_50000_01.dat
 *    20180614_HD3-2_1_LASER_PLS_81_PAPER_AGILENT_31_AS_2_50000_01.dat
 *    20180614_HD3-2_1_LASER_PLS_81_PAPER_AGILENT_32_AS_2_50000_01.dat
 *    20180614_HD3-2_1_LASER_PLS_81_PAPER_AGILENT_33_AS_2_50000_01.dat
 *    20180614_HD3-2_1_LASER_PLS_81_PAPER_AGILENT_34_AS_2_50000_01.dat
 *    20180614_HD3-2_1_LASER_PLS_81_PAPER_AGILENT_35_AS_2_50000_01.dat
 *    20180614_HD3-2_1_LASER_PLS_81_PAPER_AGILENT_36_AS_2_50000_01.dat
 *
\******************************************************************************/

// #define n_GAIN 8    // for HV = 31 ... 36 V
// #define n_MEAN 10   // for HV = 29 ... 36 V

#define n_GAIN 6    // for HV = 31 ... 36 V
#define n_MEAN 8    // for HV = 29 ... 36 V


#define h 600
#define w 1000

void OperationPoint_LASER_PLS_FitGausSum_01(){

    // ERRORS
    bool fix_error_bool = true;

    // HV
    double HV[] =  {29.00,30.00,    31.00,        32.00,        33.00,        34.00,               35.00,               36.00};
    double errHV[]={0.01, 0.01,     0.01,         0.01,         0.01,         0.01,                0.01,                0.01};
    // double HV[] =  {29.00,30.00,    31.00,        32.00,        33.00,        34.00,        34.50,        35.00,        35.50,        36.00};
    // double errHV[]={0.01, 0.01,     0.01,         0.01,         0.01,         0.01,         0.01,         0.01,         0.01,         0.01};



    // PEAK 0
    // H_peak_0
    double H_peak_0[n_GAIN] = {0.};
    double errH_peak_0[n_GAIN] = {0.};

    // Sigma_peak_0
    double Sigma_peak_0[n_GAIN] = {0.};
    double errSigma_peak_0[n_GAIN] = {0.};

    // PEAK 1
    // H_peak_1
    double H_peak_1[n_GAIN] = {0.};
    double errH_peak_1[n_GAIN] = {0.};

    // Sigma_peak_1
    double Sigma_peak_1[n_GAIN] = {0.};
    double errSigma_peak_1[n_GAIN] = {0.};

    // Mean_peak_1
    double Mean_peak_1[n_GAIN] = {0.};
    double errMean_peak_1[n_GAIN] = {0.};

    // PEAK 2
    // Mean_peak_2
    double Mean_peak_2[n_GAIN] = {0.};
    double errMean_peak_2[n_GAIN] = {0.};

    // GLOBAL HIST
    // Mean hist global
    double Mean_hg[n_MEAN] = {0.};
    double errMean_hg[n_MEAN] = {0.};

    // Standard_dev_hist_global
    double Std_hg[n_MEAN] = {0.};
    double errStd_hg[n_MEAN] = {0.};

     // Integral of the average
    double Integral[n_GAIN] = {0.};
    double errIntegral[n_GAIN] = {0.};

    // Entries
    double Entries[n_GAIN] = {0.};

    // GAIN
    double GAIN[n_GAIN] = {0.};
    double errGAIN[n_GAIN] = {0.};

    // Other variables
    double GainWeighted[n_GAIN] = {0.};
    double errGainWeighted[n_GAIN] = {0.};
    double Mean_hist_St_dev[n_MEAN] = {0.};
    double errMean_hist_St_dev[n_MEAN] = {0.};
    double Integral_GAIN[n_GAIN] = {0.};
    double errIntegral_GAIN[n_GAIN] = {0.};
    double HV_GAIN[n_GAIN] = {0.};
    double errHV_GAIN[n_GAIN] = {0.};
    double Area0[n_GAIN] = {0.};
    double Area1[n_GAIN] = {0.};
    double errArea0[n_GAIN] = {0.};
    double errArea1[n_GAIN] = {0.};
    double Prob_0pe[n_GAIN] = {0.};
    double errProb_0pe[n_GAIN] = {0.};
    double Prob_1pe[n_GAIN] = {0.};
    double errProb_1pe[n_GAIN] = {0.};
    double Prob_1peS[n_GAIN] = {0.};
    double errProb_1peS[n_GAIN] = {0.};
    double Mu[n_GAIN] = {0.};
    double errMu[n_GAIN] = {0.};
    double Prob_Cross_Talk[n_GAIN] = {0.};
    double errProb_Cross_Talk[n_GAIN] = {0.};
    double errEntries[n_GAIN] = {0.};



    // HV = 29 V


    // HV = 30 V



    // HV = 31 V [0]
    // Window for LED peak: (168, 177) ns  (minLED_amp, maxLED_amp)
    // Other: dleddt = 9; smooth_trace_bool = 0;
    //        histLED_low = -50; histLED_high = 700; histLED_binw = 0.7;
    //        fit_low = -10; fit_high = 29;
    H_peak_0[0]          = 798.911;
    errH_peak_0[0]       = 12.4535;
    Sigma_peak_0[0]      = 2.45416;
    errSigma_peak_0[0]   = 0.0310133;
    H_peak_1[0]          = 1180.98;
    errH_peak_1[0]       = 16.6779;
    Sigma_peak_1[0]      = 2.96697;
    errSigma_peak_1[0]   = 0.0387114;
    Mean_peak_1[0]       = 5.9277;
    errMean_peak_1[0]    = 0.0543554;
    Mean_peak_2[0]       = 8.38187;
    errMean_peak_2[0]    = 0.0625806;
    GAIN[0]              = 2.45416;
    errGAIN[0]           = 0.0310133;
    Mean_hg[2]           = 23.9805;
    errMean_hg[2]        = 0.0713393;
    Std_hg[2]            = 15.9518;
    errStd_hg[2]         = 0.0504445;
    Entries[0]           = 49999;


    // HV = 32 V [1]
    // Window for LED peak: (168, 177) ns  (minLED_amp, maxLED_amp)
    // Other: dleddt = 9; smooth_trace_bool = 0;
    //        histLED_low = -50; histLED_high = 700; histLED_binw = 0.7;
    //        fit_low = -10; fit_high = 34;
    H_peak_0[1]          = 602.069;
    errH_peak_0[1]       = 9.94509;
    Sigma_peak_0[1]      = 2.80037;
    errSigma_peak_0[1]   = 0.0322954;
    H_peak_1[1]          = 955.733;
    errH_peak_1[1]       = 11.8414;
    Sigma_peak_1[1]      = 3.274;
    errSigma_peak_1[1]   = 0.0365093;
    Mean_peak_1[1]       = 5.91894;
    errMean_peak_1[1]    = 0.0551977;
    Mean_peak_2[1]       = 8.71931;
    errMean_peak_2[1]    = 0.0639514;
    GAIN[1]              = 2.80037;
    errGAIN[1]           = 0.0322954;
    Mean_hg[3]           = 32.2026;
    errMean_hg[3]        = 0.0987055;
    Std_hg[3]            = 22.071;
    errStd_hg[3]         = 0.0697953;
    Entries[1]           = 49999;



    // HV = 33 V [2]
    // Window for LED peak: (168, 177) ns  (minLED_amp, maxLED_amp)
    // Other: dleddt = 9; smooth_trace_bool = 0;
    //        histLED_low = -50; histLED_high = 700; histLED_binw = 0.7;
    //        fit_low = -10; fit_high = 38;
    H_peak_0[2]          = 465.179;
    errH_peak_0[2]       = 8.49261;
    Sigma_peak_0[2]      = 3.20809;
    errSigma_peak_0[2]   = 0.0409653;
    H_peak_1[2]          = 782.935;
    errH_peak_1[2]       = 9.37982;
    Sigma_peak_1[2]      = 3.66242;
    errSigma_peak_1[2]   = 0.0447216;
    Mean_peak_1[2]       = 5.98651;
    errMean_peak_1[2]    = 0.065181;
    Mean_peak_2[2]       = 9.19459;
    errMean_peak_2[2]    = 0.0769852;
    GAIN[2]              = 3.20809;
    errGAIN[2]           = 0.0409653;
    Mean_hg[4]           = 41.6162;
    errMean_hg[4]        = 0.131494;
    Std_hg[4]            = 29.4028;
    errStd_hg[4]         = 0.0929806;
    Entries[2]           = 49999;





    // HV = 34 V [3]
    // Window for LED peak: (168, 177) ns  (minLED_amp, maxLED_amp)
    // Other: dleddt = 9; smooth_trace_bool = 0;
    //        histLED_low = -50; histLED_high = 700; histLED_binw = 0.7;
    //        fit_low = -10; fit_high = 46;
    H_peak_0[3]          = 367.064;
    errH_peak_0[3]       = 7.19149;
    Sigma_peak_0[3]      = 3.55882;
    errSigma_peak_0[3]   = 0.0455046;
    H_peak_1[3]          = 619.135;
    errH_peak_1[3]       = 7.96463;
    Sigma_peak_1[3]      = 4.17185;
    errSigma_peak_1[3]   = 0.0486095;
    Mean_peak_1[3]       = 5.51912;
    errMean_peak_1[3]    = 0.0708883;
    Mean_peak_2[3]       = 9.07794;
    errMean_peak_2[3]    = 0.0842367;
    GAIN[3]              = 3.55882;
    errGAIN[3]           = 0.0455046;
    Mean_hg[5]           = 52.3904;
    errMean_hg[5]        = 0.167837;
    Std_hg[5]            = 37.529;
    errStd_hg[5]         = 0.118678;
    Entries[3]           = 49999;




    // HV = 35.00 V [4]
    // Window for LED peak: (168, 177) ns  (minLED_amp, maxLED_amp)
    // Other: dleddt = 9; smooth_trace_bool = 0;
    //        histLED_low = -50; histLED_high = 700; histLED_binw = 0.7;
    //        fit_low = -10; fit_high = 51;
    H_peak_0[4]          = 310.268;
    errH_peak_0[4]       = 6.52916;
    Sigma_peak_0[4]      = 3.78661;
    errSigma_peak_0[4]   = 0.054867;
    H_peak_1[4]          = 517.979;
    errH_peak_1[4]       = 6.82964;
    Sigma_peak_1[4]      = 4.54551;
    errSigma_peak_1[4]   = 0.0561603;
    Mean_peak_1[4]       = 4.9444;
    errMean_peak_1[4]    = 0.0799529;
    Mean_peak_2[4]       = 8.73101;
    errMean_peak_2[4]    = 0.0969682;
    GAIN[4]              = 3.78661;
    errGAIN[4]           = 0.054867;
    Mean_hg[6]           = 63.8248;
    errMean_hg[6]        = 0.20694;
    Std_hg[6]            = 46.2728;
    errStd_hg[6]         = 0.146329;
    Entries[4]           = 49999;







    // HV = 36 V ([5])
    // Window for LED peak: (168, 177) ns  (minLED_amp, maxLED_amp)
    // Other: dleddt = 9; smooth_trace_bool = 0;
    //        histLED_low = -50; histLED_high = 700; histLED_binw = 0.7;
    //        fit_low = -10; fit_high = 60;
    H_peak_0[5]          = 253.065;
    errH_peak_0[5]       = 5.72925;
    Sigma_peak_0[5]      = 4.17716;
    errSigma_peak_0[5]   = 0.0685689;
    H_peak_1[5]          = 394.896;
    errH_peak_1[5]       = 5.7841;
    Sigma_peak_1[5]      = 5.1988;
    errSigma_peak_1[5]   = 0.0663333;
    Mean_peak_1[5]       = 4.39273;
    errMean_peak_1[5]    = 0.0982355;
    Mean_peak_2[5]       = 8.56989;
    errMean_peak_2[5]    = 0.119799;
    GAIN[5]              = 4.17716;
    errGAIN[5]           = 0.0685689;
    Mean_hg[7]           = 76.9744;
    errMean_hg[7]        = 0.252346;
    Std_hg[7]            = 56.4257;
    errStd_hg[7]         = 0.178435;
    Entries[5]           = 49999;







    //------------------------------

    // select only points for n_GAIN
    for (int i = 0; i < n_GAIN; i++) {
        HV_GAIN[i] = HV[i+n_MEAN-n_GAIN];
        errHV_GAIN[i] = errHV[i+n_MEAN-n_GAIN];
    }

    for (int i=0; i<n_GAIN; i++){

        // GainWeighted
        GainWeighted[i] = GAIN[i]/TMath::Sqrt(Sigma_peak_1[i]*Sigma_peak_1[i]-Sigma_peak_0[i]*Sigma_peak_0[i]);
        errGainWeighted[i] = errGAIN[i]*errGAIN[i]/(Sigma_peak_1[i]*Sigma_peak_1[i]-Sigma_peak_0[i]*Sigma_peak_0[i]);
        errGainWeighted[i] += errSigma_peak_1[i]*errSigma_peak_1[i]*Sigma_peak_1[i]*Sigma_peak_1[i]*GAIN[i]*GAIN[i]/TMath::Power(Sigma_peak_1[i]*Sigma_peak_1[i]-Sigma_peak_0[i]*Sigma_peak_0[i],3);
        errGainWeighted[i] += errSigma_peak_0[i]*errSigma_peak_0[i]*Sigma_peak_0[i]*Sigma_peak_0[i]*GAIN[i]*GAIN[i]/TMath::Power(Sigma_peak_1[i]*Sigma_peak_1[i]-Sigma_peak_0[i]*Sigma_peak_0[i],3);
        errGainWeighted[i] = TMath::Sqrt(errGainWeighted[i]);

        // Integral_GAIN
        // Integral_GAIN[i]=TMath::Abs(Integral[i])/(GAIN[i]*Cross_Talk[i]);
        // errIntegral_GAIN[i]=Integral[i]/(GAIN[i]*Cross_Talk[i])*TMath::Power(errGAIN[i]*errGAIN[i]/(GAIN[i]*GAIN[i])+errCross_Talk[i]*errCross_Talk[i]/(Cross_Talk[i]*Cross_Talk[i]),0.5);
        //Integral_GAIN[i]=TMath::Abs(Integral[i])/(GAIN[i]);
        //errIntegral_GAIN[i]=Integral[i]*errGAIN[i]/(GAIN[i]*GAIN[i]);

        //Probabilità CrossTalk LED
        errEntries[i]=TMath::Power(Entries[i],0.5);
        Area0[i]=H_peak_0[i]*Sigma_peak_0[i]*TMath::Power(2*TMath::Pi(),0.5)/1;
        Prob_0pe[i]=Area0[i]/Entries[i];
        Mu[i]=-TMath::Log(Prob_0pe[i]);
        Prob_1pe[i]=Mu[i]*TMath::Exp(-Mu[i]);
        Area1[i]=H_peak_1[i]*Sigma_peak_1[i]*TMath::Power(2*TMath::Pi(),0.5)/1;
        Prob_1peS[i]=Area1[i]/Entries[i];
        Prob_Cross_Talk[i]=1-(Prob_1peS[i]/Prob_1pe[i]);
        cout << "Cross Talk " << Prob_Cross_Talk[i] << endl;
        cout << "Prob_1peS[i]/Prob_1pe[i]\t" << Prob_1peS[i]/Prob_1pe[i] << endl;
        cout << "p 0 \t" << Prob_0pe[i] << endl;
        cout << "p 1 s\t" << Prob_1peS[i] << endl;
        cout << "p 1\t" << Prob_1pe[i] << endl<<endl    ;


        errArea1[i]=Area1[i]*TMath::Power(TMath::Power(errH_peak_1[i]/H_peak_1[i],2)+TMath::Power(errSigma_peak_1[i]/Sigma_peak_1[i],2),0.5);
        errProb_1peS[i]=errArea1[i]/Entries[i];
        errMu[i]=errArea0[i]/(Prob_0pe[i]*Entries[i]);
        errProb_1pe[i]=TMath::Exp(-Mu[i])*TMath::Abs(1-Mu[i])*errMu[i];
        errArea0[i]=Area0[i]*TMath::Power(TMath::Power(errH_peak_0[i]/H_peak_0[i],2)+TMath::Power(errSigma_peak_0[i]/Sigma_peak_0[i],2),0.5);
        errProb_0pe[i]=errArea0[i]/Entries[i];
        errProb_Cross_Talk[i]=TMath::Power(TMath::Power(errProb_1peS[i]/Prob_1pe[i],2)+TMath::Power(Prob_1peS[i]*errProb_1pe[i]/Prob_1pe[i],2),0.5);
    }
    double temp;
    for (int i = 0; i < n_MEAN; i++) {
      // Mean_hg/St_dev
      Mean_hist_St_dev[i]= Mean_hg[i]/Std_hg[i];


      errMean_hist_St_dev[i] = TMath::Power(Std_hg[i],-1)*TMath::Power(errMean_hg[i]*errMean_hg[i]+Mean_hg[i]*Mean_hg[i]*errStd_hg[i]*errStd_hg[i]/(Std_hg[i]*Std_hg[i]),0.5);

      // temp = errMean_hg[i]*errMean_hg[i]+Mean_hg[i]*Mean_hg[i]*errStd_hg[i]*errStd_hg[i]/(Std_hg[i]*Std_hg[i]);
      // temp = TMath::Power(temp, 0.5);
      // temp *= 1/Std_hg[i];
      //


    }

    ///////////////////////////////////////////////////////////////////////////
    //          FIXED ERROR
    ///////////////////////////////////////////////////////////////////////////

    if(fix_error_bool){
        // double err_fix_CT = 0.015;
        double err_fix_CT = 0.02;
        for(int i=0; i<n_GAIN; i++){
            errProb_Cross_Talk[i] = err_fix_CT;
        }
    }

    ///////////////////////////////////////////////////////////////////////////

    for(int i=0; i<n_GAIN; i++){
            cout<<"Error CT = "<<errProb_Cross_Talk[i]<<endl;
        }


    //------------------------------
    // Gain_Weighted
    //------------------------------

    TGraphErrors *gV_GW  = new TGraphErrors(n_GAIN, HV_GAIN, GainWeighted, errHV_GAIN, errGainWeighted);


    //------------------------------

    gV_GW->SetMarkerStyle(21);
    gV_GW->SetMarkerColor(kOrange+2);
    gV_GW->SetTitle();
    gV_GW->GetXaxis()->SetTitle("Voltage (V)");
    gV_GW->GetYaxis()->SetTitle("Weighted Gain");

    //------------------------------

    TCanvas *cV_GW = new TCanvas("cV_GW", "cV_GW",w,h);
    cV_GW->SetGrid();
    gV_GW->Draw("AP");

     //------------------------------
     // Mean/Dev_st
     //------------------------------

    TGraphErrors *gV_MS  = new TGraphErrors(n_MEAN, HV, Mean_hist_St_dev, errHV, errMean_hist_St_dev);


    //------------------------------

    gV_MS->SetMarkerStyle(21);
    gV_MS->SetMarkerColor(kOrange+2);
    gV_MS->SetTitle();
    gV_MS->GetXaxis()->SetTitle("Voltage (V)");
    gV_MS->GetYaxis()->SetTitle("Mean/St_Dev ()");

    //------------------------------

    TCanvas *cV_MS = new TCanvas("cV_MS", "cV_MS",w,h);
    cV_MS->SetGrid();
    gV_MS->Draw("AP");


    //------------------------------
    // Cross Talk
    //------------------------------
    //for(int i=0; i<n_GAIN; i++){
      //  cout<<HV_GAIN[i]<<"\t"<<Prob_Cross_Talk[i]<<endl;
    //}
    TGraphErrors *gV_PCT  = new TGraphErrors(n_GAIN, HV_GAIN, Prob_Cross_Talk, errHV_GAIN, errProb_Cross_Talk);


    //------------------------------

    gV_PCT->SetMarkerStyle(21);
    // gV_PCT->SetMarkerSize(2);
    gV_PCT->SetMarkerColor(kOrange+2);
    gV_PCT->SetTitle();
    gV_PCT->GetXaxis()->SetTitle("Voltage (V)");
    gV_PCT->GetYaxis()->SetTitle("Cross Talk");

    //------------------------------

    TCanvas *cV_PCT = new TCanvas("cV_PCT", "cV_PCT",w,h);
    cV_PCT->SetGrid();
    gV_PCT->Draw("AP");




/* The following method is not yet implemented

    //------------------------------
    // Integral_GAIN
    //------------------------------

    TGraphErrors *gV_IG  = new TGraphErrors(n_GAIN, HV_GAIN, Integral_GAIN, errHV_GAIN, errIntegral_GAIN);


    //------------------------------

    gV_IG->SetMarkerStyle(47);
    gV_IG->SetMarkerColor(kOrange+2);
    gV_IG->SetTitle();
    gV_IG->GetXaxis()->SetTitle("Voltage (V)");
    gV_IG->GetYaxis()->SetTitle("Integral/Gain");

    //------------------------------

    TCanvas *cV_IG = new TCanvas("cV_IG", "cV_IG",w,h);
    cV_IG->SetGrid();
    // gV_IG->Draw("AP");
    for (size_t i = 0; i < n_GAIN; i++) {
      cout<<HV_GAIN[i]<<endl<<endl;
      cout<<errHV_GAIN[i]<<endl<<endl;
      cout<<Integral_GAIN[i]<<endl<<endl;
      cout<<errIntegral_GAIN[i]<<endl<<endl;
    }

  */
}


/******************************************************************************\
 *
 *  REFERENCES
 *
 *  > Zorzi Filippo - Caratterizzazione di fotosensori al silicio per telescopi
 *    Cherenkov di nuova generazione
 *
 *  > Mallamaci Manuela, Mariotti Mosé - Report SiPM tests
 *
 *
 *
 ******************************************************************************/



///////////////////////////////////////////////////////////////////////////////
//      PREVIOUS DATA (calulated on 15/06/2018)
///////////////////////////////////////////////////////////////////////////////

// // HV = 29 V
// // Window for LED peak: (168, 176) ns  (minLED_amp, maxLED_amp)
// Mean_hg[0]           = 11.0454;
// errMean_hg[0]        = 0.0349697;
// Std_hg[0]            = 7.19416;
// errStd_hg[0]         = 0.0247273;
//
//
// // HV = 30 V
// // Window for LED peak: (168, 176) ns  (minLED_amp, maxLED_amp)
// Mean_hg[1]           = 17.2405;
// errMean_hg[1]        = 0.0579644;
// Std_hg[1]            = 11.6465;
// errStd_hg[1]         = 0.040987;
//
//
//
// // HV = 31 V
// // Window for LED peak: (168, 176) ns  (minLED_amp, maxLED_amp)
// H_peak_0[0]          = 1676.19;
// errH_peak_0[0]       = 24.8442;
// Sigma_peak_0[0]      = 2.49041;
// errSigma_peak_0[0]   = 0.0270042;
// H_peak_1[0]          = 1716.92;
// errH_peak_1[0]       = 27.2019;
// Sigma_peak_1[0]      = 2.78504;
// errSigma_peak_1[0]   = 0.0374966;
// Mean_peak_1[0]       = 13.2836;
// errMean_peak_1[0]    = 0.0868258;
// Mean_peak_2[0]       = 23.0244;
// errMean_peak_2[0]    = 0.101794;
// GAIN[0]              = 9.74076;
// errGAIN[0]           = 0.0531344;
// Mean_hg[2]           = 24.5025;
// errMean_hg[2]        = 0.0817773;
// Std_hg[2]            = 16.7689;
// errStd_hg[2]         = 0.0578253;
// Entries[0]           = 42048;
//
//
// // HV = 32 V
// // Window for LED peak: (168, 176) ns  (minLED_amp, maxLED_amp)
// H_peak_0[1]          = 1274.36;
// errH_peak_0[1]       = 20.3547;
// Sigma_peak_0[1]      = 2.71012;
// errSigma_peak_0[1]   = 0.0289319;
// H_peak_1[1]          = 1389.22;
// errH_peak_1[1]       = 20.9821;
// Sigma_peak_1[1]      = 3.09644;
// errSigma_peak_1[1]   = 0.0365057;
// Mean_peak_1[1]       = 15.4873;
// errMean_peak_1[1]    = 0.0900793;
// Mean_peak_2[1]       = 27.8335;
// errMean_peak_2[1]    = 0.104941;
// GAIN[1]              = 12.3462;
// errGAIN[1]           = 0.0538359;
// Mean_hg[3]           = 33.1871;
// errMean_hg[3]        = 0.111983;
// Std_hg[3]            = 23.0281;
// errStd_hg[3]         = 0.0791843;
// Entries[1]           = 42287;
//
// // HV = 33 V
// // Window for LED peak: (168, 176) ns  (minLED_amp, maxLED_amp)
// H_peak_0[2]          = 1005.74;
// errH_peak_0[2]       = 17.4816;
// Sigma_peak_0[2]      = 3.01203;
// errSigma_peak_0[2]   = 0.0341288;
// H_peak_1[2]          = 1133.25;
// errH_peak_1[2]       = 16.489;
// Sigma_peak_1[2]      = 3.4767;
// errSigma_peak_1[2]   = 0.0401889;
// Mean_peak_1[2]       = 17.6329;
// errMean_peak_1[2]    = 0.0936316;
// Mean_peak_2[2]       = 32.5561;
// errMean_peak_2[2]    = 0.106978;
// GAIN[2]              = 14.9232;
// errGAIN[2]           = 0.051744;
// Mean_hg[4]           = 43.0989;
// errMean_hg[4]        = 0.14868;
// Std_hg[4]            = 30.536;
// errStd_hg[4]         = 0.105133;
// Entries[2]           = 42181;
//
//
// // HV = 34 V
// // Window for LED peak: (168, 176) ns  (minLED_amp, maxLED_amp)
// H_peak_0[3]          = 812.22;
// errH_peak_0[3]       = 14.6323;
// Sigma_peak_0[3]      = 3.28217;
// errSigma_peak_0[3]   = 0.0397998;
// H_peak_1[3]          = 890.785;
// errH_peak_1[3]       = 14.1818;
// Sigma_peak_1[3]      = 3.91139;
// errSigma_peak_1[3]   = 0.0520462;
// Mean_peak_1[3]       = 19.5695;
// errMean_peak_1[3]    = 0.147611;
// Mean_peak_2[3]       = 37.3578;
// errMean_peak_2[3]    = 0.174144;
// GAIN[3]              = 17.7883;
// errGAIN[3]           = 0.0923972;
// Mean_hg[5]           = 54.2596;
// errMean_hg[5]        = 0.190426;
// Std_hg[5]            = 38.901;
// errStd_hg[5]         = 0.134651;
// Entries[3]           = 41732;
//
//
// // HV = 34.5 V
// // Window for LED peak: (168, 176) ns  (minLED_amp, maxLED_amp)
// // H_peak_0[4]          = 770.553;
// // errH_peak_0[4]       = 14.8191;
// // Sigma_peak_0[4]      = 3.27105;
// // errSigma_peak_0[4]   = 0.0408399;
// // H_peak_1[4]          = 756.362;
// // errH_peak_1[4]       = 12.4416;
// // Sigma_peak_1[4]      = 3.9776;
// // errSigma_peak_1[4]   = 0.0500078;
// // Mean_peak_1[4]       = 20.2041;
// // errMean_peak_1[4]    = 0.125302;
// // Mean_peak_2[4]       = 38.952;
// // errMean_peak_2[4]    = 0.144446;
// // GAIN[4]              = 18.7479;
// // errGAIN[4]           = 0.0718609;
// // Mean_hg[6]           = 59.8421;
// // errMean_hg[6]        = 0.218923;
// // Std_hg[6]            = 43.4229;
// // errStd_hg[6]         = 0.154802;
// // Entries[4]           = 39342;
//
//
//
//
// // HV = 35 V
// // Window for LED peak: (168, 176) ns  (minLED_amp, maxLED_amp)
// H_peak_0[4]          = 699.866;
// errH_peak_0[4]       = 13.313;
// Sigma_peak_0[4]      = 3.47154;
// errSigma_peak_0[4]   = 0.0456399;
// H_peak_1[4]          = 707.385;
// errH_peak_1[4]       = 11.5108;
// Sigma_peak_1[4]      = 4.2467;
// errSigma_peak_1[4]   = 0.0515783;
// Mean_peak_1[4]       = 21.2611;
// errMean_peak_1[4]    = 0.126853;
// Mean_peak_2[4]       = 41.391;
// errMean_peak_2[4]    = 0.144826;
// GAIN[4]              = 20.1299;
// errGAIN[4]           = 0.0698776;
// Mean_hg[6]           = 66.1856;
// errMean_hg[6]        = 0.239378;
// Std_hg[6]            = 48.1251;
// errStd_hg[6]         = 0.169266;
// Entries[4]           = 40418;
//
//
//
//
// // HV = 35.5 V
// // Window for LED peak: (168, 176) ns  (minLED_amp, maxLED_amp)
// // H_peak_0[6]          = 622.18;
// // errH_peak_0[6]       = 12.9319;
// // Sigma_peak_0[6]      = 3.57571;
// // errSigma_peak_0[6]   = 0.0506767;
// // H_peak_1[6]          = 586.952;
// // errH_peak_1[6]       = 11.3562;
// // Sigma_peak_1[6]      = 4.39935;
// // errSigma_peak_1[6]   = 0.0741474;
// // Mean_peak_1[6]       = 21.9082;
// // errMean_peak_1[6]    = 0.27288;
// // Mean_peak_2[6]       = 43.0821;
// // errMean_peak_2[6]    = 0.328738;
// // GAIN[6]              = 21.1738;
// // errGAIN[6]           = 0.183317;
// // Mean_hg[8]           = 73.7344;
// // errMean_hg[8]        = 0.272271;
// // Std_hg[8]            = 53.9193;
// // errStd_hg[8]         = 0.192525;
// // Entries[6]           = 39218;
//
//
//
//
// // HV = 36 V
// // Window for LED peak: (168, 176) ns  (minLED_amp, maxLED_amp)
// H_peak_0[5]          = 574.096;
// errH_peak_0[5]       = 12.0926;
// Sigma_peak_0[5]      = 3.78057;
// errSigma_peak_0[5]   = 0.0569589;
// H_peak_1[5]          = 535.369;
// errH_peak_1[5]       = 10.3089;
// Sigma_peak_1[5]      = 4.79561;
// errSigma_peak_1[5]   = 0.0746405;
// Mean_peak_1[5]       = 22.7587;
// errMean_peak_1[5]    = 0.228875;
// Mean_peak_2[5]       = 45.5447;
// errMean_peak_2[5]    = 0.270868;
// GAIN[5]              = 22.786;
// errGAIN[5]           = 0.144865;
// Mean_hg[7]           = 80.0673;
// errMean_hg[7]        = 0.294391;
// Std_hg[7]            = 58.7329;
// errStd_hg[7]         = 0.208166;
// Entries[5]           = 39803;
