/******************************************************************************\
 * OPERATION POINT and MORE
 *
 * Setup LED:
 *    > HV:         AGILENT E3641A
 *    > SUPPLY:     Kenwood Regulated DC Power Supply PW18-1.8Q
 *    > AMPLIFIER:  ADVANSID OUT 2
 *    > DIGITIZER:  DRS4 Evaluation Board
 *
 *
 * > File obtained using Ana_LED(...) function in Ana_Traces_SiPM.cxx using
 *   fit_hist_peaks_gaus_sum_012(...)
 *   Fit based on the function reported in:
 *         Mallamaci Manuela, Mariotti Mosé - Report SiPM tests
 *
 *
 * > Fit range and the initial values of the parameters are set manually for
 *   each file
 *
 * File used:
 *    20180626_HD3-2_2_LASER_PLS_75_PAPER_AGILENT_29_AS_2_50000_01.dat
 *    20180626_HD3-2_2_LASER_PLS_75_PAPER_AGILENT_30_AS_2_50000_01.dat
 *    20180626_HD3-2_2_LASER_PLS_75_PAPER_AGILENT_31_AS_2_50000_01.dat
 *    20180626_HD3-2_2_LASER_PLS_75_PAPER_AGILENT_32_AS_2_50000_01.dat
 *    20180626_HD3-2_2_LASER_PLS_75_PAPER_AGILENT_33_AS_2_50000_01.dat
 *    20180626_HD3-2_2_LASER_PLS_75_PAPER_AGILENT_34_AS_2_50000_01.dat
 *    20180626_HD3-2_2_LASER_PLS_75_PAPER_AGILENT_35_AS_2_50000_01.dat
 *    20180626_HD3-2_2_LASER_PLS_75_PAPER_AGILENT_36_AS_2_50000_01.dat
 *
\******************************************************************************/




// #define n_GAIN 8    // for HV = 31 ... 36 V
// #define n_MEAN 10   // for HV = 29 ... 36 V

#define n_GAIN 6    // for HV = 31 ... 36 V
#define n_MEAN 8    // for HV = 29 ... 36 V


#define h 600
#define w 800

void OperationPoint_LASER_PLS_FitGausSum_02(){
    // HV
    double HV[] =  {29.00,30.00,    31.00,        32.00,        33.00,        34.00,               35.00,               36.00};
    double errHV[]={0.01, 0.01,     0.01,         0.01,         0.01,         0.01,                0.01,                0.01};


    // PEAK 0
    // H_peak_0
    double H_peak_0[n_GAIN];
    double errH_peak_0[n_GAIN];

    // Sigma_peak_0
    double Sigma_peak_0[n_GAIN];
    double errSigma_peak_0[n_GAIN];

    // PEAK 1
    // H_peak_1
    double H_peak_1[n_GAIN];
    double errH_peak_1[n_GAIN];

    // Sigma_peak_1
    double Sigma_peak_1[n_GAIN];
    double errSigma_peak_1[n_GAIN];

    // Mean_peak_1
    double Mean_peak_1[n_GAIN];
    double errMean_peak_1[n_GAIN];

    // PEAK 2
    // Mean_peak_2
    double Mean_peak_2[n_GAIN];
    double errMean_peak_2[n_GAIN];

    // GLOBAL HIST
    // Mean hist global
    double Mean_hg[n_MEAN];
    double errMean_hg[n_MEAN];

    // Standard_dev_hist_global
    double Std_hg[n_MEAN];
    double errStd_hg[n_MEAN];

     // Integral of the average
    double Integral[n_GAIN];
    double errIntegral[n_GAIN];

    // Entries
    double Entries[n_GAIN];

    // GAIN
    double GAIN[n_GAIN];
    double errGAIN[n_GAIN];

    // Other variables
    double GainWeighted[n_GAIN];
    double errGainWeighted[n_GAIN];
    double Mean_hist_St_dev[n_MEAN];
    double errMean_hist_St_dev[n_MEAN];
    double Integral_GAIN[n_GAIN];
    double errIntegral_GAIN[n_GAIN];
    double HV_GAIN[n_GAIN];
    double errHV_GAIN[n_GAIN];
    double Area0[n_GAIN];
    double Area1[n_GAIN];
    double errArea0[n_GAIN];
    double errArea1[n_GAIN];
    double Prob_0pe[n_GAIN];
    double errProb_0pe[n_GAIN];
    double Prob_1pe[n_GAIN];
    double errProb_1pe[n_GAIN];
    double Prob_1peS[n_GAIN];
    double errProb_1peS[n_GAIN];
    double Mu[n_GAIN];
    double errMu[n_GAIN];
    double Prob_Cross_Talk[n_GAIN];
    double errProb_Cross_Talk[n_GAIN];
    double errEntries[n_GAIN];


    // ../../Data_DRS/20180626/20180626_HD3-2_2_LASER_PLS_75_PAPER_AGILENT_29_AS_2_50000_01.dat
    // Window for LED peak: (177, 184) ns  (minLED_amp, maxLED_amp)
    Mean_hg[0]           = 10.9;
    errMean_hg[0]        = 0.0408309;
    Std_hg[0]            = 7.16533;
    errStd_hg[0]         = 0.0288718;


    // ../../Data_DRS/20180626/20180626_HD3-2_2_LASER_PLS_75_PAPER_AGILENT_30_AS_2_50000_01.dat
    // Window for LED peak: (177, 184) ns  (minLED_amp, maxLED_amp)
    Mean_hg[1]           = 17.4443;
    errMean_hg[1]        = 0.0648047;
    Std_hg[1]            = 11.4776;
    errStd_hg[1]         = 0.0458238;


    // ../../Data_DRS/20180626/20180626_HD3-2_2_LASER_PLS_75_PAPER_AGILENT_31_AS_2_50000_01.dat
    // Window for LED peak: (177, 184) ns  (minLED_amp, maxLED_amp)
    H_peak_0[0]          = 954.545;
    errH_peak_0[0]       = 20.9612;
    Sigma_peak_0[0]      = 2.22455;
    errSigma_peak_0[0]   = 0.0338465;
    H_peak_1[0]          = 1247.23;
    errH_peak_1[0]       = 22.0836;
    Sigma_peak_1[0]      = 2.5304;
    errSigma_peak_1[0]   = 0.0374124;
    Mean_peak_1[0]       = 11.1771;
    errMean_peak_1[0]    = 0.0861384;
    Mean_peak_2[0]       = 20.3622;
    errMean_peak_2[0]    = 0.099172;
    GAIN[0]              = 9.18513;
    errGAIN[0]           = 0.0491453;
    Mean_hg[2]           = 25.1508;
    errMean_hg[2]        = 0.0931926;
    Std_hg[2]            = 16.7552;
    errStd_hg[2]         = 0.0658971;
    Entries[0]           = 32325;


    // ../../Data_DRS/20180626/20180626_HD3-2_2_LASER_PLS_75_PAPER_AGILENT_32_AS_2_50000_01.dat
    // Window for LED peak: (177, 184) ns  (minLED_amp, maxLED_amp)
    H_peak_0[1]          = 740.108;
    errH_peak_0[1]       = 17.1586;
    Sigma_peak_0[1]      = 2.48347;
    errSigma_peak_0[1]   = 0.0367411;
    H_peak_1[1]          = 910.073;
    errH_peak_1[1]       = 15.7551;
    Sigma_peak_1[1]      = 2.86582;
    errSigma_peak_1[1]   = 0.0375723;
    Mean_peak_1[1]       = 13.2364;
    errMean_peak_1[1]    = 0.086227;
    Mean_peak_2[1]       = 24.9152;
    errMean_peak_2[1]    = 0.0941384;
    GAIN[1]              = 11.6787;
    errGAIN[1]           = 0.037775;
    Mean_hg[3]           = 34.5864;
    errMean_hg[3]        = 0.132661;
    Std_hg[3]            = 23.5427;
    errStd_hg[3]         = 0.0938053;
    Entries[1]           = 31494;


    // ../../Data_DRS/20180626/20180626_HD3-2_2_LASER_PLS_75_PAPER_AGILENT_33_AS_2_50000_01.dat
    // Window for LED peak: (177, 184) ns  (minLED_amp, maxLED_amp)
    H_peak_0[2]          = 679.723;
    errH_peak_0[2]       = 17.3937;
    Sigma_peak_0[2]      = 2.4195;
    errSigma_peak_0[2]   = 0.0417604;
    H_peak_1[2]          = 684.609;
    errH_peak_1[2]       = 13.0383;
    Sigma_peak_1[2]      = 3.0766;
    errSigma_peak_1[2]   = 0.0421667;
    Mean_peak_1[2]       = 14.89;
    errMean_peak_1[2]    = 0.0964008;
    Mean_peak_2[2]       = 28.8386;
    errMean_peak_2[2]    = 0.107269;
    GAIN[2]              = 13.9486;
    errGAIN[2]           = 0.0470471;
    Mean_hg[4]           = 45.5606;
    errMean_hg[4]        = 0.176687;
    Std_hg[4]            = 31.3325;
    errStd_hg[4]         = 0.124937;
    Entries[2]           = 31447;



    // ../../Data_DRS/20180626/20180626_HD3-2_2_LASER_PLS_75_PAPER_AGILENT_34_AS_2_50000_01.dat
    // Window for LED peak: (177, 184) ns  (minLED_amp, maxLED_amp)
    H_peak_0[3]          = 558.301;
    errH_peak_0[3]       = 15.6788;
    Sigma_peak_0[3]      = 2.73936;
    errSigma_peak_0[3]   = 0.0526505;
    H_peak_1[3]          = 486.676;
    errH_peak_1[3]       = 10.2421;
    Sigma_peak_1[3]      = 3.52059;
    errSigma_peak_1[3]   = 0.0514007;
    Mean_peak_1[3]       = 16.8473;
    errMean_peak_1[3]    = 0.115863;
    Mean_peak_2[3]       = 33.3641;
    errMean_peak_2[3]    = 0.126916;
    GAIN[3]              = 16.5168;
    errGAIN[3]           = 0.0518018;
    Mean_hg[5]           = 57.8561;
    errMean_hg[5]        = 0.228594;
    Std_hg[5]            = 40.0765;
    errStd_hg[5]         = 0.161641;
    Entries[3]           = 30736;

    // ../../Data_DRS/20180626/20180626_HD3-2_2_LASER_PLS_75_PAPER_AGILENT_35_AS_2_50000_01.dat
    // Window for LED peak: (177, 184) ns  (minLED_amp, maxLED_amp)
    H_peak_0[4]          = 314.557;
    errH_peak_0[4]       = 7.79229;
    Sigma_peak_0[4]      = 4.85258;
    errSigma_peak_0[4]   = 0.106354;
    H_peak_1[4]          = 388.811;
    errH_peak_1[4]       = 8.3045;
    Sigma_peak_1[4]      = 5.22298;
    errSigma_peak_1[4]   = 0.103341;
    Mean_peak_1[4]       = 15.3968;
    errMean_peak_1[4]    = 0.244818;
    Mean_peak_2[4]       = 34.9716;
    errMean_peak_2[4]    = 0.268865;
    GAIN[4]              = 19.5748;
    errGAIN[4]           = 0.111143;
    Mean_hg[6]           = 70.5002;
    errMean_hg[6]        = 0.24596;
    Std_hg[6]            = 49.1957;
    errStd_hg[6]         = 0.17392;
    Entries[4]           = 40006;


    // ../../Data_DRS/20180626/20180626_HD3-2_2_LASER_PLS_75_PAPER_AGILENT_36_AS_2_50000_01.dat
    // Window for LED peak: (177, 184) ns  (minLED_amp, maxLED_amp)
    H_peak_0[5]          = 419.314;
    errH_peak_0[5]       = 12.7225;
    Sigma_peak_0[5]      = 3.0784;
    errSigma_peak_0[5]   = 0.0646934;
    H_peak_1[5]          = 315.359;
    errH_peak_1[5]       = 7.49133;
    Sigma_peak_1[5]      = 4.26062;
    errSigma_peak_1[5]   = 0.0643761;
    Mean_peak_1[5]       = 20.1109;
    errMean_peak_1[5]    = 0.158778;
    Mean_peak_2[5]       = 41.4653;
    errMean_peak_2[5]    = 0.175993;
    GAIN[5]              = 21.3544;
    errGAIN[5]           = 0.0759147;
    Mean_hg[7]           = 87.3688;
    errMean_hg[7]        = 0.344128;
    Std_hg[7]            = 60.7568;
    errStd_hg[7]         = 0.243335;
    Entries[5]           = 31171;




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

    //------------------------------
    // Gain_Weighted
    //------------------------------

    TGraphErrors *gV_GW  = new TGraphErrors(n_GAIN, HV_GAIN, GainWeighted, errHV_GAIN, errGainWeighted);


    //------------------------------

    gV_GW->SetMarkerStyle(20);
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

    gV_MS->SetMarkerStyle(20);
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

    gV_PCT->SetMarkerStyle(20);
    gV_PCT->SetMarkerColor(kOrange+2);
    gV_PCT->SetTitle();
    gV_PCT->GetXaxis()->SetTitle("Voltage (V)");
    gV_PCT->GetYaxis()->SetTitle("Probability Cross Talk");

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

    gV_IG->SetMarkerStyle(20);
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
