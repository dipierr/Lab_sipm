// ****************************************************************************
// *********************************   OLD   **********************************
// ****************************************************************************

/******************************************************************************\
 * OPERATION POINT
 *
 * Setup LED:
 *    > HV:         AGILENT E3641A
 *    > SUPPLY:     Kenwood Regulated DC Power Supply PW18-1.8Q
 *    > AMPLIFIER:  ADVANSID OUT 2
 *    > DIGITIZER:  DRS4 Evaluation Board
 *    > FUNC GEN:
 *
 * File obtained using Ana_LED(...) function in Ana_Traces_SiPM.cxx
 *
 *
\******************************************************************************/

#define n_GAIN 8    // for HV = 31 ... 36 V
#define n_MEAN 10   // for HV = 29 ... 36 V


#define h 600
#define w 800

void OperationPoint_LASER_PLS_01(){
    // HV
    double HV[] =  {29.00,30.00,    31.00,        32.00,        33.00,        34.00,        34.50,        35.00,        35.50,        36.00};
    double errHV[]={0.01, 0.01,     0.01,         0.01,         0.01,         0.01,         0.01,         0.01,         0.01,         0.01};


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



    // HV = 29 V
    // Window for LED peak: (168, 176) ns  (minLED_amp, maxLED_amp)
    Mean_hg[0]		= 11.0454;
    errMean_hg[0]		= 0.0349697;
    Std_hg[0]		= 7.19416;
    errStd_hg[0]		= 0.0247273;

    // HV = 30 V
    // Window for LED peak: (168, 176) ns  (minLED_amp, maxLED_amp)
    Mean_hg[1]		= 17.2405;
    errMean_hg[1]		= 0.0579644;
    Std_hg[1]		= 11.6465;
    errStd_hg[1]		= 0.040987;


    // HV = 31 V
    // Window for LED peak: (168, 176) ns  (minLED_amp, maxLED_amp)
    H_peak_0[0]		= 1195.76;
    errH_peak_0[0]		= 17.2176;
    Sigma_peak_0[0]		= 2.63489;
    errSigma_peak_0[0]	= 0.0349133;
    H_peak_1[0]		= 1212.74;
    errH_peak_1[0]		= 17.2658;
    Sigma_peak_1[0]		= 3.59014;
    errSigma_peak_1[0]	= 0.0682874;
    Mean_peak_1[0]		= 13.7287;
    errMean_peak_1[0]	= 0.0599376;
    Mean_peak_2[0]		= 23.3496;
    errMean_peak_2[0]	= 0.0948799;
    Mean_hg[2]		= 24.5025;
    errMean_hg[2]		= 0.0817773;
    Std_hg[2]		= 16.7689;
    errStd_hg[2]		= 0.0578253;
    Entries[0]		= 42048;

    // HV = 32 V
    // Window for LED peak: (168, 176) ns  (minLED_amp, maxLED_amp)
    H_peak_0[1]		= 901.284;
    errH_peak_0[1]		= 14.9443;
    Sigma_peak_0[1]		= 2.77953;
    errSigma_peak_0[1]	= 0.0350864;
    H_peak_1[1]		= 980.541;
    errH_peak_1[1]		= 14.9913;
    Sigma_peak_1[1]		= 3.52126;
    errSigma_peak_1[1]	= 0.0559572;
    Mean_peak_1[1]		= 15.8204;
    errMean_peak_1[1]	= 0.0551003;
    Mean_peak_2[1]		= 28.1085;
    errMean_peak_2[1]	= 0.0799129;
    Mean_hg[3]		= 33.1871;
    errMean_hg[3]		= 0.111983;
    Std_hg[3]		= 23.0281;
    errStd_hg[3]		= 0.0791843;
    Entries[1]		= 42287;


    // HV = 33 V
    // Window for LED peak: (168, 176) ns  (minLED_amp, maxLED_amp)
    H_peak_0[2]		= 751.372;
    errH_peak_0[2]		= 14.3412;
    Sigma_peak_0[2]		= 2.64121;
    errSigma_peak_0[2]	= 0.0357337;
    H_peak_1[2]		= 826.58;
    errH_peak_1[2]		= 13.444;
    Sigma_peak_1[2]		= 3.38943;
    errSigma_peak_1[2]	= 0.0458594;
    Mean_peak_1[2]		= 17.9284;
    errMean_peak_1[2]	= 0.0483072;
    Mean_peak_2[2]		= 33.2002;
    errMean_peak_2[2]	= 0.0768053;
    Mean_hg[4]		= 44.4914;
    errMean_hg[4]		= 0.151745;
    Std_hg[4]		= 30.9359;
    errStd_hg[4]		= 0.1073;
    Entries[2]		= 41562;




    // HV = 34 V
    // Window for LED peak: (168, 176) ns  (minLED_amp, maxLED_amp)
    H_peak_0[3]		= 575.042;
    errH_peak_0[3]		= 11.3432;
    Sigma_peak_0[3]		= 3.30848;
    errSigma_peak_0[3]	= 0.0476995;
    H_peak_1[3]		= 613.364;
    errH_peak_1[3]		= 10.4044;
    Sigma_peak_1[3]		= 4.32238;
    errSigma_peak_1[3]	= 0.0614707;
    Mean_peak_1[3]		= 19.8414;
    errMean_peak_1[3]	= 0.063687;
    Mean_peak_2[3]		= 37.7005;
    errMean_peak_2[3]	= 0.0987244;
    Mean_hg[5]		= 54.2596;
    errMean_hg[5]		= 0.190426;
    Std_hg[5]		= 38.901;
    errStd_hg[5]		= 0.134651;
    Entries[3]		= 41732;



    // HV = 34.5 V
    // Window for LED peak: (168, 176) ns  (minLED_amp, maxLED_amp)
    H_peak_0[4]		= 555.651;
    errH_peak_0[4]		= 11.2784;
    Sigma_peak_0[4]		= 3.21193;
    errSigma_peak_0[4]	= 0.0481069;
    H_peak_1[4]		= 528.17;
    errH_peak_1[4]		= 10.1674;
    Sigma_peak_1[4]		= 4.21188;
    errSigma_peak_1[4]	= 0.0738239;
    Mean_peak_1[4]		= 20.2697;
    errMean_peak_1[4]	= 0.0694128;
    Mean_peak_2[4]		= 39.6444;
    errMean_peak_2[4]	= 0.110517;
    Mean_hg[6]		= 59.8421;
    errMean_hg[6]		= 0.218923;
    Std_hg[6]		= 43.4229;
    errStd_hg[6]		= 0.154802;
    Entries[4]		= 39342;



    // HV = 35 V
    // Window for LED peak: (168, 176) ns  (minLED_amp, maxLED_amp)
    H_peak_0[5]		= 496.62;
    errH_peak_0[5]		= 10.6236;
    Sigma_peak_0[5]		= 3.46826;
    errSigma_peak_0[5]	= 0.0551058;
    H_peak_1[5]		= 494.428;
    errH_peak_1[5]		= 9.40511;
    Sigma_peak_1[5]		= 4.5143;
    errSigma_peak_1[5]	= 0.0819648;
    Mean_peak_1[5]		= 21.3839;
    errMean_peak_1[5]	= 0.0792882;
    Mean_peak_2[5]		= 42.0863;
    errMean_peak_2[5]	= 0.128901;
    Mean_hg[7]		= 66.1856;
    errMean_hg[7]		= 0.239378;
    Std_hg[7]		= 48.1251;
    errStd_hg[7]		= 0.169266;
    Entries[5]		= 40418;



    // HV = 35.5 V
    // Window for LED peak: (168, 176) ns  (minLED_amp, maxLED_amp)
    H_peak_0[6]		= 440.346;
    errH_peak_0[6]		= 9.8741;
    Sigma_peak_0[6]		= 3.60426;
    errSigma_peak_0[6]	= 0.0622195;
    H_peak_1[6]		= 404.539;
    errH_peak_1[6]		= 8.39867;
    Sigma_peak_1[6]		= 4.82882;
    errSigma_peak_1[6]	= 0.0860928;
    Mean_peak_1[6]		= 22.1967;
    errMean_peak_1[6]	= 0.0784199;
    Mean_peak_2[6]		= 44.2131;
    errMean_peak_2[6]	= 0.122782;
    Mean_hg[8]		= 73.7344;
    errMean_hg[8]		= 0.272271;
    Std_hg[8]		= 53.9193;
    errStd_hg[8]		= 0.192525;
    Entries[6]		= 39218;



    // HV = 36 V
    // Window for LED peak: (168, 176) ns  (minLED_amp, maxLED_amp)
    H_peak_0[7]		= 396.66;
    errH_peak_0[7]		= 9.2731;
    Sigma_peak_0[7]		= 3.9136;
    errSigma_peak_0[7]	= 0.0734443;
    H_peak_1[7]		= 376.185;
    errH_peak_1[7]		= 8.03764;
    Sigma_peak_1[7]		= 5.17263;
    errSigma_peak_1[7]	= 0.112446;
    Mean_peak_1[7]		= 22.9992;
    errMean_peak_1[7]	= 0.0995025;
    Mean_peak_2[7]		= 46.0661;
    errMean_peak_2[7]	= 0.148843;
    Mean_hg[9]		= 80.0673;
    errMean_hg[9]		= 0.294391;
    Std_hg[9]		= 58.7329;
    errStd_hg[9]		= 0.208166;
    Entries[7]		= 39803;



    double GAIN[n_GAIN];
    double errGAIN[n_GAIN];
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


    //------------------------------

    // select only points for n_GAIN
    for (int i = 0; i < n_GAIN; i++) {
        HV_GAIN[i] = HV[i+n_MEAN-n_GAIN];
        errHV_GAIN[i] = errHV[i+n_MEAN-n_GAIN];
    }

    for (int i=0; i<n_GAIN; i++){

         // GAIN
        GAIN[i] = Mean_peak_2[i]-Mean_peak_1[i];
        errGAIN[i] = TMath::Sqrt((errMean_peak_2[i]/Mean_peak_2[i])*(errMean_peak_2[i]/Mean_peak_2[i])+(errMean_peak_1[i]/Mean_peak_1[i])*(errMean_peak_1[i]/Mean_peak_1[i]));

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
 *
 *
 *
 *
 ******************************************************************************/
