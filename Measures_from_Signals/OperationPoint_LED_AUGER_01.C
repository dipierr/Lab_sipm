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
 *    > FUNC GEN:   IMPULSER Auger Home Made
 *    > DRIVER:     DRIVER Auger Home Made
 *
 * File obtained using Ana_LED(...) function in Ana_Traces_SiPM.cxx
 *
 * 20180522_HD3-2_1_LED_5583_AGILENT_36_AS_2_50000_04.dat
 * 20180523_HD3-2_1_LED_5583_AGILENT_35.5_AS_2_50000_01.dat
 * 20180511_HD3-2_1_LED_5583_AGILENT_35_AS_2_01-50000eve.dat
 * 20180523_HD3-2_1_LED_5583_AGILENT_34.5_AS_2_50000_01.dat
 * 20180522_HD3-2_1_LED_5583_AGILENT_34_AS_2_50000_01.dat
 * 20180522_HD3-2_1_LED_5583_AGILENT_33_AS_2_50000_02.dat
 * 20180522_HD3-2_1_LED_5583_AGILENT_32_AS_2_50000_03.dat
 * 20180528_HD3-2_1_LED_5583_AGILENT_31_AS_2_50000_01.dat
 * 20180528_HD3-2_1_LED_5583_AGILENT_30_AS_2_50000_01.dat -> No distinction beetween peaks
 * 20180528_HD3-2_1_LED_5583_AGILENT_29_AS_2_50000_01.dat -> No distinction beetween peaks
 *
\******************************************************************************/

#define n_GAIN 8    // for HV = 31 ... 36 V
#define n_MEAN 10   // for HV = 29 ... 36 V


#define h 600
#define w 800

void OperationPoint_LED_AUGER_01(){
    // HV
    double HV[] =  {29.00,30.00,    31.00,        32.00,        33.00,        34.00,        34.50,        35.00,        35.50,        36.00};
    double errHV[]={0.01, 0.01,     0.01,         0.01,         0.01,         0.01,         0.01,         0.01,         0.01,         0.01};


    // PEAK 0
    // H_peak_0
    double H_peak_0[] =            {6.33035e+02,  4.62060e+02,  3.69550e+02,  2.88937e+02,  2.53373e+02,  2.41934e+02,  1.90727e+02,  1.63868e+02};
    double errH_peak_0[] =         {1.11321e+01,  1.17808e+01,  9.97152e+00,  8.50633e+00,  1.00300e+01,  7.27642e+00,  7.84682e+00,  5.80190e+00};

    // Sigma_peak_0
    double Sigma_peak_0[] =        {1.74381e+00,  2.15692e+00,  2.18775e+00,  2.35603,      2.26689e+00,  3.37690e+00,  2.52442e+00,  3.15461e+00};
    double errSigma_peak_0[] =     {5.65270e-02,  3.70254e-02,  5.55387e-02,  7.16945e-02,  6.32656e-02,  7.33161e-02,  7.36245e-02,  1.07703e-01};

    // PEAK 1
    // H_peak_1
    double H_peak_1[] =            {7.34431e+02,  9.16189e+02,  4.44348e+02,  3.83262e+02,  3.53597e+02,  2.05660e+02,  6.62548e+02,  1.83551e+02};
    double errH_peak_1[] =         {1.26639e+01,  1.55551e+01,  9.87949e+00,  8.95463e+00,  8.27998e+00,  6.09070e+00,  6.88303e+00,  4.89531e+00};

    // Sigma_peak_1
    double Sigma_peak_1[] =        {6.44109e+00,  3.02227e+00,  5.83062e+00,  5.60473e+00,  5.13521e+00,  4.49008e+00,  5.95813e+00,  8.44510e+00};
    double errSigma_peak_1[] =     {3.08318e-01,  5.26817e-02,  1.53347e-01,  1.22846e-01,  8.78786e-02,  1.27376e-01,  1.44267e-01,  3.08829e-01};

    // Mean_peak_1
    double Mean_peak_1[] =         {1.49911e+01,  1.54211e+01,  1.75317e+01,  1.88485e+01,  2.01527e+01,  2.13082e+01,  2.16402e+01,  2.29223e+01};
    double errMean_peak_1[] =      {3.37541e-01,  5.15472e-02,  1.22079e-01,  1.09392e-01,  8.35216e-02,  1.19988e-01,  1.43441e-01,  3.58637e-01};

    // PEAK 2
    // Mean_peak_2
    double Mean_peak_2[] =         {2.32595e+01,  2.83134e+01,  3.11680e+01,  3.60932e+01,  3.88570e+01,  4.23323e+01,  4.34901e+01,  4.48226e+01};
    double errMean_peak_2[] =      {1.97267e-01,  8.83312e-02,  2.54707e-01,  1.64642e-01,  1.62696e-01,  2.51139e-01,  3.00476e-01,  6.91813e-01};

    // GLOBAL HIST
    // Mean hist global
    double Mean_hg[]=  {13.37,21.1, 31.28,        42.27,        51.58,        65.41,        73.94,        80.43,        89.72,        98.15};
    double errMean_hg[]={0.03227,0.05395,0.08042, 0.1117,       0.1441,       0.1832,       0.2069,       0.229,        0.254,        0.2756};

    // Standard_dev_hist_global
    double Std_hg[] =  {7.217,12.06,17.98,        24.58,         32.22,        40.97,        46.26,        51.21,        56.8,        69.37};
    double errStd_hg[]={0.02282,0.03815,0.05687,  0.07896,      0.1019,       0.1295,       0.1463,       0.1619,       0.1796,      0.1949};

     // Integral of the average
    double Integral[] =            {-7107.6,      -8169.74,     -9689.57,     -11058.5,     -12327.4,     -13178.8,     -14144.1,     -15025.8};
    double errIntegral[] =         {0.0,          0.0,          0.0,          0.0,          0.0,          0.0,          0.0,          0.0};

    // Entries
    double Entries[] =              {41288,       48452,        39756,        32634,        40183,        29202,        39799,        37981};


     // Cross Talk
    // double Cross_Talk[] =           {0.265269,    0.293336,   0.3389,   0.352329,   0.383516,   0.4243,   0.48297};
    // double errCross_Talk[] =        {0.0229779,   0.010678,   0.00689956,   0.00572358,   0.00477179,   0.00413051,   0.00367729};


    // HV = 31 V
    H_peak_0[0]		= 548.556;
    errH_peak_0[0]		= 12.7617;
    Sigma_peak_0[0]		= 2.19856;
    errSigma_peak_0[0]	= 0.0461278;
    H_peak_1[0]		= 967.574;
    errH_peak_1[0]		= 15.4229;
    Sigma_peak_1[0]		= 3.3765;
    errSigma_peak_1[0]	= 0.0809846;
    Mean_peak_1[0]		= 13.754;
    errMean_peak_1[0]	= 0.0836035;
    Mean_peak_2[0]		= 23.8791;
    errMean_peak_2[0]	= 0.152962;
    Mean_hg[0]		= 32.6272;
    errMean_hg[0]		= 0.0870777;
    Std_hg[0]		= 18.5324;
    errStd_hg[0]		= 0.0615732;
    Entries[0]		= 45295;


    // HV = 32 V
    H_peak_0[1]		= 483.012;
    errH_peak_0[1]		= 11.8118;
    Sigma_peak_0[1]		= 2.28763;
    errSigma_peak_0[1]	= 0.039793;
    H_peak_1[1]		= 682.096;
    errH_peak_1[1]		= 13.1632;
    Sigma_peak_1[1]		= 3.11522;
    errSigma_peak_1[1]	= 0.0700661;
    Mean_peak_1[1]		= 15.2621;
    errMean_peak_1[1]	= 0.0738854;
    Mean_peak_2[1]		= 28.0867;
    errMean_peak_2[1]	= 0.134874;
    Mean_hg[1]		= 42.3386;
    errMean_hg[1]		= 0.12626;
    Std_hg[1]		= 25.3506;
    errStd_hg[1]		= 0.0892792;
    Entries[1]		= 40313;


    // HV = 33 V
    H_peak_0[2]		= 386.91;
    errH_peak_0[2]		= 10.5396;
    Sigma_peak_0[2]		= 2.37543;
    errSigma_peak_0[2]	= 0.0446958;
    H_peak_1[2]		= 534.224;
    errH_peak_1[2]		= 10.6512;
    Sigma_peak_1[2]		= 3.53539;
    errSigma_peak_1[2]	= 0.0642718;
    Mean_peak_1[2]		= 17.4766;
    errMean_peak_1[2]	= 0.0693609;
    Mean_peak_2[2]		= 32.6962;
    errMean_peak_2[2]	= 0.109185;
    Mean_hg[2]		= 55.5792;
    errMean_hg[2]		= 0.164499;
    Std_hg[2]		= 33.6809;
    errStd_hg[2]		= 0.116318;
    Entries[2]		= 41922;


    // HV = 34 V
    H_peak_0[3]		= 324.821;
    errH_peak_0[3]		= 9.88783;
    Sigma_peak_0[3]		= 2.5421;
    errSigma_peak_0[3]	= 0.0592733;
    H_peak_1[3]		= 453.87;
    errH_peak_1[3]		= 9.83674;
    Sigma_peak_1[3]		= 3.64036;
    errSigma_peak_1[3]	= 0.0703767;
    Mean_peak_1[3]		= 19.2893;
    errMean_peak_1[3]	= 0.0727874;
    Mean_peak_2[3]		= 37.4447;
    errMean_peak_2[3]	= 0.137472;
    Mean_hg[3]		= 69.2906;
    errMean_hg[3]		= 0.201358;
    Std_hg[3]		= 42.2329;
    errStd_hg[3]		= 0.142382;
    Entries[3]		= 43991;


    // HV = 34.5 V
    H_peak_0[4]		= 303.46;
    errH_peak_0[4]		= 9.25891;
    Sigma_peak_0[4]		= 2.54088;
    errSigma_peak_0[4]	= 0.0534233;
    H_peak_1[4]		= 454.811;
    errH_peak_1[4]		= 9.94128;
    Sigma_peak_1[4]		= 3.56096;
    errSigma_peak_1[4]	= 0.0638676;
    Mean_peak_1[4]		= 20.5702;
    errMean_peak_1[4]	= 0.0646549;
    Mean_peak_2[4]		= 39.9906;
    errMean_peak_2[4]	= 0.138089;
    Mean_hg[4]		= 77.3584;
    errMean_hg[4]		= 0.224058;
    Std_hg[4]		= 47.4992;
    errStd_hg[4]		= 0.158433;
    Entries[4]		= 44942;


    // HV = 35 V
    H_peak_0[5]		= 249.284;
    errH_peak_0[5]		= 8.45175;
    Sigma_peak_0[5]		= 2.76026;
    errSigma_peak_0[5]	= 0.0691203;
    H_peak_1[5]		= 342.497;
    errH_peak_1[5]		= 8.33875;
    Sigma_peak_1[5]		= 4.01678;
    errSigma_peak_1[5]	= 0.0855678;
    Mean_peak_1[5]		= 21.3667;
    errMean_peak_1[5]	= 0.081144;
    Mean_peak_2[5]		= 42.0344;
    errMean_peak_2[5]	= 0.186486;
    Mean_hg[5]		= 86.2339;
    errMean_hg[5]		= 0.25916;
    Std_hg[5]		= 53.346;
    errStd_hg[5]		= 0.183254;
    Entries[5]		= 42371;

    // HV = 35.5 V
    H_peak_0[6]		= 242.803;
    errH_peak_0[6]		= 8.09327;
    Sigma_peak_0[6]		= 2.79065;
    errSigma_peak_0[6]	= 0.0647025;
    H_peak_1[6]		= 323.747;
    errH_peak_1[6]		= 7.98278;
    Sigma_peak_1[6]		= 4.10116;
    errSigma_peak_1[6]	= 0.0838668;
    Mean_peak_1[6]		= 22.4535;
    errMean_peak_1[6]	= 0.0801244;
    Mean_peak_2[6]		= 44.7902;
    errMean_peak_2[6]	= 0.240814;
    Mean_hg[6]		= 95.2705;
    errMean_hg[6]		= 0.281311;
    Std_hg[6]		= 58.7469;
    errStd_hg[6]		= 0.198917;
    Entries[6]		= 43611;


    // HV = 36 V
    H_peak_0[7]		= 181.508;
    errH_peak_0[7]		= 6.38476;
    Sigma_peak_0[7]		= 3.3331;
    errSigma_peak_0[7]	= 0.0834515;
    H_peak_1[7]		= 206.852;
    errH_peak_1[7]		= 5.94343;
    Sigma_peak_1[7]		= 5.10523;
    errSigma_peak_1[7]	= 0.135809;
    Mean_peak_1[7]		= 23.5037;
    errMean_peak_1[7]	= 0.122086;
    Mean_peak_2[7]		= 46.2825;
    errMean_peak_2[7]	= 0.339862;
    Mean_hg[7]		= 105.433;
    errMean_hg[7]		= 0.335317;
    Std_hg[7]		= 65.5628;
    errStd_hg[7]		= 0.237105;
    Entries[7]		= 38230;





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
