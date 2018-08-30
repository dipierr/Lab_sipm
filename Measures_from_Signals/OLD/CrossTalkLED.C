// ****************************************************************************
// *********************************   OLD   **********************************
// ****************************************************************************

/******************************************************************************\
 * CROSS TALK
 *
 * Fast cross talk evaluation from LED measures
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

void CrossTalkLED(){
    // HV
    double HV =               35;
    double errHV =            0.01;

    double H_peak_0 =         4.11691e+02;//4.11503e+02;//2.65356e+02;//2.55733e+02;//2.81738e+02;//4.20934e+02;//4.40787e+02;
    double errH_peak_0 =      1.05600e+01;//1.05713e+01;//7.41084e+00;//8.04486e+00;//9.09474e+00;//1.06258e+01;//1.08467e+01;

    double Sigma_peak_0 =     2.96736e+00;//2.96915e+00;//2.87864e+00;//2.97378e+00;//2.51727e+00;//2.89074e+00;//2.81710e+00;
    double errSigma_peak_0 =  5.47997e-02;//5.50531e-02;//5.22749e-02;//6.48433e-02;//5.91543e-02;//5.32370e-02;//4.97934e-02;

    double H_peak_1 =         6.62548e+02;//6.69509e+02;//5.82884e+02;//5.26792e+02;//5.44120e+02;//7.26267e+02;//7.45052e+02;
    double errH_peak_1 =      1.10700e+01;//1.11124e+01;//8.83242e+00;//8.76612e+00;//1.00994e+01;//1.14091e+01;//1.14332e+01;

    double Sigma_peak_1 =     3.92106e+00;//4.05929e+00;//5.77610e+00;//5.14188e+00;//4.13337e+00;//4.50766e+00;//4.43916e+00;
    double errSigma_peak_1 =  4.99643e-02;//5.30357e-02;//8.61158e-02;//7.60005e-02;//6.19769e-02;//6.35786e-02;//5.81713e-02;

    double Entries =          46735;//47217;//49999;


    //------------------------------



    double errEntries = TMath::Power(Entries,0.5);
    double Area0 = H_peak_0*Sigma_peak_0*TMath::Power(2*TMath::Pi(),0.5)/1;
    double Prob_0pe = Area0/Entries;
    double Mu = -TMath::Log(Prob_0pe);
    double Prob_1pe = Mu*TMath::Exp(-Mu);
    cout<<Prob_1pe<<endl;
    double Area1 = H_peak_1*Sigma_peak_1*TMath::Power(2*TMath::Pi(),0.5)/1;
    double Prob_1peS = Area1/Entries;
    cout<<Prob_1peS<<endl;
    double Prob_Cross_Talk = 1-(Prob_1peS/Prob_1pe);
    cout << "Cross Talk " << Prob_Cross_Talk*100 <<" \%" << endl;


    // double errArea1=Area1*TMath::Power(TMath::Power(errH_peak_1/H_peak_1,2)+TMath::Power(errSigma_peak_1/Sigma_peak_1,2),0.5);
    // double errProb_1peS=errArea1/Entries;
    // double errMu=errArea0/(Prob_0pe*Entries);
    // double errProb_1pe=TMath::Exp(-Mu)*TMath::Abs(1-Mu)*errMu;
    // double errArea0=Area0*TMath::Power(TMath::Power(errH_peak_0/H_peak_0,2)+TMath::Power(errSigma_peak_0/Sigma_peak_0,2),0.5);
    // double errProb_0pe=errArea0/Entries;
    // double errProb_Cross_Talk=TMath::Power(TMath::Power(errProb_1peS/Prob_1pe,2)+TMath::Power(Prob_1peS*errProb_1pe/Prob_1pe,2),0.5);
  }
