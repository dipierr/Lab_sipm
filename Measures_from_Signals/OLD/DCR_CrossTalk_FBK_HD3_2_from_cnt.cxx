// DCR_CrossTalk_FBK_HD3_2_from_cnt.cxx

// ****************************************************************************
// *********************************   OLD   **********************************
// ****************************************************************************

/********************************************************************************
 *  DCR_CrossTalk_FBK_HD3_2_from_cnt.cxx                                        *
 *                                                                              *
 *  Plots DCR and CrossTalk for the following SiPMs:                            *
 *       FBK HD3-2 # 1                                                          *
 *       FBK HD3-2 # 2                                                          *
 *       FBK HD3-2 # 3                                                          *
 *                                                                              *
 *  > from Ana_Traces_SiPM.cxx, 23/07/2018                                      *
 *  > dleddt = 8;                                                               *
 *  > pe_0_5 and pe_1_5 set manually                                            *
 *  > smoot_trace_4()                                                           *
 *                                                                              *
 *  To run:                                                                     *
 *    $ root -l                                                                 *
 *    root[0] .x DCR_CrossTalk_FBK_HD3_2_from_cnt.cxx                           *
 *                                                                              *
 *  Davide Depaoli 2018                                                         *
 *                                                                              *
 ********************************************************************************/


void DCR_CrossTalk_FBK_HD3_2_from_cnt(){
    //----------------------------------------------------
    //---------------[   SiPM 1 (HD3_2)   ]---------------
    //----------------------------------------------------
    double HV1[] = {34.00, 35.00, 36.00};
    double errHV1[] = {0.01, 0.01, 0.01};

    // Files analyzed:
    // ../../Data_DRS/20180305/20180305_HD3-2_1_DARK_34_AS_2_01   pe_0_5 = 10 mV; pe_1_5 = 25 mV
    // ../../Data_DRS/20180305/20180305_HD3-2_1_DARK_35_AS_2_01   pe_0_5 = 10 mV; pe_1_5 = 27 mV
    // ../../Data_DRS/20180305/20180305_HD3-2_1_DARK_36_AS_2_01   pe_0_5 = 10 mV; pe_1_5 = 30 mV
    double DCR1[] =          {14.8465, 17.0197, 19.0359};
    double errDCR1[] =       {0.0135578, 0.0147189, 0.0157625};
    double CrossTalk1[] =    {0.264503, 0.298612, 0.327726};
    double errCrossTalk1[] = {0.000497763, 0.00050644, 0.000512896};


//********************************************************************************

    //----------------------------------------------------
    //---------------[   SiPM 2 (HD3_2)   ]---------------
    //----------------------------------------------------
    double HV2[] = {34.00, 35.00, 36.00};
    double errHV2[] = {0.01, 0.01, 0.01};

    // Files analyzed:
    // ../../Data_DRS/20180305/20180305_HD3-2_2_DARK_34_AS_2_01   pe_0_5 = 10 mV; pe_1_5 = 23 mV
    // ../../Data_DRS/20180305/20180305_HD3-2_2_DARK_35_AS_2_01   pe_0_5 = 10 mV; pe_1_5 = 25 mV
    // ../../Data_DRS/20180305/20180305_HD3-2_2_DARK_36_AS_2_01   pe_0_5 = 10 mV; pe_1_5 = 28 mV
    double DCR2[] =          {14.6634, 17.4487, 19.7036};
    double errDCR2[] =       {0.013458, 0.0149434, 0.0161022};
    double CrossTalk2[] =    {0.272783, 0.304457, 0.332619};
    double errCrossTalk2[] = {0.000510578, 0.000507324, 0.000510367};


//********************************************************************************

    //----------------------------------------------------
    //---------------[   SiPM 3 (HD3_2)   ]---------------
    //----------------------------------------------------
    double HV3[] = {34.00, 35.00, 36.00};
    double errHV3[] = {0.01, 0.01, 0.01};

    // Files analyzed:
    // ../../Data_DRS/20180305/20180305_HD3-2_3_DARK_34_AS_2_01   pe_0_5 = 10 mV; pe_1_5 = 25 mV
    // ../../Data_DRS/20180305/20180305_HD3-2_3_DARK_35_AS_2_01   pe_0_5 = 10 mV; pe_1_5 = 27 mV
    // ../../Data_DRS/20180305/20180305_HD3-2_3_DARK_36_AS_2_01   pe_0_5 = 10 mV; pe_1_5 = 30 mV
    double DCR3[] =          {18.0903, 20.7759, 23.2972};
    double errDCR3[] =       {0.0152766, 0.0166421, 0.017888};
    double CrossTalk3[] =    {0.270357, 0.304444, 0.333352};
    double errCrossTalk3[] = {0.000461877, 0.000470001, 0.00047586};


    //------------------------------

    TGraphErrors *gDCR_1 = new TGraphErrors(3, HV1, DCR1,errHV1, errDCR1);
    TGraphErrors *gDCR_2 = new TGraphErrors(3, HV2, DCR2,errHV2, errDCR2);
    TGraphErrors *gDCR_3 = new TGraphErrors(3, HV3, DCR3,errHV3, errDCR3);

    TGraphErrors *gCT_1 = new TGraphErrors(3, HV1, CrossTalk1,errHV1, errCrossTalk1);
    TGraphErrors *gCT_2 = new TGraphErrors(3, HV2, CrossTalk2,errHV2, errCrossTalk2);
    TGraphErrors *gCT_3 = new TGraphErrors(3, HV3, CrossTalk3,errHV3, errCrossTalk3);

    //------------------------------

    gDCR_1->SetMarkerStyle(20);
    gDCR_2->SetMarkerStyle(20);
    gDCR_3->SetMarkerStyle(20);

    gCT_1->SetMarkerStyle(20);
    gCT_2->SetMarkerStyle(20);
    gCT_3->SetMarkerStyle(20);

    gDCR_1->SetMarkerColor(kOrange+1);
    gDCR_2->SetMarkerColor(kRed);
    gDCR_3->SetMarkerColor(kMagenta);

    gCT_1->SetMarkerColor(kOrange+1);
    gCT_2->SetMarkerColor(kRed);
    gCT_3->SetMarkerColor(kMagenta);

    gDCR_1->SetMarkerSize(3);
    gDCR_2->SetMarkerSize(3);
    gDCR_3->SetMarkerSize(3);

    gCT_1->SetMarkerSize(3);
    gCT_2->SetMarkerSize(3);
    gCT_3->SetMarkerSize(3);

    //------------------------------

    TMultiGraph *DCR = new TMultiGraph("DCR", ";HV (V); DCR (MHz)");
    DCR->Add(gDCR_1);
    DCR->Add(gDCR_2);
    DCR->Add(gDCR_3);

    TMultiGraph *CT = new TMultiGraph("CrossTalk", ";HV (V); Cross Talk");
    CT->Add(gCT_1);
    CT->Add(gCT_2);
    CT->Add(gCT_3);

    //------------------------------

    TCanvas *cDCR = new TCanvas("DCR", "DCR");
    cDCR->SetGrid();
    DCR->Draw("AP");
    auto legendDCR = new TLegend(0.15,0.75,0.35,0.9);
    legendDCR->AddEntry(gDCR_1,"HD3_2 (1)","p");
    legendDCR->AddEntry(gDCR_2,"HD3_2 (2)","p");
    legendDCR->AddEntry(gDCR_3,"HD3_2 (3)","p");
    legendDCR->Draw();

    TCanvas *cCT = new TCanvas("CrossTalk", "CrossTalk");
    cCT->SetGrid();
    CT->Draw("AP");
    auto legendCT = new TLegend(0.15,0.75,0.35,0.9);
    legendCT->AddEntry(gDCR_1,"HD3_2 (1)","p");
    legendCT->AddEntry(gDCR_2,"HD3_2 (2)","p");
    legendCT->AddEntry(gDCR_3,"HD3_2 (3)","p");
    legendCT->Draw();

}
