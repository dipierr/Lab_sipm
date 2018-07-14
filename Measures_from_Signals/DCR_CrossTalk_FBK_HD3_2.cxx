//DCR_CrossTalk.cxx

/********************************************************************************
 *  DCR_CrossTalk.cxx                                                           *
 *                                                                              *
 *  Plots DCR and CrossTalk for the following SiPMs:                            *
 *       FBK HD3-2 # 1                                                          *
 *       FBK HD3-2 # 2                                                          *
 *       FBK HD3-2 # 3                                                          *
 *                                                                              *
 *  > from Ana_Traces_SiPM.cxx, 14/07/2018                                      *
 *  > pe_0_5 and pe_1_5 set manually                                            *
 *  > smoot_trace_3()                                                           *
 *                                                                              *
 *  To run:                                                                     *
 *    $ root -l                                                                 *
 *    root[0] .x DCR_CrossTalk.cxx                                              *
 *                                                                              *
 *  Davide Depaoli 2018                                                         *
 *                                                                              *
 ********************************************************************************/


void DCR_CrossTalk_FBK_HD3_2(){
    //----------------------------------------------------
    //---------------[   SiPM 1 (HD3_2)   ]---------------
    //----------------------------------------------------
    double HV1[] = {34.00, 35.00, 36.00};
    double errHV1[] = {0.01, 0.01, 0.01};
    // Files analyzed:
    // ../../Data_DRS/20180305/20180305_HD3-2_1_DARK_34_AS_2_01   pe_0_5 = 10 mV; pe_1_5 = 25 mV
    // ../../Data_DRS/20180305/20180305_HD3-2_1_DARK_35_AS_2_01   pe_0_5 = 10 mV; pe_1_5 = 28 mV
    // ../../Data_DRS/20180305/20180305_HD3-2_1_DARK_36_AS_2_01   pe_0_5 = 10 mV; pe_1_5 = 30 mV
    double DCR1[] =         {15.7224, 17.7358, 19.6714};
    double errDCR1[] =      {0.061682, 0.0602938, 0.0595938};
    double CrossTalk1[] =   {0.327841, 0.333764, 0.3734};
    double errCrossTalk1[] ={0.00863639, 0.00625938, 0.00488089};



//********************************************************************************

    //----------------------------------------------------
    //---------------[   SiPM 2 (HD3_2)   ]---------------
    //----------------------------------------------------
    double HV2[] = {34.00, 35.00, 36.00};
    double errHV2[] = {0.01, 0.01, 0.01};
    // Files analyzed:
    // ../../Data_DRS/20180305/20180305_HD3-2_2_DARK_34_AS_2_01   pe_0_5 = 10 mV; pe_1_5 = 25 mV
    // ../../Data_DRS/20180305/20180305_HD3-2_2_DARK_35_AS_2_01   pe_0_5 = 10 mV; pe_1_5 = 28 mV
    // ../../Data_DRS/20180305/20180305_HD3-2_2_DARK_36_AS_2_01   pe_0_5 = 10 mV; pe_1_5 = 30 mV
    double DCR2[] =         {15.6179, 18.0689, 20.1908};
    double errDCR2[] =      {0.0618584, 0.0599533, 0.0596247};
    double CrossTalk2[] =   {0.304479, 0.336474, 0.358623};
    double errCrossTalk2[] ={0.0086457, 0.00607813, 0.00465934};






//********************************************************************************

    //----------------------------------------------------
    //---------------[   SiPM 3 (HD3_2)   ]---------------
    //----------------------------------------------------
    double HV3[] = {34.00, 35.00, 36.00};
    double errHV3[] = {0.01, 0.01, 0.01};

    // Files analyzed:
    // ../../Data_DRS/20180305/20180305_HD3-2_3_DARK_34_AS_2_01   pe_0_5 = 10 mV; pe_1_5 = 25 mV
    // ../../Data_DRS/20180305/20180305_HD3-2_3_DARK_35_AS_2_01   pe_0_5 = 10 mV; pe_1_5 = 28 mV
    // ../../Data_DRS/20180305/20180305_HD3-2_3_DARK_36_AS_2_01   pe_0_5 = 10 mV; pe_1_5 = 30 mV
    double DCR3[] =         {18.7316, 21.1419, 23.6828};
    double errDCR3[] =      {0.0592112, 0.0589229, 0.060225};
    double CrossTalk3[] =   {0.309768, 0.328183, 0.354642};
    double errCrossTalk3[] ={0.00606509, 0.00451975, 0.00354651};


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
