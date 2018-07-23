// DCR_CrossTalk.cxx

/********************************************************************************
 *  DCR_CrossTalk.cxx                                                           *
 *                                                                              *
 *  Plots DCR and CrossTalk for the following SiPMs:                            *
 *       FBK HD3-2 # 1                                                          *
 *       FBK HD3-2 # 2                                                          *
 *       FBK HD3-2 # 3                                                          *
 *                                                                              *
 *  To run:                                                                     *
 *    $ root -l                                                                 *
 *    root[0] .x DCR_CrossTalk.cxx                                              *
 *                                                                              *
 *  Davide Depaoli 2018                                                         *
 *                                                                              *
 ********************************************************************************/

void DCR_CrossTalk(){
    //----------------------------------------------------
    //---------------[   SiPM 1 (HD3_2)   ]---------------
    //----------------------------------------------------
    double HV1[] = {34.00, 35.00, 36.00};
    double errHV1[] = {0.01, 0.01, 0.01};
//--------------------------------------------------------------------------------
//--------------------------------------------------------------------------------
//  from files 20180221_HD3-2_1_DARK_34_AS_2_01.txt, 20180221_HD3-2_1_DARK_35_AS_2_01.txt, 20180221_HD3-2_1_DARK_36_AS_2_01.txt:

    //FROM PEAKS NUM (DCR = #peaks/time)
//         double DCR1[] = {10.4116, 11.7785, 13.0908};
//         double errDCR1[] = {0.0, 0.0, 0.0};
//         double CrossTalk1[] = {0.324796, 0.366729, 0.4047};
//         double errCrossTalk1[] = {0., 0., 0.};
    //FROM EXP FIT (exp fit of delays distribution)
//         double DCR1[] =         {14.145, 16.1991, 18.123};
//         double errDCR1[] =      {0.149264, 0.142905, 0.141507};
//         double CrossTalk1[] =   {0.319225, 0.381651, 0.438408};
//         double errCrossTalk1[] ={0.0216727, 0.0158506, 0.0121084};

//--------------------------------------------------------------------------------
//--------------------------------------------------------------------------------
//  from files 20180305_HD3-2_1_DARK_34_AS_2_01, 20180305_HD3-2_1_DARK_35_AS_2_01, 20180305_HD3-2_1_DARK_36_AS_2_01, setting automatically pe_0_5 and pe_1_5:
//     double DCR1[] =         {15.4513, 17.2266, 18.7854};
//     double errDCR1[] =      {0.188006, 0.182623, 0.179374};
//     double CrossTalk1[] =   {0.341012, 0.359794, 0.402453};
//     double errCrossTalk1[] ={0.0263061, 0.0194315, 0.0155378};

//  from files 20180305_HD3-2_1_DARK_34_AS_2_01, 20180305_HD3-2_1_DARK_35_AS_2_01, 20180305_HD3-2_1_DARK_36_AS_2_01, setting manually pe_0_5 and pe_1_5:
    double DCR1[] =         {15.4943, 17.5239, 19.4647};
    double errDCR1[] =      {0.0596534, 0.0581081, 0.0571528};
    double CrossTalk1[] =   {0.332033, 0.358656, 0.411316};
    double errCrossTalk1[] ={0.00822453, 0.00584069, 0.00433039};




//********************************************************************************

    //----------------------------------------------------
    //---------------[   SiPM 2 (HD3_2)   ]---------------
    //----------------------------------------------------
    double HV2[] = {34.00, 35.00, 36.00};
    double errHV2[] = {0.01, 0.01, 0.01};
//--------------------------------------------------------------------------------
//--------------------------------------------------------------------------------
//  from files 20180221_HD3-2_2_DARK_34_AS_2_02.txt, 20180221_HD3-2_2_DARK_35_AS_2_02.txt, 20180221_HD3-2_2_DARK_36_AS_2_02.txt
    //FROM PEAKS NUM (DCR = #peaks/time)
//         double DCR2[] = {12.3351, 13.4354, 14.5025};
//         double errDCR2[] = {0.0, 0.0, 0.0};
//         double CrossTalk2[] ={0.33728, 0.384085, 0.42422};
//         double errCrossTalk2[] = {0., 0., 0.};
    //FROM EXP FIT (exp fit of delays distribution)
//         double DCR2[] =         {17.0027, 18.7497, 20.4274};
//         double errDCR2[] =      {0.141721, 0.139356, 0.14036};
//         double CrossTalk2[] =   {0.360867, 0.381071, 0.40627};
//         double errCrossTalk2[] ={0.0152721, 0.0120017, 0.00990401};

//--------------------------------------------------------------------------------
//--------------------------------------------------------------------------------
//  from files 20180305_HD3-2_2_DARK_34_AS_2_01, 20180305_HD3-2_2_DARK_35_AS_2_01, 20180305_HD3-2_2_DARK_36_AS_2_01, setting automatically pe_0_5 and pe_1_5::
//     double DCR2[] =         {15.5559, 17.9235, 20.4624};
//     double errDCR2[] =      {0.188033, 0.181892, 0.180436};
//     double CrossTalk2[] =   {0.316189, 0.36327, 0.367274};
//     double errCrossTalk2[] ={0.0260671, 0.0182476, 0.0137298};

//  from files 20180305_HD3-2_2_DARK_34_AS_2_01, 20180305_HD3-2_2_DARK_35_AS_2_01, 20180305_HD3-2_2_DARK_36_AS_2_01, setting manually pe_0_5 and pe_1_5::
    double DCR2[] =         {15.4412, 17.8785, 20.1119};
    double errDCR2[] =      {0.0595269, 0.0576375, 0.0569482};
    double CrossTalk2[] =   {0.335797, 0.366711, 0.375398};
    double errCrossTalk2[] ={0.00800791, 0.0057306, 0.0043311};



//********************************************************************************

    //----------------------------------------------------
    //---------------[   SiPM 3 (HD3_2)   ]---------------
    //----------------------------------------------------

//--------------------------------------------------------------------------------
//--------------------------------------------------------------------------------
//  from files 20180221_HD3-2_3_DARK_34_AS_2_01.txt, 20180221_HD3-2_3_DARK_35_AS_2_01.txt, 20180221_HD3-2_3_DARK_36_AS_2_01.txt
    double HV3[] = {34.00, 35.00, 36.00};
    double errHV3[] = {0.01, 0.01, 0.01};
    //FROM PEAKS NUM (DCR = #peaks/time)
//         double DCR3[] = {13.2867, 14.7243, 16.0225};
//         double errDCR3[] = {0.0, 0.0, 0.0};
//         double CrossTalk3[] = {0.351083, 0.402813, 0.444291};
//         double errCrossTalk3[] = {0., 0., 0.};
    //FROM EXP FIT (exp fit of delays distribution)
//         double DCR3[] =         {18.5443, 20.8988, 23.2489};
//         double errDCR3[] =      {0.140562, 0.138892, 0.140885};
//         double CrossTalk3[] =   {0.35188, 0.377134, 0.408868};
//         double errCrossTalk3[] ={0.0130153, 0.00995042, 0.00793973};
//--------------------------------------------------------------------------------
//--------------------------------------------------------------------------------
//  from files 20180305_HD3-2_3_DARK_34_AS_2_01, 20180305_HD3-2_3_DARK_35_AS_2_01, 20180305_HD3-2_3_DARK_36_AS_2_01, setting automatically pe_0_5 and pe_1_5:
//     double DCR3[] =         {18.7002, 21.3137, 23.7438};
//     double errDCR3[] =      {0.180412, 0.178593, 0.180637};
//     double CrossTalk3[] =   {0.313792, 0.343309, 0.362197};
//     double errCrossTalk3[] ={0.0183106, 0.0134591, 0.0106564};

//  from files 20180305_HD3-2_3_DARK_34_AS_2_01, 20180305_HD3-2_3_DARK_35_AS_2_01, 20180305_HD3-2_3_DARK_36_AS_2_01, setting manually pe_0_5 and pe_1_5:
    double DCR3[] =         {18.5218, 20.9923, 23.328};
    double errDCR3[] =      {0.0570749, 0.0563976, 0.0570021};
    double CrossTalk3[] =   {0.317707, 0.346962, 0.375962};
    double errCrossTalk3[] ={0.00579663, 0.00424867, 0.00324151};


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
