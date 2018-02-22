//DCR_CrossTalk.C

void DCR_CrossTalk(){
    /* SiPM 1 (HD3_2)
     *      20180221_HD3-2_1_DARK_34_AS_2_01.txt
     *      20180221_HD3-2_1_DARK_35_AS_2_01.txt
     *      20180221_HD3-2_1_DARK_36_AS_2_01.txt
     */
    double HV1[] = {34.00, 35.00, 36.00};
    double errHV1[] = {0.01, 0.01, 0.01};
    double DCR1[] = {10.4116, 11.7785, 13.0908};
    double errDCR1[] = {0.0, 0.0, 0.0};
    double CrossTalk1[] = {0.324796, 0.366729, 0.4047};
    double errCrossTalk1[] = {0., 0., 0.};
    
    /* SiPM 2 (HD3_2)
     *      20180221_HD3-2_2_DARK_34_AS_2_02.txt
     *      20180221_HD3-2_2_DARK_35_AS_2_02.txt
     *      20180221_HD3-2_2_DARK_36_AS_2_02.txt
     */
    double HV2[] = {34.00, 35.00, 36.00};
    double errHV2[] = {0.01, 0.01, 0.01};
    double DCR2[] = {12.3351, 13.4354, 14.5025};
    double errDCR2[] = {0.0, 0.0, 0.0};
    double CrossTalk2[] ={0.33728, 0.384085, 0.42422};
    double errCrossTalk2[] = {0., 0., 0.};
    
    /* SiPM 3 (HD3_2)
     *      20180221_HD3-2_3_DARK_34_AS_2_01.txt
     *      20180221_HD3-2_3_DARK_35_AS_2_01.txt
     *      20180221_HD3-2_3_DARK_36_AS_2_01.txt
     */
    double HV3[] = {34.00, 35.00, 36.00};
    double errHV3[] = {0.01, 0.01, 0.01};
    double DCR3[] = {13.2867, 14.7243, 16.0225};
    double errDCR3[] = {0.0, 0.0, 0.0};
    double CrossTalk3[] = {0.351083, 0.402813, 0.444291};
    double errCrossTalk3[] = {0., 0., 0.};
    
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
