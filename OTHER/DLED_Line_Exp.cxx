// DLED_Line_Exp.cxx

/********************************************************************************
 *  DLED_Line_Exp.cxx                                                           *
 *                                                                              *
 *  Plots some examples of signals analyzed using the DLED techniques           *
 *                                                                              *
 *  To RUN:                                                                     *
 *  $ root -l                                                                   *
 *  root[0] .x DLED_Line_Exp.cxx                                                *
 *                                                                              *
 *                                                                              *
 *  Davide Depaoli 2018                                                         *
 *                                                                              *
 ********************************************************************************/


#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <cstdio>
#include <fcntl.h>
#include <unistd.h>
#include <math.h>
#include "TMath.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TAxis.h"
#include "TPad.h"
#include "TLine.h"
#include "TSpectrum.h"
#include "TVirtualFitter.h"

#define timeUp 5
#define timeDown 200
#define k_fall 0.01

#define label_size 0.05

const int tot = timeDown+timeUp+1;


void DLED_Line_Exp(){

    double x[tot], x_shift[tot], y[tot], y_shift[tot], y_DLED[tot], y1[tot], ysum1[tot], ysum1_DLED[tot], y2[tot], ysum2[tot], ysum2_DLED[tot], y3[tot], ysum3[tot], ysum3_DLED[tot],y4[tot], ysum4[tot], ysum4_DLED[tot],y5[tot], ysum5[tot], ysum5_DLED[tot], y_DLED_2[tot], y_DLED_3[tot], y_DLED_4[tot], y_DLED_5[tot];
    double min = 0;
    double max = 100;
    int dleddt = timeUp;
    double dyUp, dyDown;
    int d1 = (int) (0.5*timeUp); //DLED_0_timeUp
    int d2 = timeUp; //DLED_timeUp
    int d3 = (int) (1.5*timeUp); //DLED_timeUp_2timeUp
    int d4 = 2*timeUp; //DLED_2timeUp
    int d5 = 10*timeUp; //DLED_2timeUp_inf
    int k=0;
    int dleddt2 = dleddt-4;
    int dleddt3 = dleddt-2;
    int dleddt4 = dleddt+2;
    int dleddt5 = dleddt+4;

    cout<<"timeUp  = "<<timeUp<<endl;
    cout<<"d1      = "<<d1<<endl;
    cout<<"d2      = "<<d2<<endl;
    cout<<"d3      = "<<d3<<endl;
    cout<<"d4      = "<<d4<<endl;
    cout<<"d5      = "<<d5<<endl;
    cout<<"dleddt  = "<<dleddt<<endl;
    cout<<"dleddt2 = "<<dleddt2<<endl;
    cout<<"dleddt3 = "<<dleddt3<<endl;
    cout<<"dleddt4 = "<<dleddt4<<endl;
    cout<<"dleddt5 = "<<dleddt5<<endl;


    dyUp = (max-min)/timeUp;
    dyDown = (max-min)/timeDown;

    for(int i=0; i<tot; i++){
        x[i] = x_shift[i] =  y[i] = y_DLED[i] = y1[i] = ysum1[i]= ysum1_DLED[i] = y2[i] = ysum2[i]= ysum2_DLED[i] = y3[i] = ysum3[i]= ysum3_DLED[i] = y4[i] = ysum4[i]= ysum4_DLED[i] = y5[i] = ysum5[i]= ysum5_DLED[i] = ysum5_DLED[i] = y_DLED_2[i] = y_DLED_3[i] = y_DLED_4[i] = y_DLED_5[i]= 0.;
    }

    //x
    for(int i=0; i<=timeUp+timeDown; i++){
        x[i] = i;
        x_shift[i] = dleddt+i;
    }

    //y
    for(int i=1; i<=timeUp; i++){
        y[i] = y[i-1]+dyUp;
    }
    // for(int i=timeUp+1; i<=timeUp+timeDown; i++){
    //     x[i] = i;
    //     y[i] = y[i-1]-dyDown;
    // }
    for(int i=timeUp+1; i<=timeUp+timeDown; i++){
        y[i] = y[timeUp] *TMath::Exp(- k_fall * (i - timeUp));
    }

    //yi
    for(int i=d1; i<tot; i++){
        y1[i] = y[i - d1];
    }
    for(int i=d2; i<tot; i++){
        y2[i] = y[i - d2];
    }
    for(int i=d3; i<tot; i++){
        y3[i] = y[i - d3];
    }
    for(int i=d4; i<tot; i++){
        y4[i] = y[i - d4];
    }
    for(int i=d5; i<tot; i++){
        y5[i] = y[i - d5];
    }

    //ysumi
    for(int i=0; i<tot; i++){
        ysum1[i] = y[i] + y1[i];
        ysum2[i] = y[i] + y2[i];
        ysum3[i] = y[i] + y3[i];
        ysum4[i] = y[i] + y4[i];
        ysum5[i] = y[i] + y5[i];
    }

    //y_shift
    k=dleddt;
    for(int i=0; i<tot; i++){
        y_shift[i] = -y[k - dleddt];
        k++;
    }

    //DLED
    for(int i=dleddt; i<tot; i++){
        y_DLED[i] = y[i]-y[i - dleddt];
        ysum1_DLED[i] = ysum1[i]-ysum1[i - dleddt];
        ysum2_DLED[i] = ysum2[i]-ysum2[i - dleddt];
        ysum3_DLED[i] = ysum3[i]-ysum3[i - dleddt];
        ysum4_DLED[i] = ysum4[i]-ysum4[i - dleddt];
        ysum5_DLED[i] = ysum5[i]-ysum5[i - dleddt];
    }

    // DLED with different dleddt
    for(int i=dleddt2; i<tot; i++){
        y_DLED_2[i] = y[i]-y[i - dleddt2];
    }
    for(int i=dleddt3; i<tot; i++){
        y_DLED_3[i] = y[i]-y[i - dleddt3];
    }
    for(int i=dleddt4; i<tot; i++){
        y_DLED_4[i] = y[i]-y[i - dleddt4];
    }
    for(int i=dleddt5; i<tot; i++){
        y_DLED_5[i] = y[i]-y[i - dleddt5];
    }

    //y
    TCanvas *DLED_Original = new TCanvas("DLED_Original","DLED_Original");
    DLED_Original->Divide(1,2);

    TGraph* gr = new TGraph(tot,x,y);
    gr->SetTitle("Signal");
    gr->SetLineColor(kBlue);
    DLED_Original->cd(1);
    gr->Draw("al");

    TGraph* gr_DLED = new TGraph(tot,x,y_DLED);
    gr_DLED->SetTitle("DLED");
    gr_DLED->SetLineColor(kRed);
    DLED_Original->cd(2);
    gr_DLED->Draw("al");


    // different dleddt
    TCanvas *DLED_Different_dleddt = new TCanvas("DLED_Different_dleddt","DLED_Different_dleddt");

    TGraph* gr_DLED_2 = new TGraph(tot,x,y_DLED_2);
    gr_DLED_2->SetLineColor(kAzure+1);

    TGraph* gr_DLED_3 = new TGraph(tot,x,y_DLED_3);
    gr_DLED_3->SetLineColor(kGreen+1);

    TGraph* gr_DLED_4 = new TGraph(tot,x,y_DLED_4);
    gr_DLED_4->SetLineColor(kMagenta);

    TGraph* gr_DLED_5 = new TGraph(tot,x,y_DLED_5);
    gr_DLED_5->SetLineColor(kViolet-1);

    TMultiGraph *gDLED_Different_dleddt = new TMultiGraph("gDLED_Different_dleddt", ";Time (a.u.); Amplitude (a.u.)");
    gDLED_Different_dleddt->Add(gr);
    gDLED_Different_dleddt->Add(gr_DLED);
    gDLED_Different_dleddt->Add(gr_DLED_2);
    gDLED_Different_dleddt->Add(gr_DLED_3);
    gDLED_Different_dleddt->Add(gr_DLED_4);
    gDLED_Different_dleddt->Add(gr_DLED_5);

    DLED_Different_dleddt->cd();
    gDLED_Different_dleddt->Draw("al");


    //y and y_shift
    TCanvas *DLED_Original_Shift = new TCanvas("DLED_Original_Shift","DLED_Original_Shift");
    DLED_Original_Shift->Divide(1,2);

    TGraph* gr_shift = new TGraph(tot,x_shift,y_shift);
    gr_shift->SetLineColor(kGreen+1);
    DLED_Original_Shift->cd(1);
    gr->Draw("al");
    gr_shift->Draw("same");

    DLED_Original_Shift->cd(2);
    gr_DLED->Draw("al");



    //ysum1
    TCanvas *DLED_0_timeUp = new TCanvas("DLED_0_timeUp","DLED_0_timeUp");
    DLED_0_timeUp->Divide(1,2);

    TGraph* gr1 = new TGraph(tot,x,y1);
    gr1->SetTitle("Signal");
    gr1->SetLineColor(kGreen+1);

    TGraph* grsum1 = new TGraph(tot,x,ysum1);
    grsum1->SetTitle("Signal");
    grsum1->SetLineColor(kBlack);
    DLED_0_timeUp->cd(1);
    grsum1->Draw("al");
    gr->Draw("same");
    gr1->Draw("same");

    TGraph* grsum1_DLED = new TGraph(tot,x,ysum1_DLED);
    grsum1_DLED->SetTitle("DLED");
    grsum1_DLED->SetLineColor(kRed);
    DLED_0_timeUp->cd(2);
    grsum1_DLED->Draw("al");

    //ysum2
    TCanvas *DLED_timeUp = new TCanvas("DLED_timeUp","DLED_timeUp");
    DLED_timeUp->Divide(1,2);

    TGraph* gr2 = new TGraph(tot,x,y2);
    gr2->SetTitle("Signal");
    gr2->SetLineColor(kGreen+1);

    TGraph* grsum2 = new TGraph(tot,x,ysum2);
    grsum2->SetTitle("Signal");
    grsum2->SetLineColor(kBlack);
    DLED_timeUp->cd(1);
    grsum2->Draw("al");
    gr->Draw("same");
    gr2->Draw("same");

    TGraph* grsum2_DLED = new TGraph(tot,x,ysum2_DLED);
    grsum2_DLED->SetTitle("DLED");
    grsum2_DLED->SetLineColor(kRed);
    DLED_timeUp->cd(2);
    grsum2_DLED->Draw("al");

    //ysum3
    TCanvas *DLED_timeUp_2timeUp = new TCanvas("DLED_timeUp_2timeUp","DLED_timeUp_2timeUp");
    DLED_timeUp_2timeUp->Divide(1,2);

    TGraph* gr3 = new TGraph(tot,x,y3);
    gr3->SetTitle("Signal");
    gr3->SetLineColor(kGreen+1);

    TGraph* grsum3 = new TGraph(tot,x,ysum3);
    grsum3->SetTitle("Signal");
    grsum3->SetLineColor(kBlack);
    DLED_timeUp_2timeUp->cd(1);
    grsum3->Draw("al");
    gr->Draw("same");
    gr3->Draw("same");

    TGraph* grsum3_DLED = new TGraph(tot,x,ysum3_DLED);
    grsum3_DLED->SetTitle("DLED");
    grsum3_DLED->SetLineColor(kRed);
    DLED_timeUp_2timeUp->cd(2);
    grsum3_DLED->Draw("al");


    //ysum4
    TCanvas *DLED_2timeUp = new TCanvas("DLED_2timeUp","DLED_2timeUp");
    DLED_2timeUp->Divide(1,2);

    TGraph* gr4 = new TGraph(tot,x,y4);
    gr4->SetTitle("Signal");
    gr4->SetLineColor(kGreen+1);

    TGraph* grsum4 = new TGraph(tot,x,ysum4);
    grsum4->SetTitle("Signal");
    grsum4->SetLineColor(kBlack);
    DLED_2timeUp->cd(1);
    grsum4->Draw("al");
    gr->Draw("same");
    gr4->Draw("same");

    TGraph* grsum4_DLED = new TGraph(tot,x,ysum4_DLED);
    grsum4_DLED->SetTitle("DLED");
    grsum4_DLED->SetLineColor(kRed);
    DLED_2timeUp->cd(2);
    grsum4_DLED->Draw("al");


    //ysum5
    TCanvas *DLED_2timeUp_inf = new TCanvas("DLED_2timeUp_inf","DLED_2timeUp_inf");
    DLED_2timeUp_inf->Divide(1,2);

    TGraph* gr5 = new TGraph(tot,x,y5);
    gr5->SetTitle("Signal");
    gr5->SetLineColor(kGreen+1);

    TGraph* grsum5 = new TGraph(tot,x,ysum5);
    grsum5->SetTitle("Signal");
    grsum5->SetLineColor(kBlack);
    DLED_2timeUp_inf->cd(1);
    grsum5->Draw("al");
    gr->Draw("same");
    gr5->Draw("same");

    TGraph* grsum5_DLED = new TGraph(tot,x,ysum5_DLED);
    grsum5_DLED->SetTitle("DLED");
    grsum5_DLED->SetLineColor(kRed);
    DLED_2timeUp_inf->cd(2);
    grsum5_DLED->Draw("al");

    //////////////////////
    ///   AXIS TITLE   ///
    //////////////////////

    gr->GetXaxis()->SetTitle("Time (a.u.)");
    gr1->GetXaxis()->SetTitle("Time (a.u.)");
    gr2->GetXaxis()->SetTitle("Time (a.u.)");
    gr3->GetXaxis()->SetTitle("Time (a.u.)");
    gr4->GetXaxis()->SetTitle("Time (a.u.)");
    gr5->GetXaxis()->SetTitle("Time (a.u.)");

    grsum1->GetXaxis()->SetTitle("Time (a.u.)");
    grsum2->GetXaxis()->SetTitle("Time (a.u.)");
    grsum3->GetXaxis()->SetTitle("Time (a.u.)");
    grsum4->GetXaxis()->SetTitle("Time (a.u.)");
    grsum5->GetXaxis()->SetTitle("Time (a.u.)");

    gr_DLED->GetXaxis()->SetTitle("Time (a.u.)");
    gr_DLED_2->GetXaxis()->SetTitle("Time (a.u.)");
    gr_DLED_3->GetXaxis()->SetTitle("Time (a.u.)");
    grsum1_DLED->GetXaxis()->SetTitle("Time (a.u.)");
    grsum2_DLED->GetXaxis()->SetTitle("Time (a.u.)");
    grsum3_DLED->GetXaxis()->SetTitle("Time (a.u.)");
    grsum4_DLED->GetXaxis()->SetTitle("Time (a.u.)");
    grsum5_DLED->GetXaxis()->SetTitle("Time (a.u.)");


    gr->GetYaxis()->SetTitle("Amplitude (a.u.)");
    gr1->GetYaxis()->SetTitle("Amplitude (a.u.)");
    gr2->GetYaxis()->SetTitle("Amplitude (a.u.)");
    gr3->GetYaxis()->SetTitle("Amplitude (a.u.)");
    gr4->GetYaxis()->SetTitle("Amplitude (a.u.)");
    gr5->GetYaxis()->SetTitle("Amplitude (a.u.)");

    grsum1->GetYaxis()->SetTitle("Amplitude (a.u.)");
    grsum2->GetYaxis()->SetTitle("Amplitude (a.u.)");
    grsum3->GetYaxis()->SetTitle("Amplitude (a.u.)");
    grsum4->GetYaxis()->SetTitle("Amplitude (a.u.)");
    grsum5->GetYaxis()->SetTitle("Amplitude (a.u.)");

    grsum1_DLED->GetYaxis()->SetTitle("Amplitude (a.u.)");
    grsum2_DLED->GetYaxis()->SetTitle("Amplitude (a.u.)");
    grsum3_DLED->GetYaxis()->SetTitle("Amplitude (a.u.)");
    grsum4_DLED->GetYaxis()->SetTitle("Amplitude (a.u.)");
    grsum5_DLED->GetYaxis()->SetTitle("Amplitude (a.u.)");

    ///////////////////////////////////////////////////////////////////////////


    //////////////////////////
    ///   AXIS TITLE SIZE  ///
    //////////////////////////



    // gr->GetXaxis()->SetTitleSize(label_size);
    // gr1->GetXaxis()->SetTitleSize(label_size);
    // gr2->GetXaxis()->SetTitleSize(label_size);
    // gr3->GetXaxis()->SetTitleSize(label_size);
    // gr4->GetXaxis()->SetTitleSize(label_size);
    // gr5->GetXaxis()->SetTitleSize(label_size);
    //
    // grsum1->GetXaxis()->SetTitleSize(label_size);
    // grsum2->GetXaxis()->SetTitleSize(label_size);
    // grsum3->GetXaxis()->SetTitleSize(label_size);
    // grsum4->GetXaxis()->SetTitleSize(label_size);
    // grsum5->GetXaxis()->SetTitleSize(label_size);
    //
    // gr_DLED->GetXaxis()->SetTitleSize(label_size);
    // grsum1_DLED->GetXaxis()->SetTitleSize(label_size);
    // grsum2_DLED->GetXaxis()->SetTitleSize(label_size);
    // grsum3_DLED->GetXaxis()->SetTitleSize(label_size);
    // grsum4_DLED->GetXaxis()->SetTitleSize(label_size);
    // grsum5_DLED->GetXaxis()->SetTitleSize(label_size);
    //
    // gr->GetYaxis()->SetTitleSize(label_size);
    // gr1->GetYaxis()->SetTitleSize(label_size);
    // gr2->GetYaxis()->SetTitleSize(label_size);
    // gr3->GetYaxis()->SetTitleSize(label_size);
    // gr4->GetYaxis()->SetTitleSize(label_size);
    // gr5->GetYaxis()->SetTitleSize(label_size);
    //
    // grsum1->GetYaxis()->SetTitleSize(label_size);
    // grsum2->GetYaxis()->SetTitleSize(label_size);
    // grsum3->GetYaxis()->SetTitleSize(label_size);
    // grsum4->GetYaxis()->SetTitleSize(label_size);
    // grsum5->GetYaxis()->SetTitleSize(label_size);
    //
    // grsum1_DLED->GetYaxis()->SetTitleSize(label_size);
    // grsum2_DLED->GetYaxis()->SetTitleSize(label_size);
    // grsum3_DLED->GetYaxis()->SetTitleSize(label_size);
    // grsum4_DLED->GetYaxis()->SetTitleSize(label_size);
    // grsum5_DLED->GetYaxis()->SetTitleSize(label_size);


    ///////////////////////////////////////////////////////////////////////////


    // gr->GetXaxis()->SetLabelSize(label_size);
    // gr1->GetXaxis()->SetLabelSize(label_size);
    // gr2->GetXaxis()->SetLabelSize(label_size);
    // gr3->GetXaxis()->SetLabelSize(label_size);
    // gr4->GetXaxis()->SetLabelSize(label_size);
    // gr5->GetXaxis()->SetLabelSize(label_size);

    // grsum1->GetXaxis()->SetLabelSize(label_size);
    // grsum2->GetXaxis()->SetLabelSize(label_size);
    // grsum3->GetXaxis()->SetLabelSize(label_size);
    // grsum4->GetXaxis()->SetLabelSize(label_size);
    // grsum5->GetXaxis()->SetLabelSize(label_size);

    // gr_DLED->GetXaxis()->SetLabelSize(label_size);
    // grsum1_DLED->GetXaxis()->SetLabelSize(label_size);
    // grsum2_DLED->GetXaxis()->SetLabelSize(label_size);
    // grsum3_DLED->GetXaxis()->SetLabelSize(label_size);
    // grsum4_DLED->GetXaxis()->SetLabelSize(label_size);
    // grsum5_DLED->GetXaxis()->SetLabelSize(label_size);

    // gr->GetYaxis()->SetLabelSize(label_size);
    // gr1->GetYaxis()->SetLabelSize(label_size);
    // gr2->GetYaxis()->SetLabelSize(label_size);
    // gr3->GetYaxis()->SetLabelSize(label_size);
    // gr4->GetYaxis()->SetLabelSize(label_size);
    // gr5->GetYaxis()->SetLabelSize(label_size);

    // grsum1->GetYaxis()->SetLabelSize(label_size);
    // grsum2->GetYaxis()->SetLabelSize(label_size);
    // grsum3->GetYaxis()->SetLabelSize(label_size);
    // grsum4->GetYaxis()->SetLabelSize(label_size);
    // grsum5->GetYaxis()->SetLabelSize(label_size);

    // grsum1_DLED->GetYaxis()->SetLabelSize(label_size);
    // grsum2_DLED->GetYaxis()->SetLabelSize(label_size);
    // grsum3_DLED->GetYaxis()->SetLabelSize(label_size);
    // grsum4_DLED->GetYaxis()->SetLabelSize(label_size);
    // grsum5_DLED->GetYaxis()->SetLabelSize(label_size);

    ///////////////////////////////////////////////////////////////////////////


    grsum1->SetLineWidth(2);
    grsum2->SetLineWidth(2);
    grsum3->SetLineWidth(2);
    grsum4->SetLineWidth(2);
    grsum5->SetLineWidth(2);


    double minx = 0;
    double maxx = timeUp + timeDown;

    double maxy = y[timeUp]+20;
    double miny = -maxy;

    gr->GetXaxis()->SetRangeUser(minx,maxx);
    grsum1->GetXaxis()->SetRangeUser(minx,maxx);
    grsum2->GetXaxis()->SetRangeUser(minx,maxx);
    grsum3->GetXaxis()->SetRangeUser(minx,maxx);
    grsum4->GetXaxis()->SetRangeUser(minx,maxx);
    grsum5->GetXaxis()->SetRangeUser(minx,maxx);
    gr_shift->GetXaxis()->SetRangeUser(minx,maxx);

    gr_DLED->GetXaxis()->SetRangeUser(minx,maxx);
    gr_DLED_2->GetXaxis()->SetRangeUser(minx,maxx);
    gr_DLED_3->GetXaxis()->SetRangeUser(minx,maxx);
    grsum1_DLED->GetXaxis()->SetRangeUser(minx,maxx);
    grsum2_DLED->GetXaxis()->SetRangeUser(minx,maxx);
    grsum3_DLED->GetXaxis()->SetRangeUser(minx,maxx);
    grsum4_DLED->GetXaxis()->SetRangeUser(minx,maxx);
    grsum5_DLED->GetXaxis()->SetRangeUser(minx,maxx);

    gr_shift->GetYaxis()->SetRangeUser(miny,maxy);
    gr->GetYaxis()->SetRangeUser(miny,maxy);


    DLED_Original->cd(1)->SetGridy();
    DLED_Original_Shift->cd(1)->SetGridy();
    DLED_0_timeUp->cd(1)->SetGridy();
    DLED_timeUp->cd(1)->SetGridy();
    DLED_timeUp_2timeUp->cd(1)->SetGridy();
    DLED_2timeUp->cd(1)->SetGridy();
    DLED_2timeUp_inf->cd(1)->SetGridy();
    DLED_Original->cd(2)->SetGridy();
    DLED_Original_Shift->cd(2)->SetGridy();
    DLED_0_timeUp->cd(2)->SetGridy();
    DLED_timeUp->cd(2)->SetGridy();
    DLED_timeUp_2timeUp->cd(2)->SetGridy();
    DLED_2timeUp->cd(2)->SetGridy();
    DLED_2timeUp_inf->cd(2)->SetGridy();
    DLED_Different_dleddt->SetGridy();


    DLED_Original->Update();
    DLED_0_timeUp->Update();
    DLED_timeUp->Update();
    DLED_timeUp_2timeUp->Update();
    DLED_2timeUp->Update();
    DLED_2timeUp_inf->Update();

    DLED_Original_Shift->Update();

    double maxx_diff_dleddt = (int)(maxx*0.375);
    double maxy_diff_dleddt = maxy;
    double miny_diff_dleddt = -20;
    gDLED_Different_dleddt->GetXaxis()->SetRangeUser(minx,maxx_diff_dleddt);
    gDLED_Different_dleddt->GetYaxis()->SetRangeUser(miny_diff_dleddt,maxy_diff_dleddt);
    DLED_Different_dleddt->Update();


}
