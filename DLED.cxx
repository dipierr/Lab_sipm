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

#define timeUp 10
#define timeDown 200

const int tot = timeDown+timeUp+1;


void DLED(){

    double x[tot], y[tot], y_DLED[tot], y1[tot], ysum1[tot], ysum1_DLED[tot], y2[tot], ysum2[tot], ysum2_DLED[tot], y3[tot], ysum3[tot], ysum3_DLED[tot];
    double min = 0;
    double max = 100;
    int dleddt = timeUp;
    double dyUp, dyDown;
    int d1 = 2*timeUp;
    int d2 = 0.5*timeUp;
    int d3 = 2*timeUp - 3;
    
    dyUp = (max-min)/timeUp;
    dyDown = (max-min)/timeDown;
    
    for(int i=0; i<tot; i++){ 
        x[i] = y[i] = y_DLED[i] = y1[i] = ysum1[i]= ysum1_DLED[i] = y2[i] = ysum2[i]= ysum2_DLED[i] = y3[i] = ysum3[i]= ysum3_DLED[i] = 0.;
    }
    
    //y
    for(int i=1; i<=timeUp; i++){
        x[i] = i;
        y[i] = y[i-1]+dyUp;
    }
    for(int i=timeUp+1; i<=timeUp+timeDown; i++){
        x[i] = i;
        y[i] = y[i-1]-dyDown;
    }
    
    //y1, y2
    for(int i=d1; i<tot; i++){
        y1[i] = y[i - d1];
    }
    for(int i=d2; i<tot; i++){
        y2[i] = y[i - d2];
    }
    for(int i=d3; i<tot; i++){
        y3[i] = y[i - d3];
    }
    
    //ysum1, ysum2
    for(int i=0; i<tot; i++){
        ysum1[i] = y[i] + y1[i];
        ysum2[i] = y[i] + y2[i];
        ysum3[i] = y[i] + y3[i];
    }
    
    //DLED
    for(int i=dleddt; i<tot; i++){
        y_DLED[i] = y[i]-y[i - dleddt];
        ysum1_DLED[i] = ysum1[i]-ysum1[i - dleddt];
        ysum2_DLED[i] = ysum2[i]-ysum2[i - dleddt];
        ysum3_DLED[i] = ysum3[i]-ysum3[i - dleddt];
    }
    
    //y
    TCanvas *c = new TCanvas("c","Original");
    c->Divide(1,2);
    
    TGraph* gr = new TGraph(tot,x,y);
    gr->SetTitle("Signal");
    gr->SetLineColor(kBlue);
    c->cd(1);
    gr->Draw();
    
    TGraph* gr_DLED = new TGraph(tot,x,y_DLED);
    gr_DLED->SetTitle("DLED");
    gr_DLED->SetLineColor(kRed);
    c->cd(2);
    gr_DLED->Draw();
    
    //ysum1
    TCanvas *c1 = new TCanvas("c1","DLED_NoBlind");
    c1->Divide(1,2);
    
    TGraph* gr1 = new TGraph(tot,x,y1);
    gr1->SetTitle("Signal");
    gr1->SetLineColor(kGreen+1);
    
    TGraph* grsum1 = new TGraph(tot,x,ysum1);
    grsum1->SetTitle("Signal");
    grsum1->SetLineColor(kBlack);
    c1->cd(1);
    grsum1->Draw();
    gr->Draw("same");
    gr1->Draw("same");
    
    TGraph* grsum1_DLED = new TGraph(tot,x,ysum1_DLED);
    grsum1_DLED->SetTitle("DLED");
    grsum1_DLED->SetLineColor(kRed);
    c1->cd(2);
    grsum1_DLED->Draw();
    
    //ysum2
    TCanvas *c2 = new TCanvas("c2","DLED_Blind");
    c2->Divide(1,2);
    
    TGraph* gr2 = new TGraph(tot,x,y2);
    gr2->SetTitle("Signal");
    gr2->SetLineColor(kGreen+1);
    
    TGraph* grsum2 = new TGraph(tot,x,ysum2);
    grsum2->SetTitle("Signal");
    grsum2->SetLineColor(kBlack);
    c2->cd(1);
    grsum2->Draw();
    gr->Draw("same");
    gr2->Draw("same");
    
    TGraph* grsum2_DLED = new TGraph(tot,x,ysum2_DLED);
    grsum2_DLED->SetTitle("DLED");
    grsum2_DLED->SetLineColor(kRed);
    c2->cd(2);
    grsum2_DLED->Draw();
    
     //ysum3
    TCanvas *c3 = new TCanvas("c3","DLED_3");
    c3->Divide(1,2);
    
    TGraph* gr3 = new TGraph(tot,x,y3);
    gr3->SetTitle("Signal");
    gr3->SetLineColor(kGreen+1);
    
    TGraph* grsum3 = new TGraph(tot,x,ysum3);
    grsum3->SetTitle("Signal");
    grsum3->SetLineColor(kBlack);
    c3->cd(1);
    grsum3->Draw();
    gr->Draw("same");
    gr3->Draw("same");
    
    TGraph* grsum3_DLED = new TGraph(tot,x,ysum3_DLED);
    grsum3_DLED->SetTitle("DLED");
    grsum3_DLED->SetLineColor(kRed);
    c3->cd(2);
    grsum3_DLED->Draw();
    
    
    
    
    
}