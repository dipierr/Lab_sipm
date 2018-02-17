#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <TMath.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TLegend.h>
#include <TStyle.h>
#include "TMath.h"
#include "TGraphErrors.h"
#include "TAxis.h"
#include "TPad.h"
#include <stdlib.h>     //for using the function sleep

using namespace std;
using namespace TMath;

void LineParser(ifstream& file_to_parse, vector<string>& lines);

void Parse(ifstream& file_to_parse,vector<vector<string> >& spreadsheet);


//GLOBAL VARIABLES
double **trace;
double **trace_DLED;
double **trace_AVG;
double *peak;

int bins_Volt = 180;

//===============================================================================
int DLED(int trace_lenght, int dleddt){
    int trace_DLED_lenght = trace_lenght - dleddt;
    //CREATE TRACE DLED
    trace_DLED = new double*[2];
    for(int i = 0; i < 2; i++) {
        trace_DLED[i] = new double[trace_DLED_lenght];
    }
    for(int i=0; i<trace_DLED_lenght; i++){
        trace_DLED[0][i] = trace[0][i + dleddt];
        trace_DLED[1][i] = trace[1][i + dleddt]-trace[1][i];
    }
    return trace_DLED_lenght;
}

//===============================================================================
int find_peak_fix_time(int mintp, int maxtp){
       
    double max = -100;
    int index=0;
    for(int i=mintp; i<maxtp; i++){
        if(trace_DLED[1][i]>max){
            max=trace_DLED[1][i];
            index=i;
        }
    }
    //cout<<index<<endl;
    return index;
}

//===============================================================================
void average_func(int trace_lenght, int tot_ev, bool last){
    int i;
    for (i=0; i<trace_lenght; i++){
        trace_AVG[1][i] = trace_AVG[1][i] + trace[1][i];
        if(last)
            trace_AVG[1][i] = trace_AVG[1][i]/tot_ev;
    }
}

//===============================================================================   
void show_trace(TCanvas* canv, double *x, double *y,int trace_lenght){
    canv->cd();
    TGraphErrors *graph = new TGraphErrors(trace_lenght,x,y,0,0);
    graph->SetTitle("");
    graph->GetXaxis()->SetTitle("Time (s)");
    graph->GetYaxis()->SetTitle("Amplitude (V)");
    graph->GetYaxis()->SetTitleOffset(1.2);
    graph->GetXaxis()->SetTitleOffset(1.2);
    graph->Draw("");
    canv->Update();
}



//=============================================================================== 
//  MAIN FUNCTION
//===============================================================================   
int Analysis(string file){
    gROOT->Reset();
    
    //VARIABLES TO BE CHANGED
    int min_ind_offset = 0;
    int max_ind_offset = 80;
    double noise_level = 0.;
    int mintp = 450; //min_time_peak
    int maxtp = 550; //min_time_peak
    bool dled = true;
    int dleddt = 9;
    double maxy = .2;
    bool display = true;
    bool average = true;
    
    
    int last_event_n = 5000;
    
    ifstream OpenFile (file.c_str());
    
    //Local variables
    char temp[20];
    int n_ev, index;
    bool reading = true;
    bool last_event_flag = false;
    
    int trace_DLED_lenght = 0;
    int trace_lenght = 0;
    
    TCanvas *c = new TCanvas("Trace","Trace");
    c->SetGrid();
    TCanvas *cDLED = new TCanvas("DLED","DLED");
    cDLED->SetGrid();
    TCanvas *cAVG = new TCanvas("AVG","AVG");
    cAVG->SetGrid();
    TCanvas *cHist = new TCanvas("hist_GAIN","hist_GAIN");
    cHist->SetGrid();
    TH1D *ptrHist = new TH1D("hist","hist",bins_Volt,0,maxy);
    
    
    
    n_ev=0;
    while(!OpenFile.eof() and (reading)){ // reading file
        
        if(n_ev%500==0)
            cout<<"Read ev\t"<<n_ev<<endl;
        
        OpenFile>>temp; 	 
        OpenFile>>temp;
        trace_lenght = atoi(temp);
        OpenFile>>temp>>temp;
        if(n_ev==0)
            cout<<"Number points "<<trace_lenght<<endl;
        
        
        //CREATE TRACE
        trace = new double*[2];
        for(int i = 0; i < 2; i++) {
            trace[i] = new double[trace_lenght];
        }
         if(average==true and n_ev==0){
            trace_AVG = new double*[2];
            for(int i = 0; i < 2; i++) {
                trace_AVG[i] = new double[trace_lenght];
            }
        }
        //CREATE PEAK
        double *peak = new double[2];        
        
        //READ TRACE FROM FILE
        for(int i=0; i<trace_lenght; i++){
            OpenFile>>temp;
            trace[0][i] = atof(temp);
            OpenFile>>temp;
            trace[1][i]  = -atof(temp);
        }
        
        //DLED
        if(dled){
            trace_DLED_lenght = DLED(trace_lenght,dleddt);
        }
        
        //AVERAGE
        if(average){
            if(n_ev==0){
                for(int i=0; i< trace_lenght; i++){
                    trace_AVG[0][i]=trace[0][i];
                    trace_AVG[1][i]=0;
                }
            }
            else{
                if(n_ev==last_event_n-1)
                    last_event_flag = true;
            }
            average_func(trace_lenght,last_event_n, last_event_flag);
        }
        
        index = find_peak_fix_time(mintp, maxtp);
        peak[0] = trace[0][index];
        peak[1] = trace[1][index];
        ptrHist->Fill(peak[1]);
        
        
        
        if(display and n_ev==1){
            show_trace(c,trace[0], trace[1], trace_lenght);
            show_trace(cDLED,trace_DLED[0], trace_DLED[1], trace_DLED_lenght);
        }
        
        
        
        if(n_ev==last_event_n-1)
            reading=false;
        
        delete trace;
        n_ev++;
    }
    if(average){
        show_trace(cAVG,trace_AVG[0], trace_AVG[1], trace_lenght);
    }
    cHist->cd();
    ptrHist->Draw();
    
    
        
    return 0;
}


