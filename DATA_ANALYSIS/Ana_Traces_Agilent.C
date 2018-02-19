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
#include <TGraphErrors.h>
#include <TAxis.h>
#include <TPad.h>
#include <TLine.h>
#include <stdlib.h>     //for using the function sleep

using namespace std;
using namespace TMath;

void LineParser(ifstream& file_to_parse, vector<string>& lines);

void Parse(ifstream& file_to_parse,vector<vector<string> >& spreadsheet);


//------------------------------------------------------------------------------
//---------------------------[   GLOBAL VARIABLES   ]---------------------------
//------------------------------------------------------------------------------
double **trace;
double **trace_DLED;
double **trace_AVG;
double *peak;

int bins_Volt = 204;

/* VALUES:
 * 
 * LED from Digitizer_CAEN: bins_Volt = 204;
 * 
 */
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------


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
void show_trace(TCanvas* canv, double *x, double *y, int trace_lenght, double miny, double maxy){
    canv->cd();
    for(int i=0; i<trace_lenght; i++){
        x[i] = x[i]*TMath::Power(10,6);
        y[i] = y[i]*TMath::Power(10,3);
    }
    TGraphErrors *graph = new TGraphErrors(trace_lenght,x,y,0,0);
    graph->SetTitle("");
    graph->GetXaxis()->SetTitle("Time (#mus)");
    graph->GetYaxis()->SetTitle("Amplitude (mV)");
    graph->GetYaxis()->SetTitleOffset(1.2);
    graph->GetXaxis()->SetTitleOffset(1.2);
    graph->GetYaxis()-> SetRangeUser(miny,maxy);
    graph->Draw("");
    canv->Update();
}

//===============================================================================   
void show_line(TCanvas *canv, double *x, int mintp, int maxtp, double miny, double maxy){
    canv->cd();
    TLine *lmin = new TLine (x[mintp], miny, x[mintp], maxy);
    TLine *lmax = new TLine (x[maxtp], miny, x[maxtp], maxy);
    lmin->SetLineColor(kBlue);
    lmax->SetLineColor(kBlue);
    lmin->Draw("same");
    lmax->Draw("same");
}



//=============================================================================== 
//  MAIN FUNCTION
//===============================================================================   
int Analysis(string file, int last_event_n){
    gROOT->Reset();
    
//-------------------------------------------------------------------------------
//---------------------------[   SETTING VARIABLES   ]---------------------------
//-------------------------------------------------------------------------------
    int min_ind_offset = 0;
    int max_ind_offset = 80;
    double noise_level = 0.;
    int mintp = 140; //min_time_peak
    int maxtp = 280; //max_time_peak
    int dleddt = 10;
    double maxyhist = .2;
    bool display = true;
    bool average = true;
    bool Agilent_MSO6054A = false; //true if data taken by Agilent_MSO6054A, false otherwise
    bool Digitizer_CAEN = true;  //true if data taken by Digitizer_CAEN, false otherwise
    bool SetLogyHist = false;
    
/* VALUES:
 * 
 * DARK from Agilent: mintp = 160 (or 200);  maxtp = 320; dleddt = 39;
 *                   
 * LED from Agilent: approx in the middle, mintp = 420; maxtp = 500;  dleddt = 9; maxyhist = .2;
 * 
 * LED from Digitizer_CAEN: mintp = 140; maxtp = 280; dleddt = 10;
 * 
 */

//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------

   
    ifstream OpenFile (file.c_str());
    
    //Local variables
    char temp[20];
    int n_ev, index, i;
    bool reading = true;
    bool last_event_flag = false;
    
    int trace_DLED_lenght = 0;
    int trace_lenght = 0;
    
    TCanvas *c = new TCanvas("Trace","Trace");
    c->SetGrid();
    TCanvas *cDLED = new TCanvas("DLED","DLED");
    cDLED->SetGrid();
    TCanvas *cHist = new TCanvas("hist_GAIN","hist_GAIN");
    cHist->SetGrid();
    TH1D *ptrHist = new TH1D("hist","",bins_Volt,0,maxyhist);
    
    double miny=0;
    double maxy=0;
    
    n_ev=0;
    while(!OpenFile.eof() and (reading)){ // reading file
        
        if(n_ev%1000==0)
            cout<<"Read ev\t"<<n_ev<<endl;

//***** READ HEADER FROM FILE
        //select device
        if(Agilent_MSO6054A and not Digitizer_CAEN){
            OpenFile>>temp; 	 
            OpenFile>>temp;
            trace_lenght = atoi(temp);
            OpenFile>>temp>>temp;
        }
        else{
            if(Digitizer_CAEN and not Agilent_MSO6054A){
                OpenFile>>temp>>temp; 
                OpenFile>>temp;
                trace_lenght = atoi(temp);
                OpenFile>>temp>>temp>>temp>>temp>>temp>>temp>>temp>>temp>>temp>>temp>>temp>>temp>>temp>>temp>>temp>>temp>>temp;
        }else{
            cout<<"ERROR: check acquisition device"<<endl;
            return 1;
        }
        }
        if(n_ev==0){
            cout<<"Acquisition device: ";
            if(Agilent_MSO6054A) cout<<"Agilent_MSO6054A"<<endl;
            if(Digitizer_CAEN)   cout<<"Digitizer_CAEN"<<endl;
            cout<<"Number points "<<trace_lenght<<endl;
        }
        
        //cout<<n_ev<<"\t"<<trace_lenght<<endl;
        
        
//***** CREATE TRACE
        //trace
        trace = new double*[2];
        for(i = 0; i < 2; i++) {
            trace[i] = new double[trace_lenght];
        }
        //trace_AVG
        if(average==true and n_ev==0){
            trace_AVG = new double*[2];
            for(i = 0; i < 2; i++) {
                trace_AVG[i] = new double[trace_lenght];
            }
        }

//***** CREATE PEAK
        double *peak = new double[2];        
        
//***** READ TRACE FROM FILE
        if(Agilent_MSO6054A){
                for(i=0; i<trace_lenght; i++){
                    OpenFile>>temp;
                    trace[0][i] = atof(temp);
                    OpenFile>>temp;
                    trace[1][i]  = -atof(temp);
        }
        }
        else{
            if(Digitizer_CAEN){
                for(i=0; i<trace_lenght; i++){
                    trace[0][i] = i*TMath::Power(10,-9);
                    OpenFile>>temp;
                    trace[1][i]  = -atof(temp)/1024;
                }
        }
        }
        
//***** DLED
        trace_DLED_lenght = DLED(trace_lenght,dleddt);
        //Now: trace_DLED
        
//***** AVERAGE
        if(average){
            if(n_ev==0){
                for(i=0; i< trace_lenght; i++){
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
        peak[0] = trace_DLED[0][index];
        peak[1] = trace_DLED[1][index];
        ptrHist->Fill(peak[1]);
        
        
        
        if(display and n_ev == 1){
            if(Agilent_MSO6054A) {miny=10; maxy=180;}
            if(Digitizer_CAEN)   {miny=-900;  maxy=-820;}
            show_trace(c,trace[0], trace[1], trace_lenght,miny,maxy);
            show_line(c, trace[0], mintp, maxtp,miny,maxy);
            if(Agilent_MSO6054A) {miny=-30; maxy=90;}
            if(Digitizer_CAEN)   {miny=-30; maxy=90;}
            show_trace(cDLED,trace_DLED[0], trace_DLED[1], trace_DLED_lenght, miny, maxy);
            show_line(cDLED, trace[0], mintp, maxtp, miny, maxy);
        }
        
        
        
        if(n_ev==last_event_n-1)
            reading=false;
        
        delete []trace[0];
        delete []trace[1];
        delete []peak;
        n_ev++;
    }
    if(average){
        TCanvas *cAVG = new TCanvas("AVG","AVG");
        cAVG->SetGrid();
        if(Agilent_MSO6054A){miny=50; maxy=90;}
        if(Digitizer_CAEN)  {miny=-855; maxy=-825;}
        show_trace(cAVG,trace_AVG[0], trace_AVG[1], trace_lenght, miny, maxy);
        show_line(cAVG, trace_AVG[0], mintp, maxtp, miny, maxy);
    }
    cHist->cd();
    if(SetLogyHist) cHist->SetLogy();
    ptrHist->Draw();
    
    
        
    return 0;
}


