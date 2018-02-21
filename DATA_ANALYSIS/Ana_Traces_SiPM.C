//Ana_Traces_SiPM.C

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include "TMath.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TGraphErrors.h"
#include "TAxis.h"
#include "TPad.h"
#include "TLine.h"

//------------------------------------------------------------------------------
//-------------------------------[   FUNCTIONS   ]------------------------------
//------------------------------------------------------------------------------
void DLED(int trace_lenght, int dleddt);
int find_peak_fix_time(int mintp, int maxtp);
int find_peaks(double noise_level, int jump);
void average_func(int trace_lenght);
void show_trace(TCanvas* canv, double *x, double *y, int trace_lenght, double miny, double maxy, int mintp, int maxtp, bool line_bool, bool delete_bool, bool reverse);
void help();
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
//---------------------------[   GLOBAL VARIABLES   ]---------------------------
//------------------------------------------------------------------------------
double **trace;
double **trace_DLED;
double **trace_AVG;
double *peak;

int trace_DLED_lenght;

int ii=0;
double max_func;
int index_func = 0;


double w = 1000;
double h = 800;

int bins_Volt = 204;
int bins_DCR = 300;

double maxyhistAllPeaks = .2; 
double maxyhistDCR = 200;


TH1D *ptrHistAllPeaks = new TH1D("histAllPeaks","",bins_DCR,0,maxyhistAllPeaks);

/* VALUES:
 * 
 * LED from Digitizer_CAEN: bins_Volt = 204;
 * 
 */
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

void help(){
    cout<<"USAGE:"<<endl;
    cout<<"Analysis(string file, int last_event_n, bool display)\n"<<endl;
}


//------------------------------------------------------------------------------
//----------------------------[   MAIN FUNCTION   ]-----------------------------
//------------------------------------------------------------------------------   
int Analysis(string file, int last_event_n, bool display){
    gROOT->Reset();
    
//-------------------------------------------------------------------------------
//---------------------------[   SETTING VARIABLES   ]---------------------------
//-------------------------------------------------------------------------------
    int min_ind_offset = 0;
    int max_ind_offset = 80;
    int mintp = 250; //min_time_peak
    int maxtp = 280; //max_time_peak
    int dleddt = 9;
    double noise_level = 0.010; //noise_level, as seen in DLED trace (in V)
    int jump = 12; //used for find_peaks
    double maxyhist = .2;
   
    
    bool DCR_bool = true;
    bool average = false;
    
    bool Agilent_MSO6054A = false; //true if data taken by Agilent_MSO6054A, false otherwise
    bool Digitizer_CAEN = true;  //true if data taken by Digitizer_CAEN, false otherwise
    bool SetLogyHist = false;
    
    bool running_graph = false;
    
/* VALUES:
 * 
 * DARK from Agilent: mintp = 160 (or 200);  maxtp = 320; dleddt = 39;
 *                   
 * LED from Agilent: approx in the middle, mintp = 420; maxtp = 500;  dleddt = 9; maxyhist = .2;
 * 
 * LED from Digitizer_CAEN: depends on the wave
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
    
    int trace_lenght = 0;
    int DCR_cnt = 0;
    
    //TCanvas
    TCanvas *c = new TCanvas("Trace","Trace",w,h);
    c->SetGrid();
    TCanvas *cDLED = new TCanvas("DLED","DLED",w,h);
    cDLED->SetGrid();
    TCanvas *cHist = new TCanvas("hist_GAIN","hist_GAIN",w,h);
    cHist->SetGrid();
    TH1D *ptrHist = new TH1D("hist","",bins_Volt,0,maxyhist);
    
    TCanvas *cDCR = new TCanvas("hist_DCR","hist_DCR",w,h);
    
    double miny=0;
    double maxy=0;
    
    double DCR_time = 0.;
    double DCR = 0.;
    
    n_ev=0;

//***** READ FILE
    while(!OpenFile.eof() and (reading)){
        
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
                    trace[1][i]  = -atof(temp); //AdvanSid is an inverting amplifier, so I have negative signals. In order to analyze them, I reverse the signals. For the plot I can re-reverse them
        }
        }
        else{
            if(Digitizer_CAEN){
                for(i=0; i<trace_lenght; i++){
                    trace[0][i] = i*TMath::Power(10,-9); //1point=1ns
                    OpenFile>>temp;
                    trace[1][i]  = -atof(temp)/1024; //1024 channels from 0 V to 1 V
                }
        }
        }
        
//***** DLED
        DLED(trace_lenght,dleddt);
        //Now: trace_DLED
        
//***** AVERAGE
        if(average){
            if(n_ev==0){
                for(i=0; i< trace_lenght; i++){
                    trace_AVG[0][i]=trace[0][i];
                    trace_AVG[1][i]=0;
                }
            }
            average_func(trace_lenght);
        }
        
        index = find_peak_fix_time(mintp, maxtp);
        peak[0] = trace_DLED[0][index];
        peak[1] = trace_DLED[1][index];
        ptrHist->Fill(peak[1]);
        
//***** DCR
        if(DCR_bool){
            DCR_cnt = DCR_cnt + find_peaks(noise_level,jump);
            DCR_time = DCR_time + trace_lenght*TMath::Power(10,-9);
            //cout<<find_peaks(trace_lenght,noise_level,jump)<<endl;
        }
        
        
//***** DISPLAY        
        if(display){
            if(Agilent_MSO6054A) {miny=10; maxy=180;}
            if(Digitizer_CAEN)   {miny=700;  maxy=900;}
            show_trace(c,trace[0], trace[1], trace_lenght,miny,maxy,mintp,maxtp,true,true,true);
            if(Agilent_MSO6054A) {miny=-30; maxy=90;}
            if(Digitizer_CAEN)   {miny=-30; maxy=90;}
            show_trace(cDLED,trace_DLED[0], trace_DLED[1], trace_DLED_lenght, miny, maxy,mintp,maxtp,true,true,false);
            if(!running_graph)getchar();
        }
        
        
        
        if(n_ev==last_event_n-1)
            reading=false;
        
        delete []trace[0];
        delete []trace[1];
        delete []trace_DLED[0];
        delete []trace_DLED[1];
        delete []peak;
        n_ev++;
    }//file is closed
    
    cout<<"Last event "<<n_ev<<endl;
    
//***** AVERAGE
    if(average){
        average_func(trace_lenght);
        for(i=0; i<trace_lenght; i++){
            trace_AVG[1][i] = trace_AVG[1][i]/n_ev;
        }
        TCanvas *cAVG = new TCanvas("AVG","AVG");
        cAVG->SetGrid();
        if(Agilent_MSO6054A){miny=50; maxy=90;}
        if(Digitizer_CAEN)  {miny=700; maxy=855;}
        show_trace(cAVG,trace_AVG[0], trace_AVG[1], trace_lenght, miny, maxy,mintp,maxtp,true,false,true);
    }
    
//***** DCR    
    if(DCR_bool){
        DCR = DCR_cnt/DCR_time;
        cout<<"\nDCR = "<<DCR*TMath::Power(10,-6)<<" MHz"<<endl;
        cDCR->cd();
        
        TH1D *ptrHistDCRthr = new TH1D("histDCRthr","",bins_DCR,0,maxyhistDCR);
        double temp, binCenter;
        for(int i=0; i<bins_DCR; i++){
            temp=0;
            for(int j=i; j<bins_DCR; j++){
                temp = temp + ptrHistAllPeaks->GetBinContent(j);
            }
            binCenter = ptrHistAllPeaks->GetBinCenter(i);
            ptrHistDCRthr -> Fill(binCenter*1000,temp);
        }
        cDCR->cd();
        ptrHistDCRthr->GetXaxis()->SetTitle("Threshold (mV)");
        ptrHistDCRthr->Draw("hist");
        
    }
    
    cHist->cd();
    if(SetLogyHist) cHist->SetLogy();
    ptrHist->Draw();
        
    return 0;
}



//------------------------------------------------------------------------------
//---------------------------[   OTHER FUNCTIONS   ]----------------------------
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
void DLED(int trace_lenght, int dleddt){
    trace_DLED_lenght = trace_lenght - dleddt;
    //CREATE TRACE DLED
    trace_DLED = new double*[2];
    for(ii = 0; ii < 2; ii++) {
        trace_DLED[ii] = new double[trace_DLED_lenght];
    }
    for(ii=0; ii<trace_DLED_lenght; ii++){
        trace_DLED[0][ii] = trace[0][ii + dleddt];
        trace_DLED[1][ii] = trace[1][ii + dleddt]-trace[1][ii];
    }
}

//------------------------------------------------------------------------------
int find_peak_fix_time(int mintp, int maxtp){
    max_func = -100;
    index_func=0;
    for( ii=mintp; ii<maxtp; ii++){
        if(trace_DLED[1][ii]>max_func){
            max_func=trace_DLED[1][ii];
            index_func=ii;
        }
    }
    return index_func;
}

//------------------------------------------------------------------------------
int find_peaks(double noise_level, int jump){
    ii=0;
    int index_peak;
    int DCR_cnt_temp = 0;
    while(ii<trace_DLED_lenght){
        if(trace_DLED[1][ii]>noise_level){
            DCR_cnt_temp++;
            if(ii+jump<trace_DLED_lenght) index_peak = find_peak_fix_time(ii, ii+jump);
            else index_peak = find_peak_fix_time(ii, trace_DLED_lenght);
            ptrHistAllPeaks->Fill(trace_DLED[1][index_peak]);
            ii=ii+jump;
        }else{
            ii++;
        }
    }
    //cout<<DCR_cnt_temp<<endl;
    return DCR_cnt_temp;
}


//------------------------------------------------------------------------------
void average_func(int trace_lenght){
    for (ii=0; ii<trace_lenght; ii++){
        trace_AVG[1][ii] = trace_AVG[1][ii] + trace[1][ii];
    }
}

//------------------------------------------------------------------------------
void show_trace(TCanvas* canv, double *x, double *y, int trace_lenght, double miny, double maxy, int mintp, int maxtp, bool line_bool, bool delete_bool, bool reverse){
    canv->cd();
    for(ii=0; ii<trace_lenght; ii++){
        x[ii] = x[ii]*TMath::Power(10,9);
        y[ii] = y[ii]*TMath::Power(10,3);
        
        if(reverse)
            y[ii] = -y[ii];
    }
    TGraphErrors *graph = new TGraphErrors(trace_lenght,x,y,0,0);
    graph->SetTitle("");
    graph->GetXaxis()->SetTitle("Time (Ns)");
    graph->GetYaxis()->SetTitle("Amplitude (mV)");
    graph->GetYaxis()->SetTitleOffset(1.2);
    graph->GetXaxis()->SetTitleOffset(1.2);
    graph->GetYaxis()-> SetRangeUser(miny,maxy);
    graph->Draw("aplsame");
    if(line_bool){
        TLine *lmin = new TLine (x[mintp], miny, x[mintp], maxy);
        TLine *lmax = new TLine (x[maxtp], miny, x[maxtp], maxy);
        lmin->SetLineColor(kBlue);
        lmax->SetLineColor(kBlue);
        lmin->Draw("aplsame");
        lmax->Draw("aplsame");
    }
    canv->Update();
    if(delete_bool) {delete graph;}
}
