//Ana_Traces_SiPM.C

//Read Ana_Traces_SiPM_ReadMe.txt before use

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
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TAxis.h"
#include "TPad.h"
#include "TLine.h"
#include "TSpectrum.h"
#include "TVirtualFitter.h"

#define nfilemax 3

//------------------------------------------------------------------------------
//--------------------------[   READ BIN DRS4 INTRO   ]-------------------------
//------------------------------------------------------------------------------
typedef struct {
   char           tag[3];
   char           version;
} FHEADER;

typedef struct {
   char           time_header[4];
} THEADER;

typedef struct {
   char           bn[2];
   unsigned short board_serial_number;
} BHEADER;

typedef struct {
   char           event_header[4];
   unsigned int   event_serial_number;
   unsigned short year;
   unsigned short month;
   unsigned short day;
   unsigned short hour;
   unsigned short minute;
   unsigned short second;
   unsigned short millisecond;
   unsigned short range;
} EHEADER;

typedef struct {
   char           tc[2];
   unsigned short trigger_cell;
} TCHEADER;

typedef struct {
   char           c[1];
   char           cn[3];
} CHEADER;

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------






//------------------------------------------------------------------------------
//-------------------------------[   FUNCTIONS   ]------------------------------
//------------------------------------------------------------------------------

//MAIN
int Analysis(string file, int last_event_n, bool display);

//PREDEFINED
int DCR_CT_1SiPM_1HVs(string file1, int last_event_n);
int DCR_CT_1SiPM_3HVs(string file1, string file2, string file3, int last_event_n);
int Ana1(string file1, int last_event_n, bool display_one_ev_param, bool LED_bool);

//SECONDARY
void DLED(int trace_lenght, int dleddt);
int find_peak_fix_time(int mintp, int maxtp);
int find_peaks(double thr_to_find_peaks, int max_peak_width, int min_peak_width,int blind_gap, bool DCR_DELAYS_bool);
void average_func(int trace_lenght);
void fit_hist_del(double expDelLow, double expDelHigh);
void fit_hist_all_peaks(TCanvas *c, TH1D *hist, double fit1Low, double fit1High, double fit2Low, double fit2High);
TGraphErrors* DCR_func(string file1, int last_event_n, int tot_files);
void Get_DCR_temp_and_errDCR_temp();
void show_trace(TCanvas* canv, double *x, double *y, int trace_lenght, double miny, double maxy, int mintp, int maxtp, bool line_bool, bool delete_bool, bool reverse);
void show_trace2(TCanvas* canv, double *x1, double *y1, double *x2, double *y2, int trace_lenght1, int trace_lenght2, double miny1, double maxy1, double miny2, double maxy2, int mintp, int maxtp, bool line_bool, bool delete_bool, bool reverse);

//READ FILE
void Read_Agilent_CAEN(string file, int last_event_n, bool display);
int ReadBin(string filename, bool display);
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
//--------------------------------[   DEVICE   ]--------------------------------
//------------------------------------------------------------------------------

bool Agilent_MSO6054A = false; //true if data taken by Agilent_MSO6054A, false otherwise
bool Digitizer_CAEN = true;  //true if data taken by Digitizer_CAEN, false otherwise
bool DRS4_Evaluation_Board = false; //true if data taken by DRS4_Evaluation_Board, false otherwise

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------



//-------------------------------------------------------------------------------
//-----------------------[   SETTING GLOBAL VARIABLES   ]------------------------
//-------------------------------------------------------------------------------
//peaks related
int mintp = 250; //min_time_peak
int maxtp = 280; //max_time_peak
int dleddt = 9; //10ns is approx the rise time used for HD3_2 on AS out 2
int blind_gap = 20; //ns
int max_peak_width = 20; //used for find_peaks
int min_peak_width =  5; //used for find_peaks
double pe_0_5 = 10; // !!! CHECK !!! (see comments below for values)
double pe_1_5 = 25; // !!! CHECK !!! (see comments below for values)
double thr_to_find_peaks = 10; //thr_to_find_peaks, as seen in DLED trace (in V); it should be similar to pe_0_5

//hist related
double maxyhist = 200;
double maxyhistAllPeaks = 200; 
double maxyhistDCR = 200;
double maxyhistDelays = 200;
double w = 1000;
double h = 800;
int bins_Volt = 204;
int bins_DCR = 206;
int bins_Delays = 100;

//fit related
double expDelLow_max= 60.;
double expDelHigh_max = 160.;
double fit1Low = 0;
double fit1High = 0;
double fit2Low = 0;
double fit2High = 0;

//display related
int ev_to_display = 7;


/* VALUES:
 *  DARK from Agilent: mintp = 160 (or 200);  maxtp = 320; dleddt = 39;
 *  LED from Agilent: approx in the middle, mintp = 420; maxtp = 500;  dleddt = 9; maxyhist = .2;
 *  LED from Digitizer_CAEN: depends on the wave
 */

//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------



//------------------------------------------------------------------------------
//---------------------------[   GLOBAL VARIABLES   ]---------------------------
//------------------------------------------------------------------------------
double **trace;
double **trace_DLED;
double **trace_AVG;
double **DCR;
double **errDCR;
double **DCR_thr;
double *peak;

TGraphErrors *gDCR_1;
TGraphErrors *gDCR_2;
TGraphErrors *gDCR_3;

double DCR_temp[] = {0., 0., 0.}; //I consider 3 files
double errDCR_temp[] = {0., 0., 0.};
double DCR_pe_0_5_vect[] = {0., 0., 0.};
double DCR_pe_1_5_vect[] = {0., 0., 0.};
double errDCR_pe_0_5_vect[] = {0., 0., 0.};
double errDCR_pe_1_5_vect[] = {0., 0., 0.};

int trace_DLED_lenght; int ii=0; int i=0; int index_func = 0; int nfile = 0; int n_DCR = 0; int DCR_cnt = 0; int trace_lenght = 0; int n_ev, index_for_peak; int one_window=0;
int nfiletot = 1;

double miny=0; double maxy=0; double miny1=0; double maxy1=0; double miny2=0; double maxy2=0; double gain, errgain; double DCR_time = 0.; double DCR_from_cnt = 0.; double max_func;

bool first_time_main_called = true; bool reading = true; bool last_event_flag = false;

char temp[20];
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
//-------------------------------[   OPTIONS   ]--------------------------------
//------------------------------------------------------------------------------

bool find_peak_in_a_selected_window = false; //to find a peak in a selected window (e.g. for LED measures)
bool average = false; //calculate the average of traces (useful for LED measures)

bool DCR_DELAYS_bool = false; //DCR from delays
bool CROSS_TALK_bool = false; //DCR must be true

bool drawHistAllPeaks = false; // to draw hist of all peaks in traces
bool fitHistAllPeaks = false; // fit hist of all peaks -> for GAIN
bool drawHistAllPeaksAll = false; // to draw hist of all peaks in traces for the 3 files (superimpose)
bool show_hists_DCR_DELAYS  = false;
bool showHist_bool = false; 
bool SetLogyHist = false;
bool running_graph = false;// to see traces in an osc mode (display must be true)
bool display_one_ev = false;
bool line_bool = false;

bool all_events_same_window = false; //all the events (from 0 to last_event_n) are joined in a single event

bool DO_NOT_DELETE_HIST_LED = false; //If set true, run only ONE TIME Analysis!!!

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------




//------------------------------------------------------------------------------
//--------------------[   GLOBAL HISTs, CANVs and FUNCs   ]---------------------
//------------------------------------------------------------------------------

TH1D *ptrHist = new TH1D("hist","",bins_Volt,0,maxyhist);

TH1D *ptrHistAllPeaks[nfilemax];
TH1D *ptrHistDCRthr[nfilemax];
TH1D *ptrHistDelays[nfilemax];

TF1 *expDel = new TF1("expDel","[1]*TMath::Exp(-[0]*x)",expDelLow_max,expDelHigh_max);
TF1 *gausFit1 = new TF1("gausFit1","gaus",-100,100);
TF1 *gausFit2 = new TF1("gausFit2","gaus",-100,100);

//I want to create ONLY one time the Canvas below... this is not the smartest way but it shuold work...
TCanvas *c = new TCanvas("Trace","Trace",w,h);
TCanvas *cDCR = new TCanvas("hist_DCR","hist_DCR",w,h);
TCanvas *cAllPeaks = new TCanvas("AllPeaks","AllPeaks",w,h);

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------



//------------------------------------------------------------------------------
//----------------------------[   MAIN FUNCTION   ]-----------------------------
//------------------------------------------------------------------------------   
int Analysis(string file, int last_event_n, bool display){
    gROOT->Reset();
    if(Agilent_MSO6054A or Digitizer_CAEN) 
        Read_Agilent_CAEN(file, last_event_n, display);
    if(DRS4_Evaluation_Board)
        ReadBin(file, display);
    
    
    
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
    
//***** HIST ALL PEAKS
    if(drawHistAllPeaks){
        if(nfile == 0){
            TCanvas *cAllPeaks1 = new TCanvas("AllPeaks1","AllPeaks1",w,h);
            cAllPeaks1-> SetGrid();
            cAllPeaks1->cd();
            ptrHistAllPeaks[nfile]->GetXaxis()->SetTitle("mV");
            ptrHistAllPeaks[nfile]->GetYaxis()->SetTitle("Counts");
            ptrHistAllPeaks[nfile]->Draw("hist");
            
            if(fitHistAllPeaks){
                if(file == "20180221_HD3-2_1_DARK_34_AS_2_01.txt" or file == "20180221_HD3-2_2_DARK_34_AS_2_02.txt" or file == "20180221_HD3-2_3_DARK_34_AS_2_01.txt")
                {fit1Low = 12; fit1High = 26; fit2Low = 28; fit2High = 42;}else{
                if(file == "20180221_HD3-2_1_DARK_35_AS_2_01.txt" or file == "20180221_HD3-2_2_DARK_35_AS_2_02.txt" or file == "20180221_HD3-2_3_DARK_35_AS_2_01.txt")
                {fit1Low = 12; fit1High = 28; fit2Low = 32; fit2High = 46;}else{
                if(file == "20180221_HD3-2_1_DARK_36_AS_2_01.txt" or file == "20180221_HD3-2_2_DARK_36_AS_2_02.txt" or file == "20180221_HD3-2_3_DARK_36_AS_2_01.txt")
                {fit1Low = 12; fit1High = 28; fit2Low = 36; fit2High = 50;}else{
                    fit1Low = 12; fit1High = 28; fit2Low = 36; fit2High = 50;
                }}}
                fit_hist_all_peaks(cAllPeaks1, ptrHistAllPeaks[nfile], fit1Low, fit1High, fit2Low, fit2High);
            }
            
        }
        if(drawHistAllPeaksAll){
            cAllPeaks->cd();
            
            if(nfile==1) ptrHistAllPeaks[nfile]->SetLineColor(kGreen+1);
            if(nfile==2) ptrHistAllPeaks[nfile]->SetLineColor(kRed+1);
            ptrHistAllPeaks[nfile]->Draw("histsame");
            
            auto legend = new TLegend(0.7,0.7,0.9,0.9);
            for(int k=0; k<nfilemax; k++){
                if(nfile==0) legend->AddEntry(ptrHistDCRthr[nfile],"HV = 34.00 V","l");
                if(nfile==1) legend->AddEntry(ptrHistDCRthr[nfile],"HV = 35.00 V","l");
                if(nfile==2) legend->AddEntry(ptrHistDCRthr[nfile],"HV = 36.00 V","l");
            }
            legend->Draw();
        }
    }
    
//***** DCR from DELAYS
    if(DCR_DELAYS_bool){
        if(show_hists_DCR_DELAYS){
            new TCanvas();
            ptrHistDelays[nfile]->GetXaxis()->SetTitle("Time (ns)");
            ptrHistDelays[nfile]->GetYaxis()->SetTitle("");
            ptrHistDelays[nfile]->Draw("hist");
        }
        
        cout<<"333"<<endl;
        fit_hist_del(expDelLow_max, expDelHigh_max);
        cout<<"222"<<endl;
    }
    
    if(showHist_bool){
        TCanvas *cHist = new TCanvas("hist_GAIN","hist_GAIN",w,h);
        cHist->SetGrid();
        cHist->cd();
        if(SetLogyHist) cHist->SetLogy();
        ptrHist->Draw("hist");
    }
    
    
        
    if(!DO_NOT_DELETE_HIST_LED)delete ptrHist;
    
    first_time_main_called = false;
    
    return 0;
}


//------------------------------------------------------------------------------
//-------------------------[   PREDEFINED FUNCTIONS   ]-------------------------
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
int DCR_CT_1SiPM_1HVs(string file1, int last_event_n){
    //TRUE:
    DCR_DELAYS_bool = true; //DCR from delays
    CROSS_TALK_bool = true; //DCR must be true
    DO_NOT_DELETE_HIST_LED = true;
       
    
    nfile = 0;
    
    pe_0_5 = 10; //mV
    pe_1_5 = 26; //mV
    
    ptrHistDelays[0]    = new TH1D("histDelays","",bins_Delays,0,maxyhistDelays);
    ptrHistAllPeaks[0]  = new TH1D("histAllPeaks","",bins_DCR,0,maxyhistAllPeaks);
    ptrHistDCRthr[0]    = new TH1D("histDCRthr","",bins_DCR,0,maxyhistDCR);
    
    
    gDCR_1 = DCR_func(file1,last_event_n, 1);
    
    TMultiGraph *DCR_mg = new TMultiGraph("DCR_mg", ";THR (mV); DCR (Hz)");
    DCR_mg->Add(gDCR_1);
   
    TCanvas *cDCR_loop = new TCanvas("cDCR_loop", "cDCR_loop");
    
    cDCR_loop->SetGrid();
    cDCR_loop->SetLogy();
    DCR_mg->Draw("A3L");
    
    cout<<endl<<endl;
    cout<<"-------------------------"<<endl;
    cout<<"-------[ RESULTS ]-------"<<endl;
    cout<<"-------------------------"<<endl<<endl;

    double CrossTalk = 0;
    double errCrossTalk = 0;
    
    CrossTalk = DCR_pe_1_5_vect[0]/DCR_pe_0_5_vect[0];
    errCrossTalk= CrossTalk * TMath::Sqrt( (errDCR_pe_0_5_vect[0]/DCR_pe_0_5_vect[0])*(errDCR_pe_0_5_vect[0]/DCR_pe_0_5_vect[0]) + (errDCR_pe_1_5_vect[0]/DCR_pe_1_5_vect[0])*(errDCR_pe_1_5_vect[0]/DCR_pe_1_5_vect[0]) );
    
    cout<<"File analyzed: "<<file1<<endl;
    
    cout<<"   DCR at 0.5 pe = ("<<DCR_pe_0_5_vect[0]*TMath::Power(10,-6)<<" +- "<<errDCR_pe_0_5_vect[0]*TMath::Power(10,-6)<<") MHz"<<endl;
    cout<<"   DCR at 1.5 pe = ("<<DCR_pe_1_5_vect[0]*TMath::Power(10,-6)<<" +- "<<errDCR_pe_1_5_vect[0]*TMath::Power(10,-6)<<") MHz"<<endl;
    cout<<"   Cross Talk    = ("<<CrossTalk<<" +- "<<errCrossTalk<<")"<<endl;

    return 0;
   
    
/* ---------------------------------------------------------------------
 * ----------------------------[   HD3-2   ]---------------------------- 
 * ---------------------------------------------------------------------
            HV      0.5pe   1.5pe   FILE
    --------------------------------------------------------------------
    SiPM1   34 V    10 mV   24 mV   20180221_HD3-2_1_DARK_34_AS_2_01.txt
            35 V    10 mV   26 mV   20180221_HD3-2_1_DARK_35_AS_2_01.txt
            36 V    10 mV   26 mV   20180221_HD3-2_1_DARK_36_AS_2_01.txt
    --------------------------------------------------------------------
    SiPM2   34 V    10 mV   24 mV   20180221_HD3-2_2_DARK_34_AS_2_02.txt
            35 V    10 mV   26 mV   20180221_HD3-2_2_DARK_34_AS_2_02.txt
            36 V    10 mV   26 mV   20180221_HD3-2_2_DARK_34_AS_2_02.txt
    --------------------------------------------------------------------
    SiPM3   34 V    10 mV   24 mV   20180221_HD3-2_3_DARK_34_AS_2_01.txt
            35 V    10 mV   26 mV   20180221_HD3-2_3_DARK_35_AS_2_01.txt
            36 V    10 mV   26 mV   20180221_HD3-2_3_DARK_36_AS_2_01.txt
*/
/*
 * CARLO's pe_0_5 = 7;
 */    

}


//------------------------------------------------------------------------------
int DCR_CT_1SiPM_3HVs(string file1, string file2, string file3, int last_event_n){
    
    //TRUE:
    DCR_DELAYS_bool = true; //DCR from delays
    CROSS_TALK_bool = true; //DCR must be true
    DO_NOT_DELETE_HIST_LED = true;
    
    nfiletot = 3;
    
        
    for(int k=0; k<nfiletot; k++){
        ptrHistDelays[k]   = new TH1D("histDelays","",bins_Delays,0,maxyhistDelays);
        ptrHistAllPeaks[k] = new TH1D("histAllPeaks","",bins_DCR,0,maxyhistAllPeaks);
        ptrHistDCRthr[k]   = new TH1D("histDCRthr","",bins_DCR,0,maxyhistDCR);
    }
    
    //file1:
    pe_0_5 = 10; pe_1_5 = 24;
    nfile = 0;
    gDCR_1 = DCR_func(file1,last_event_n, 3);
    //file2:
    nfile = 1;
    pe_0_5 = 10;   pe_1_5 = 26;
    gDCR_2 = DCR_func(file2,last_event_n, 3);
    //file3:
    nfile = 2;
    pe_0_5 = 10;   pe_1_5 = 26;
    gDCR_3 = DCR_func(file3,last_event_n, 3);  
    
    
    TMultiGraph *DCR_mg = new TMultiGraph("DCR_mg", ";THR (mV); DCR (Hz)");
    DCR_mg->Add(gDCR_1);
    DCR_mg->Add(gDCR_2);
    DCR_mg->Add(gDCR_3);
    
    TCanvas *cDCR_loop = new TCanvas("cDCR_loop", "cDCR_loop");
    
    cDCR_loop->SetGrid();
    cDCR_loop->SetLogy();
    DCR_mg->Draw("A3L");
    auto legendDCR_loop = new TLegend(0.75,0.75,0.9,0.9);
    legendDCR_loop->AddEntry(gDCR_1,"HV = 34.00 V","l");
    legendDCR_loop->AddEntry(gDCR_2,"HV = 35.00 V","l");
    legendDCR_loop->AddEntry(gDCR_3,"HV = 36.00 V","l");
    legendDCR_loop->Draw();
    
    
    cout<<endl<<endl;
    cout<<"-------------------------"<<endl;
    cout<<"-------[ RESULTS ]-------"<<endl;
    cout<<"-------------------------"<<endl<<endl;

    double CrossTalk[] = {0., 0., 0.};
    double errCrossTalk[] = {0., 0., 0.};
    for(int i=0; i<3; i++){
        CrossTalk[i] = DCR_pe_1_5_vect[i]/DCR_pe_0_5_vect[i];
        errCrossTalk[i]= CrossTalk[i] * TMath::Sqrt( (errDCR_pe_0_5_vect[i]/DCR_pe_0_5_vect[i])*(errDCR_pe_0_5_vect[i]/DCR_pe_0_5_vect[i]) + (errDCR_pe_1_5_vect[i]/DCR_pe_1_5_vect[i])*(errDCR_pe_1_5_vect[i]/DCR_pe_1_5_vect[i]) );
    }
    
    cout<<"Files analyzed: "<<file1<<", "<<file2<<", "<<file3<<endl;
    
    cout<<"double DCR[] =         {"<<DCR_pe_0_5_vect[0]*TMath::Power(10,-6)<<", "<<DCR_pe_0_5_vect[1]*TMath::Power(10,-6)<<", "<<DCR_pe_0_5_vect[2]*TMath::Power(10,-6)<<"};"<<endl;
    
    cout<<"double errDCR[] =      {"<<errDCR_pe_0_5_vect[0]*TMath::Power(10,-6)<<", "<<errDCR_pe_0_5_vect[1]*TMath::Power(10,-6)<<", "<<errDCR_pe_0_5_vect[2]*TMath::Power(10,-6)<<"};"<<endl;
    
    cout<<"double CrossTalk[] =   {"<<CrossTalk[0]<<", "<<CrossTalk[1]<<", "<<CrossTalk[2]<<"};"<<endl;

    cout<<"double errCrossTalk[] ={"<<errCrossTalk[0]<<", "<<errCrossTalk[1]<<", "<<errCrossTalk[2]<<"};"<<endl;
    
    return 0;
}

//------------------------------------------------------------------------------
int Ana1(string file1, int last_event_n, bool display_one_ev_param, bool LED_bool){
    //VARIABLES:
    //TRUE:
    drawHistAllPeaks = true; // to draw hist of all peaks in traces
    fitHistAllPeaks = true; // fit hist of all peaks -> for GAIN
    DCR_DELAYS_bool = true; //DCR from delays
    show_hists_DCR_DELAYS  = true;
    display_one_ev = display_one_ev_param;
    DO_NOT_DELETE_HIST_LED = true;
    
    pe_0_5 = 10;
    thr_to_find_peaks = pe_0_5;
    
    if(LED_bool){
        find_peak_in_a_selected_window = true;
        drawHistAllPeaks = false; // to draw hist of all peaks in traces
        fitHistAllPeaks = false; // fit hist of all peaks -> for GAIN
        DCR_DELAYS_bool = false; //DCR from delays
        show_hists_DCR_DELAYS  = false;
        showHist_bool = true;
        DO_NOT_DELETE_HIST_LED = true;
    }
    
    ptrHistAllPeaks[0]  = new TH1D("histAllPeaks","",bins_DCR,0,maxyhistAllPeaks);
    ptrHistDelays[0]    = new TH1D("histDelays","",bins_Delays,0,maxyhistDelays);
    ptrHistDCRthr[0]    = new TH1D("histDCRthr","",bins_DCR,0,maxyhistDCR);
    
    Analysis(file1, last_event_n, true);
    
    Get_DCR_temp_and_errDCR_temp();
    
    cout<<endl<<endl;
    cout<<"-------------------------"<<endl;
    cout<<"-------[ RESULTS ]-------"<<endl;
    cout<<"-------------------------"<<endl<<endl;
    
    cout<<"File analyzed: "<<file1<<endl;
    cout<<"   GAIN = ("<<gain<<" +- "<<errgain<<") mV"<<endl;
    cout<<"   DCR at 0.5 pe = ("<<DCR_pe_0_5_vect[0]*TMath::Power(10,-6)<<" +- "<<errDCR_pe_0_5_vect[0]*TMath::Power(10,-6)<<") MHz"<<endl;
    
    return 0;
}





//------------------------------------------------------------------------------
//-------------------------[   SECONDARY FUNCTIONS   ]--------------------------
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
    max_func = -10000;
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
int find_peaks(double thr_to_find_peaks, int max_peak_width, int min_peak_width,int blind_gap,  bool DCR_DELAYS_bool){ //I look for every peaks in the trace, only if I'm in DARK mode
    ii=2;
    int index_peak;
    int DCR_cnt_temp = 0;
    int index_old = 0;
    int index_new = 0;
    int peak_width = max_peak_width;
    while(ii<trace_DLED_lenght){//I find peaks after the DLED procedure
        if((trace_DLED[1][ii]>thr_to_find_peaks) and (trace_DLED[1][ii-2]<thr_to_find_peaks) and (ii+max_peak_width<trace_DLED_lenght)){//I only consider points above thr_to_find_peaks on the rising edge
            DCR_cnt_temp++; //I've seen a peak; if I'm in dark mode it's DCR
            
            //Now I want to see the peak amplitude.
            for(int k=ii+min_peak_width; (k<ii+max_peak_width); k++){//I open a window and look for a point below thr
                if(trace_DLED[1][k]<thr_to_find_peaks){
                    peak_width = k-ii;
                    break;
                }
            }
            
            //Now I look for the peak in that window
            if(ii+peak_width<trace_DLED_lenght)
                index_peak = find_peak_fix_time(ii, ii+peak_width);
            else 
                index_peak = find_peak_fix_time(ii, trace_DLED_lenght);
            
            if((index_peak-index_old)>blind_gap){
                index_new=index_peak;
                //DCR from the delay (Itzler Mark - Dark Count Rate Measure (pag 5 ss))
                if(DCR_DELAYS_bool){
                    //I fill the hist with delays between 1 pe peaks (or higher):
                    if(index_old>0){
                        ptrHistDelays[nfile] -> Fill(index_new - index_old); 
                    }
                }
                index_old = index_new;  
                ptrHistAllPeaks[nfile]->Fill(trace_DLED[1][index_peak]);
                
                ii=ii+peak_width;  
            }else{
                ii++;
            }
            }else{
                ii++;
        }
    }
    return DCR_cnt_temp;
}

//------------------------------------------------------------------------------
TGraphErrors *DCR_func(string file1, int last_event_n, int tot_files){ 
    
    bool display = false;
    int control = 0;
    
    thr_to_find_peaks = 7; //mV
    double max_thr_to_find_peaks = 41; //mV
    double gap_between_thr = 1; //mV
    
    n_DCR = (int)((max_thr_to_find_peaks - thr_to_find_peaks)/gap_between_thr);
    
    DCR = new double*[tot_files];
        for(int i = 0; i < 3; i++) {
            DCR[i] = new double[n_DCR];
    }
    errDCR = new double*[tot_files]; 
        for(int i = 0; i < 3; i++) {
            errDCR[i] = new double[n_DCR];
    }
    
    int h = 0;
    double *DCR_thr = new double[n_DCR]; 
    
    while(thr_to_find_peaks <= max_thr_to_find_peaks){
        control = Analysis(file1,last_event_n,display);
        cout<<"444"<<endl;
        Get_DCR_temp_and_errDCR_temp();
        cout<<"111"<<endl;
        ptrHistDelays[nfile]->Reset();
        cout<<"00"<<endl;
        DCR_thr[h] = thr_to_find_peaks; //mV
        DCR[nfile][h] = DCR_temp[nfile];
        errDCR[nfile][h] = errDCR_temp[nfile];
        thr_to_find_peaks = thr_to_find_peaks + gap_between_thr; //I jump to the next thr in order to evaluate the new DRC
        h++;
    }
    
    
    
    TGraphErrors *gDCR = new TGraphErrors(n_DCR, DCR_thr, DCR[nfile],NULL, errDCR[nfile]);
    
     gDCR->SetLineWidth(2);
    
    if(nfile==0){gDCR->SetLineColor(kGreen+1);gDCR->SetFillColorAlpha(kGreen+1, 0.3);}
    if(nfile==1){gDCR->SetLineColor(kBlue);   gDCR->SetFillColorAlpha(kBlue, 0.3);   }
    if(nfile==2){gDCR->SetLineColor(kRed+1);  gDCR->SetFillColorAlpha(kRed+1, 0.3);  }
    
    cout<<endl<<endl;
    
    return gDCR;
    
}

//------------------------------------------------------------------------------
void average_func(int trace_lenght){
    for (ii=0; ii<trace_lenght; ii++){
        trace_AVG[1][ii] = trace_AVG[1][ii] + trace[1][ii];
    }
}

//------------------------------------------------------------------------------
void show_trace(TCanvas* canv, double *x, double *y, int trace_lenght, double miny, double maxy, int mintp, int maxtp, bool line_bool, bool delete_bool, bool reverse, int section){
    if(section == 0) canv->cd();
    if(section == 1) canv->cd(1);
    if(section == 2) canv->cd(2);
    for(ii=0; ii<trace_lenght; ii++){
        x[ii] = x[ii]*TMath::Power(10,9);
        y[ii] = y[ii];
        
        if(reverse)
            y[ii] = -y[ii];
    }
    TGraphErrors *graph = new TGraphErrors(trace_lenght,x,y,0,0);
    graph->SetTitle("");
    graph->GetXaxis()->SetTitle("Time (ns)");
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

//------------------------------------------------------------------------------
void show_trace2(TCanvas* canv, double *x1, double *y1, double *x2, double *y2, int trace_lenght1, int trace_lenght2, double miny1, double maxy1, double miny2, double maxy2, int mintp, int maxtp, bool line_bool, bool delete_bool, bool reverse){
    
    for(ii=0; ii<trace_lenght1; ii++){
        x1[ii] = x1[ii]*TMath::Power(10,9);
        if(reverse)
            y1[ii] = -y1[ii];
    }
    for(ii=0; ii<trace_lenght2; ii++){
        x2[ii] = x2[ii]*TMath::Power(10,9);
    }
    canv->cd(1);
    TGraphErrors *graph1 = new TGraphErrors(trace_lenght1,x1,y1,0,0);
    graph1->SetTitle("");
    graph1->GetXaxis()->SetTitle("Time (ns)");
    graph1->GetYaxis()->SetTitle("Amplitude (mV)");
    graph1->GetYaxis()->SetTitleOffset(1.2);
    graph1->GetXaxis()->SetTitleOffset(1.2);
    graph1->GetYaxis()-> SetRangeUser(miny1,maxy1);
    graph1->Draw("aplsame");
    if(line_bool){
        TLine *lmin1 = new TLine (x1[mintp], miny1, x1[mintp], maxy1);
        TLine *lmax1 = new TLine (x1[maxtp], miny1, x1[maxtp], maxy1);
        lmin1->SetLineColor(kBlue);
        lmax1->SetLineColor(kBlue);
        lmin1->Draw("aplsame");
        lmax1->Draw("aplsame");
    }
    canv->Update();
    canv->cd(2);
    TGraphErrors *graph2 = new TGraphErrors(trace_lenght2,x2,y2,0,0);
    graph2->SetTitle("");
    graph2->GetXaxis()->SetTitle("Time (ns)");
    graph2->GetYaxis()->SetTitle("Amplitude (mV)");
    graph2->GetYaxis()->SetTitleOffset(1.2);
    graph2->GetXaxis()->SetTitleOffset(1.2);
    graph2->GetYaxis()-> SetRangeUser(miny2,maxy2);
    graph2->Draw("aplsame");
    if(line_bool){
        TLine *lmin2 = new TLine (x2[mintp], miny2, x2[mintp], maxy2);
        TLine *lmax2 = new TLine (x2[maxtp], miny2, x2[maxtp], maxy2);
        lmin2->SetLineColor(kBlue);
        lmax2->SetLineColor(kBlue);
        lmin2->Draw("aplsame");
        lmax2->Draw("aplsame");
    }
    
    canv->Update();
    if(delete_bool) {delete graph1; delete graph2;}
}

//------------------------------------------------------------------------------
void Get_DCR_temp_and_errDCR_temp(){
    DCR_temp[nfile] = expDel->GetParameter(0)*TMath::Power(10,9);
    errDCR_temp[nfile] = expDel->GetParError(0)*TMath::Power(10,9);
    
    if(CROSS_TALK_bool){
        if((int)(thr_to_find_peaks*10000)==(int)(pe_0_5*10000)){
            DCR_pe_0_5_vect[nfile] = DCR_temp[nfile];
            errDCR_pe_0_5_vect[nfile] = errDCR_temp[nfile];
        }
        else{
            if((int)(thr_to_find_peaks*10000)==(int)(pe_1_5*10000)){
                DCR_pe_1_5_vect[nfile] = DCR_temp[nfile];
                errDCR_pe_1_5_vect[nfile] = errDCR_temp[nfile];
            }
        }
    }
    
}

//------------------------------------------------------------------------------
void fit_hist_del(double expDelLow, double expDelHigh){ //fit hists filled with time delays in order to find DCR
    cout<<"Fit hist delays file "<<nfile<<endl;
    expDel-> SetParName(0,"DCR"); expDel-> SetParName(1,"MultCost");
    expDel->SetParameter(0,0.001);
    ptrHistDelays[nfile] -> Fit(expDel, "", "", expDelLow, expDelHigh);
    expDel->Draw("same");
}

//------------------------------------------------------------------------------
void fit_hist_all_peaks(TCanvas *c, TH1D *hist, double fit1Low, double fit1High, double fit2Low, double fit2High){
    
    double mean1, mean2;
    double errmean1, errmean2;
    
    c->cd();
    
    //1st
    hist -> Fit("gausFit1", "+", "", fit1Low, fit1High); //gaus fit of the 1st peak
    mean1 = gausFit1->GetParameter(1);
    errmean1 = gausFit1->GetParError(1);
    cout<<"mean1 = "<<mean1<<" +- "<<errmean1<<endl<<endl;
    gausFit1->Draw("same");
    
    //2nd
    hist -> Fit("gausFit2", "+", "", fit2Low, fit2High); //gaus fit of the 2nd peak
    mean2 = gausFit2->GetParameter(1);
    errmean2 = gausFit2->GetParError(1);
    cout<<"mean2 = "<<mean2<<" +- "<<errmean2<<endl<<endl;
    gausFit2->Draw("same");
    
    gain = mean2 - mean1; //gain is the difference between 1pe peak and 2pe peak
    errgain = TMath::Sqrt( errmean1*errmean1 + errmean2*errmean2 ); //error propagation
    
    cout<<endl;
    cout<<"GAIN = "<<gain<<" +- "<<errgain<<endl;
    
}


//------------------------------------------------------------------------------
//------------------------------[   READ FILES   ]------------------------------
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
void Read_Agilent_CAEN(string file, int last_event_n, bool display){
    ifstream OpenFile (file.c_str());
    reading = true;
    
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
                if(all_events_same_window){
                    trace_lenght = atoi(temp)*last_event_n;
                    one_window = atoi(temp);
                }
                OpenFile>>temp>>temp>>temp>>temp>>temp>>temp>>temp>>temp>>temp>>temp>>temp>>temp>>temp>>temp>>temp>>temp>>temp;
        }else{
            cout<<"ERROR: check acquisition device"<<endl;
        }
        }
        if(n_ev==0){
            cout<<"Acquisition device: ";
            if(Agilent_MSO6054A) cout<<"Agilent_MSO6054A"<<endl;
            if(Digitizer_CAEN)   cout<<"Digitizer_CAEN"<<endl;
            cout<<"Number points "<<trace_lenght<<endl;
        }        
        
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
                    trace[1][i]  = -atof(temp)/1024 * 1000; //1024 channels from 0 V to 1 V, expressed in mV
                    
                    if(i==(one_window-1)){
                        OpenFile>>temp>>temp; 
                        OpenFile>>temp;
                        OpenFile>>temp>>temp>>temp>>temp>>temp>>temp>>temp>>temp>>temp>>temp>>temp>>temp>>temp>>temp>>temp>>temp>>temp;
                    }
                    
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
        
        if(find_peak_in_a_selected_window){
            double *peak = new double[2];
            index_for_peak = find_peak_fix_time(mintp, maxtp);
            peak[0] = trace_DLED[0][index_for_peak];
            peak[1] = trace_DLED[1][index_for_peak];
            ptrHist->Fill(peak[1]);
        }
        
//***** DCR
        if(DCR_DELAYS_bool){
           find_peaks(thr_to_find_peaks,max_peak_width, min_peak_width,blind_gap,DCR_DELAYS_bool);
        }
        else{
            if(drawHistAllPeaks){//all peaks but not DCR
                thr_to_find_peaks = 10; //mV
                find_peaks(thr_to_find_peaks,max_peak_width, min_peak_width,blind_gap,DCR_DELAYS_bool);
            }
        }
        

//***** DISPLAY     
//void show_trace2(TCanvas* canv, double *x1, double *y1, double *x2, double *y2, int trace_lenght1, int trace_lenght2, double miny1, double maxy1, double miny2, double maxy2, int mintp, int maxtp, bool line_bool, bool delete_bool, bool reverse);

        if(display){
            if(!display_one_ev){
                if(n_ev==0){
                    c->Divide(1,2);
                    c->SetGrid();
                }
                if(Agilent_MSO6054A) {miny1=10; maxy1=180;}
                if(Digitizer_CAEN)   {miny1=700;  maxy1=900;}
                if(Agilent_MSO6054A) {miny2=-30; maxy2=90;}
                if(Digitizer_CAEN)   {miny2=-30; maxy2=90;}
                show_trace2(c, trace[0], trace[1], trace_DLED[0], trace_DLED[1], trace_lenght, trace_DLED_lenght, miny1, maxy1, miny2, maxy2, mintp, maxtp, line_bool, true, true);
                if(!running_graph)getchar();
            }else{
                if(n_ev==ev_to_display){
                    c->Divide(1,2);
                    c->SetGrid();
                    if(Agilent_MSO6054A) {miny1=10; maxy1=180;}
                    if(Digitizer_CAEN)   {miny1=700;  maxy1=900;}
                    if(Agilent_MSO6054A) {miny2=-30; maxy2=90;}
                    if(Digitizer_CAEN)   {miny2=-30; maxy2=90;}
                    show_trace2(c, trace[0], trace[1], trace_DLED[0], trace_DLED[1], trace_lenght, trace_DLED_lenght, miny1, maxy1, miny2, maxy2, mintp, maxtp, line_bool, false, true);
                }
            }
            
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
}

//------------------------------------------------------------------------------
int ReadBin(string filename, bool display)
{
   FHEADER  fh;
   THEADER  th;
   BHEADER  bh;
   EHEADER  eh;
   TCHEADER tch;
   CHEADER  ch;
   
   
   unsigned int scaler;
   unsigned short voltage[1024];
   double waveform[16][4][1024], time[16][4][1024];
   float bin_width[16][4][1024];
   int i, j, b, chn, n, chn_index, n_boards;
   double t1, t2, dt;

   double threshold;
   n_ev = 0;
   
//    if (argc > 1)
//       strcpy(filename, argv[1]);
//    else {
//       printf("Usage: read_binary <filename>\n");
//       return 0;
//    }
   
   int len = filename.length();
   char file_for_fopen[len];
   
   for(int k=0; k< len; k++){
       file_for_fopen[k] = filename[k];
    }   
   
   // open the binary waveform file
   FILE *f = fopen(file_for_fopen, "rb");
   if (f == NULL) {
      //printf("Cannot find file \'%s\'\n", filename);
      return 0;
   }

   // read file header
   fread(&fh, sizeof(fh), 1, f);
   if (fh.tag[0] != 'D' || fh.tag[1] != 'R' || fh.tag[2] != 'S') {
      //printf("Found invalid file header in file \'%s\', aborting.\n", filename);
      return 0;
   }
   
   if (fh.version != '2') {
      //printf("Found invalid file version \'%c\' in file \'%s\', should be \'2\', aborting.\n", fh.version, filename);
      return 0;
   }

   // read time header
   fread(&th, sizeof(th), 1, f);
   if (memcmp(th.time_header, "TIME", 4) != 0) {
      //printf("Invalid time header in file \'%s\', aborting.\n", filename);
      return 0;
   }

   for (b = 0 ; ; b++) {
      // read board header
      fread(&bh, sizeof(bh), 1, f);
      if (memcmp(bh.bn, "B#", 2) != 0) {
         // probably event header found
         fseek(f, -4, SEEK_CUR);
         break;
      }
      
      //printf("Found data for board #%d\n", bh.board_serial_number);
      
      // read time bin widths
      memset(bin_width[b], sizeof(bin_width[0]), 0);
      for (chn=0 ; chn<5 ; chn++) {
         fread(&ch, sizeof(ch), 1, f);
         if (ch.c[0] != 'C') {
            // event header found
            fseek(f, -4, SEEK_CUR);
            break;
         }
         i = ch.cn[2] - '0' - 1;
         //printf("Found timing calibration for channel #%d\n", i+1);
         fread(&bin_width[b][i][0], sizeof(float), 1024, f);
         // fix for 2048 bin mode: double channel
         if (bin_width[b][i][1023] > 10 || bin_width[b][i][1023] < 0.01) {
            for (j=0 ; j<512 ; j++)
               bin_width[b][i][j+512] = bin_width[b][i][j];
         }
      }
   }
   n_boards = b;
   
   
   // loop over all events in the data file
   for (n=0 ; ; n++) {
      // read event header
      i = (int)fread(&eh, sizeof(eh), 1, f);
      if (i < 1)
         break;
      
      //printf("Found event #%d %d %d\n", eh.event_serial_number, eh.second, eh.millisecond);
      
      // loop over all boards in data file
      for (b=0 ; b<n_boards ; b++) {
         
         // read board header
         fread(&bh, sizeof(bh), 1, f);
         if (memcmp(bh.bn, "B#", 2) != 0) {
            //printf("Invalid board header in file \'%s\', aborting.\n", filename);
            return 0;
         }
         
         // read trigger cell
         fread(&tch, sizeof(tch), 1, f);
         if (memcmp(tch.tc, "T#", 2) != 0) {
            //printf("Invalid trigger cell header in file \'%s\', aborting.\n", filename);
            return 0;
         }

    
         
         // reach channel data
         for (chn=0 ; chn<4 ; chn++) {
            // read channel header
            fread(&ch, sizeof(ch), 1, f);
            if (ch.c[0] != 'C') {
               // event header found
               fseek(f, -4, SEEK_CUR);
               break;
            }
            chn_index = ch.cn[2] - '0' - 1;
            fread(&scaler, sizeof(int), 1, f);
            fread(voltage, sizeof(short), 1024, f);
            
            for (i=0 ; i<1024 ; i++) {
               // convert data to volts
               waveform[b][chn_index][i] = (voltage[i] / 65536. + eh.range/1000.0 - 0.5);
               
               // calculate time for this cell
                // calculate time for this cell
               for (j=0,time[b][chn_index][i]=0 ; j<i ; j++)
                  time[b][chn_index][i] += bin_width[b][chn_index][(j+tch.trigger_cell) % 1024];
            }
         }
         
         // align cell #0 of all channels
         t1 = time[b][0][(1024-tch.trigger_cell) % 1024];
         for (chn=1 ; chn<4 ; chn++) {
            t2 = time[b][chn][(1024-tch.trigger_cell) % 1024];
            dt = t1 - t2;
            for (i=0 ; i<1024 ; i++)
               time[b][chn][i] += dt;
         }
         
        
//====================================================================================        
         
         trace_lenght = 1024;
         
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
    
        for(int k=0; k<trace_lenght; k++){
             trace[0][k] = time[0][0][k];
             trace[1][k] = waveform[0][0][k]*1000; //to convert in mV
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
        
        if(find_peak_in_a_selected_window){
            double *peak = new double[2];
            index_for_peak = find_peak_fix_time(mintp, maxtp);
            peak[0] = trace_DLED[0][index_for_peak];
            peak[1] = trace_DLED[1][index_for_peak];
            ptrHist->Fill(peak[1]);
        }
        
//***** DCR
        if(DCR_DELAYS_bool){
            find_peaks(thr_to_find_peaks,max_peak_width, min_peak_width,blind_gap,DCR_DELAYS_bool);
        }
        else{
            if(drawHistAllPeaks){//all peaks but not DCR
                thr_to_find_peaks = 10; //mV
                find_peaks(thr_to_find_peaks,max_peak_width, min_peak_width,blind_gap,DCR_DELAYS_bool);
            }
        }
        
        
        
//***** DISPLAY     
//void show_trace2(TCanvas* canv, double *x1, double *y1, double *x2, double *y2, int trace_lenght1, int trace_lenght2, double miny1, double maxy1, double miny2, double maxy2, int mintp, int maxtp, bool line_bool, bool delete_bool, bool reverse);

        if(display){
            if(!display_one_ev){
                if(n_ev==0){
                    c->Divide(1,2);
                    c->SetGrid();
                }
                miny1 = -100; maxy1 = 100; miny2 = -100; maxy2 = 100;
                show_trace2(c, trace[0], trace[1], trace_DLED[0], trace_DLED[1], trace_lenght, trace_DLED_lenght, miny1, maxy1, miny2, maxy2, mintp, maxtp, line_bool, true, true);
                if(!running_graph)getchar();
            }else{
                if(n_ev==ev_to_display){
                    c->Divide(1,2);
                    c->SetGrid();
                    miny1 = -100; maxy1 = 100; miny2 = -100; maxy2 = 100;
                    show_trace2(c, trace[0], trace[1], trace_DLED[0], trace_DLED[1], trace_lenght, trace_DLED_lenght, miny1, maxy1, miny2, maxy2, mintp, maxtp, line_bool, false, true);
                }
            }
            
        }
        
        
        delete []trace[0];
        delete []trace[1];
        delete []trace_DLED[0];
        delete []trace_DLED[1];
        delete []peak;
        n_ev++;
        
//====================================================================================        

      }//end loop boards
   }//loop events
   cout<<"Last event "<<n_ev<<endl;
   
   return 1;
}
