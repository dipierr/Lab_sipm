//Ana_Traces_SiPM.cxx

/********************************************************************************
 *  Ana_Traces_SiPM.cxx                                                         *
 *                                                                              *
 *  Read Ana_Traces_SiPM_ReadMe.md and/or the code before use.                  *
 *                                                                              *
 *  Key points:                                                                 *
 *  (1) Open root:                                                              *
 *          $ root -l                                                           *
 *  (2) Compile the macro:                                                      *
 *          root[0] .L Ana_Traces_SiPM.cxx++                                    *
 *  (3) Run the desired function:                                               *
 *          root [1] DCR_CT_1SiPM_1HV(string file1, int last_event_n);          *
 *          root [1] DCR_CT_1SiPM_3HVs(string file1, string file2, string file3,*
 *                    int last_event_n)                                         *
 *          root [1] Ana1(string file1, int last_event_n,                       *
 *                    bool display_one_ev_param);                               *
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
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TAxis.h"
#include "TPad.h"
#include "TLine.h"
#include "TSpectrum.h"
#include "TVirtualFitter.h"
#include "TFile.h"
#include "TTree.h"
#include "TSystem.h"

#define nfilemax 3
#define max_peak_num 50
#define max_peaks 5000000


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


//PREDEFINED
void DCR_CT_1SiPM_1HV(string file1, int last_event_n);
void DCR_CT_1SiPM_3HVs(string file1, string file2, string file3, int last_event_n);
void Ana1(string file1, int last_event_n, bool display_one_ev_param);
void Ana_LED(string file1, int last_event_n);
void Ana_Ped(string file1, int last_event_n);

//SECONDARY
void Analysis(string file, int last_event_n, bool display);
void DLED(int trace_length, int dleddt);
int find_peak_fix_time(int mintp, int maxtp);
void find_peaks(float thr_to_find_peaks, int max_peak_width, int min_peak_width,int blind_gap, bool DCR_DELAYS_bool);
void average_func(int trace_length);
void fit_hist_del(float expDelLow, float expDelHigh);
void fit_hist_all_peaks(TCanvas *c, TH1D *hist, float fit1Low, float fit1High, float fit2Low, float fit2High);
TGraphErrors* DCR_func(string file1, int last_event_n, int tot_files);
void Get_DCR_temp_and_errDCR_temp();
void show_trace(TCanvas* canv, float *x, float *y, int trace_length, float miny, float maxy, bool line_bool, bool delete_bool);
void show_trace2(TCanvas* canv, float *x1, float *y1, float *x2, float *y2, int trace_length1, int trace_length2, float miny1, float maxy1, float miny2, float maxy2, bool line_bool, bool delete_bool);
void find_peaks_from_vector();
void find_DCR_0_5_pe_and_1_5_pe();
void find_offset();
void find_offset_mod();
void find_offset_mod_2();
void find_offset_mod_3();
void find_offset_mod_4();
void subtract_offset();
void remove_peak_0_half();
void remove_peak_0_all();
void find_charge_selected_window(int mintp, int maxtp);

//READ FILE
void Read_Agilent_CAEN(string file, int last_event_n, bool display);
void ReadBin(string filename, int last_event_n, bool display);
void ReadRootFile(string filename, int last_event_n, bool display);


//HELP
void help();
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------



//-------------------------------------------------------------------------------
//-----------------------[   SETTING GLOBAL VARIABLES   ]------------------------
//-------------------------------------------------------------------------------

//-----------------
//---[ OPTIONS ]---
//-----------------

// DEVICE
bool Agilent_MSO6054A = false; //true if data taken by Agilent_MSO6054A, false otherwise
bool Digitizer_CAEN = false;  //true if data taken by Digitizer_CAEN, false otherwise
bool DRS4_Evaluation_Board = true; //true if data taken by DRS4_Evaluation_Board, false otherwise
bool DRS4_Evaluation_Board_Mod = false; //true if data taken by DRS4_Evaluation_Board mod version, false otherwise

// TRACK related options
bool reverse_bool = true; //true if the signal is negative
bool DLED_bool = true; //true to use the DLED technique

//-----------------
//-----------------


//-------------------------
//---[ TRACE SELECTION ]---
//-------------------------
int start_blind_gap = 20;
int end_bling_gap = 100;

//--------------
//---[ GSPS ]---
//--------------
double GSPS = 1;

//---------------
//---[ PEAKS ]---
//---------------

// DLED and PEAKS FINDING
int dleddt = 9*GSPS; //10ns is approx the rise time used for HD3_2 on AS out 2. Expressed in points: 9 @ 1GSPS
int blind_gap = 2*dleddt; //ns
int max_peak_width = 20; //used for find_peaks
int min_peak_width =  0; //used for find_peaks

// ONLY for DCR_CT_1SiPM_1HV and DCR_CT_1SiPM_3HVs:
float min_thr_to_find_peaks = 6;  //first thr value in the DCR vs thr plot (mV)
float max_thr_to_find_peaks = 50; //last thr value in the DCR vs thr plot (mV)
float gap_between_thr = 0.1; //gap between thresholds in the DCR vs thr plot (mV)
float min_pe_0_5 = 7;  //min value for 0.5pe threshold (mV)
float max_pe_0_5 = 15; //max value for 0.5pe threshold (mV)
float min_pe_1_5 = 28; //min value for 1.5pe threshold (mV)
float max_pe_1_5 = 33; //max value for 1.5pe threshold (mV)
int n_mean = 10; //number of points used for smoothing the DCR vs thr plot

// ONLY for Ana1:
float thr_to_find_peaks = 10; //thr_to_find_peaks, as seen in DLED trace (in V); it should be similar to pe_0_5. Only Ana1 does NOT change this values

// ONLY for LED measures
int minLED = 110; //charge window: min time for peak (ns)
int maxLED = 130; //charge window: max time for peak (ns)
int min_time_offset = 20; //min time for offset (ns)
int max_time_offset = 40; //max time for offset (ns)

//---------------
//---------------


//---------------
//---[ HISTS ]---
//---------------

float maxyhist = 200;
float maxyHistCharge = 50;
float maxyhistAllPeaks = 200; 
float maxyhistDCR = 200;
float maxyhistDelays = 200;
float w = 1000;
float h = 800;
int bins_Volt = 204;
int bins_DCR = 206;
int bins_Delays = 100;
int bins_Charge = 100;

//---------------
//---------------


//-------------
//---[ FIT ]---
//-------------

// DELAYS distribution
float expDelLow_max= 60.;
float expDelHigh_max = 160.;

//-------------
//-------------


//-----------------
//---[ DISPLAY ]---
//-----------------

int ev_to_display = 5;

//-----------------
//-----------------


/* VALUES
 *  HD3-2 dleddt = 9; min_peak_width = 5; max_peak_width = 20; maxyhistAllPeaks = 200; thr_to_find_peaks = 10;
 *  MPPC  dleddt = ;
 * 
 * 
 *  DARK from Agilent: mintp = 160 (or 200);  maxtp = 320; dleddt = 39;
 *  LED from Agilent: approx in the middle, mintp = 420; maxtp = 500;  dleddt = 9; maxyhist = .2;
 *  LED from Digitizer_CAEN: depends on the wave
 */

//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------



//------------------------------------------------------------------------------
//---------------------------[   GLOBAL VARIABLES   ]---------------------------
//------------------------------------------------------------------------------
float **trace;
float **trace_DLED;
float **trace_AVG;
float **DCR;
float **errDCR;
float **DCR_thr;
float *peak;
float *peaks;
float **thr_to_find_peaks_vect;
float **der_DCR;
float **thr_to_find_peaks_vect_mean;
float **DCR_mean;
float **errDCR_mean;

TGraphErrors *gDCR_1;
TGraphErrors *gDCR_2;
TGraphErrors *gDCR_3;

float pe_0_5_vect[3] = {1.,1.,1.};
float pe_1_5_vect[3] = {1.,1.,1.};

float DCR_temp[] = {0., 0., 0.}; //I consider 3 files
float errDCR_temp[] = {0., 0., 0.};
float DCR_pe_0_5_vect[] = {0., 0., 0.};
float DCR_pe_1_5_vect[] = {0., 0., 0.};
float errDCR_pe_0_5_vect[] = {0., 0., 0.};
float errDCR_pe_1_5_vect[] = {0., 0., 0.};

int trace_DLED_length; int ii=0; int i=0; int index_func = 0; int nfile = 0; int n_DCR = 0; int DCR_cnt = 0; int trace_length = 0; int n_ev, index_for_peak; int one_window=0;
int nfiletot = 1; int n_smooth = 0;

float miny=0; float maxy=0; float miny1=0; float maxy1=0; float miny2=0; float maxy2=0; float gain, errgain; float DCR_time = 0.; float DCR_from_cnt = 0.; float max_func;

bool first_time_main_called = true; bool reading = true; bool last_event_flag = false;
bool first_time_DCR_called = true;
char temp[20];

int num_peaks=0;
int index_vect[max_peak_num];
int mintp = 0; int maxtp = 0;

float peaks_all_delay[2][max_peaks];
int ind_peaks_all_delay = 0;
int n_ev_tot = 0;

float fit1Low = 0;
float fit1High = 0;
float fit2Low = 0;
float fit2High = 0;

float offset = 0.; float charge = 0.;

int max_bin, n_offset;

double min_peak_window, max_peak_window;
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
//-------------------------------[   OPTIONS   ]--------------------------------
//------------------------------------------------------------------------------
//set the following to false if you want to use DCR_CT_1SiPM_1HV, DCR_CT_1SiPM_3HVs or Ana1, since they will be modified by the functions.
//otherwise, do whatever you want 
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
bool display_peaks = false;
bool display_peaks_now = false;

bool all_events_same_window = false; //all the events (from 0 to last_event_n) are joined in a single event

bool DO_NOT_DELETE_HIST_LED = false; //If set true, run only ONE TIME Analysis!!!

bool find_peaks_bool = false;
bool find_offset_bool = false;
bool find_charge_window_bool = true;
bool remove_0_peak_bool = false;

bool fill_hist = false;

bool ptrAllTrace_bool = false;
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------




//------------------------------------------------------------------------------
//--------------------[   GLOBAL HISTs, CANVs and FUNCs   ]---------------------
//------------------------------------------------------------------------------

TH1D *ptrHist = new TH1D("hist","",140,0,200);

TH1F *ptrAll = new TH1F("histAll","",500,-100,100);
TH1F *ptrAllTrace = new TH1F("histAllT","",86,-10.0,20.0);

TH1F *offset_hist = new TH1F("offset_hist","Offset histogram;V [mV];freq",1000,-20.,20.);

TH1D *ptrHistCharge = new TH1D("HistCharge","",bins_Charge,-maxyHistCharge,maxyHistCharge);

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
//-------------------------[   PREDEFINED FUNCTIONS   ]-------------------------
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
void DCR_CT_1SiPM_1HV(string file1, int last_event_n){
    //TRUE:
    find_peaks_bool = true;
    DCR_DELAYS_bool = true; //DCR from delays
    CROSS_TALK_bool = true; //DCR must be true
    DO_NOT_DELETE_HIST_LED = true;
    
    nfile = 0; //I only consider 1 file
    
    ptrHistDelays[0]    = new TH1D("histDelays","",bins_Delays,0,maxyhistDelays);
    ptrHistAllPeaks[0]  = new TH1D("histAllPeaks","",bins_DCR,0,maxyhistAllPeaks);
    ptrHistDCRthr[0]    = new TH1D("histDCRthr","",bins_DCR,0,maxyhistDCR);
    
    // DCR_func
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

    float CrossTalk = 0;
    float errCrossTalk = 0;
    
    //CrossTalk = (DCR @ 1.5pe) / (DCR @ 0.5pe)
    CrossTalk = DCR_pe_1_5_vect[0]/DCR_pe_0_5_vect[0];
    errCrossTalk= CrossTalk * TMath::Sqrt( (errDCR_pe_0_5_vect[0]/DCR_pe_0_5_vect[0])*(errDCR_pe_0_5_vect[0]/DCR_pe_0_5_vect[0]) + (errDCR_pe_1_5_vect[0]/DCR_pe_1_5_vect[0])*(errDCR_pe_1_5_vect[0]/DCR_pe_1_5_vect[0]) );
    
    cout<<"File analyzed: "<<file1<<endl;
    cout<<"pe_0_5 = "<<pe_0_5_vect[0]<<" mV; pe_1_5 = "<<pe_1_5_vect[0]<<" mV"<<endl;
    cout<<"   DCR at 0.5 pe = ("<<DCR_pe_0_5_vect[0]*TMath::Power(10,-6)<<" +- "<<errDCR_pe_0_5_vect[0]*TMath::Power(10,-6)<<") MHz"<<endl;
    cout<<"   DCR at 1.5 pe = ("<<DCR_pe_1_5_vect[0]*TMath::Power(10,-6)<<" +- "<<errDCR_pe_1_5_vect[0]*TMath::Power(10,-6)<<") MHz"<<endl;
    cout<<"   Cross Talk    = ("<<CrossTalk<<" +- "<<errCrossTalk<<")"<<endl;

   
    
/* FEW THRs EXAMPLES:
 * 
 * ---------------------------------------------------------------------
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

 *           
 * CARLO's pe_0_5 = 7;
 */    

}


//------------------------------------------------------------------------------
void DCR_CT_1SiPM_3HVs(string file1, string file2, string file3, int last_event_n){
    
    //TRUE:
    find_peaks_bool = true;
    DCR_DELAYS_bool = true; //DCR from delays
    CROSS_TALK_bool = true; //DCR must be true
    DO_NOT_DELETE_HIST_LED = true;
    
    //I have 3 files
    nfiletot = 3;
  
        
    for(int k=0; k<nfiletot; k++){
        //In order to set 3 different titles:
        char h1[20], h2[20], h3[20];
        char k_temp[2];
        sprintf(h1, "histDelays");
        sprintf(h2, "histAllPeaks");
        sprintf(h3, "histDCRthr");
        sprintf(k_temp, "%d", k);
        
        strcat(h1,k_temp);
        strcat(h2,k_temp);
        strcat(h3,k_temp);
        
        //new hists:
        ptrHistDelays[k]   = new TH1D(strcat(h1,k_temp),"",bins_Delays,0,maxyhistDelays);
        ptrHistAllPeaks[k] = new TH1D(strcat(h2,k_temp),"",bins_DCR,0,maxyhistAllPeaks);
        ptrHistDCRthr[k]   = new TH1D(strcat(h3,k_temp),"",bins_DCR,0,maxyhistDCR);
    }
    
    //file1:
    nfile = 0;
    first_time_main_called = true; //will be set to false after the Analysis function is called
    ind_peaks_all_delay = 0;
    for(int i=0; i<max_peaks; i++){peaks_all_delay[0][i] = 0; peaks_all_delay[1][i] = 0;}
    // DCR_func
    gDCR_1 = DCR_func(file1,last_event_n, 3);
    
    
    //file2:
    nfile = 1;
    first_time_main_called = true; //will be set to false after the Analysis function is called
    ind_peaks_all_delay = 0;
    for(int i=0; i<max_peaks; i++){peaks_all_delay[0][i] = 0; peaks_all_delay[1][i] = 0;}
    // DCR_func
    gDCR_2 = DCR_func(file2,last_event_n, 3);
    
    
    //file3:
    nfile = 2;
    first_time_main_called = true; //will be set to false after the Analysis function is called
    ind_peaks_all_delay = 0;
    for(int i=0; i<max_peaks; i++){peaks_all_delay[0][i] = 0; peaks_all_delay[1][i] = 0;}
    // DCR_func
    gDCR_3 = DCR_func(file3,last_event_n, 3);  
    
    // graph DCR vs threshold evaluated at 3 different HVs
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

    float CrossTalk[] = {0., 0., 0.};
    float errCrossTalk[] = {0., 0., 0.};
    
    //CrossTalk = (DCR @ 1.5pe) / (DCR @ 0.5pe)
    for(int i=0; i<3; i++){
        CrossTalk[i] = DCR_pe_1_5_vect[i]/DCR_pe_0_5_vect[i];
        errCrossTalk[i]= CrossTalk[i] * TMath::Sqrt( (errDCR_pe_0_5_vect[i]/DCR_pe_0_5_vect[i])*(errDCR_pe_0_5_vect[i]/DCR_pe_0_5_vect[i]) + (errDCR_pe_1_5_vect[i]/DCR_pe_1_5_vect[i])*(errDCR_pe_1_5_vect[i]/DCR_pe_1_5_vect[i]) );
    }
    
    cout<<"Files analyzed:"<<endl;
    cout<<file1<<"   pe_0_5 = "<<pe_0_5_vect[0]<<" mV; pe_1_5 = "<<pe_1_5_vect[0]<<" mV"<<endl;
    cout<<file2<<"   pe_0_5 = "<<pe_0_5_vect[1]<<" mV; pe_1_5 = "<<pe_1_5_vect[1]<<" mV"<<endl;
    cout<<file3<<"   pe_0_5 = "<<pe_0_5_vect[2]<<" mV; pe_1_5 = "<<pe_1_5_vect[2]<<" mV"<<endl;
    cout<<"float DCR[] =         {"<<DCR_pe_0_5_vect[0]*TMath::Power(10,-6)<<", "<<DCR_pe_0_5_vect[1]*TMath::Power(10,-6)<<", "<<DCR_pe_0_5_vect[2]*TMath::Power(10,-6)<<"};"<<endl;
    
    cout<<"float errDCR[] =      {"<<errDCR_pe_0_5_vect[0]*TMath::Power(10,-6)<<", "<<errDCR_pe_0_5_vect[1]*TMath::Power(10,-6)<<", "<<errDCR_pe_0_5_vect[2]*TMath::Power(10,-6)<<"};"<<endl;
    
    cout<<"float CrossTalk[] =   {"<<CrossTalk[0]<<", "<<CrossTalk[1]<<", "<<CrossTalk[2]<<"};"<<endl;

    cout<<"float errCrossTalk[] ={"<<errCrossTalk[0]<<", "<<errCrossTalk[1]<<", "<<errCrossTalk[2]<<"};"<<endl;
    
}

//------------------------------------------------------------------------------
void Ana1(string file1, int last_event_n, bool display_one_ev_param){
    //VARIABLES:
    //TRUE:
    find_peaks_bool = true;
    drawHistAllPeaks = true; // to draw hist of all peaks in traces
    fitHistAllPeaks = true; // fit hist of all peaks -> for GAIN
    DCR_DELAYS_bool = true; //DCR from delays
    show_hists_DCR_DELAYS  = true;
    display_one_ev = display_one_ev_param;
    DO_NOT_DELETE_HIST_LED = true;
    display_peaks = true;

    //Charge:
    line_bool = true;
    find_charge_window_bool = true;
    minLED = 110; //ns
    maxLED = 130; //ns

        
    nfile = 0; //I only consider 1 file   
    
    ptrHistAllPeaks[0]  = new TH1D("histAllPeaks","",bins_DCR,0,maxyhistAllPeaks);
    ptrHistDelays[0]    = new TH1D("histDelays","",bins_Delays,0,maxyhistDelays);
    ptrHistDCRthr[0]    = new TH1D("histDCRthr","",bins_DCR,0,maxyhistDCR);
    
    //Analysis
    Analysis(file1, last_event_n, true);

new TCanvas();
    ptrHistCharge->Draw();
    ptrHist->Draw();
   
    
    //Get DCR (only @ threshold, set in the 'SETTING GLOBAL VARIABLES' section)
    Get_DCR_temp_and_errDCR_temp();
    DCR_pe_0_5_vect[nfile] = DCR_temp[nfile];
    errDCR_pe_0_5_vect[nfile] = errDCR_temp[nfile];
    
    
    cout<<endl<<endl;
    cout<<"-------------------------"<<endl;
    cout<<"-------[ RESULTS ]-------"<<endl;
    cout<<"-------------------------"<<endl<<endl;
    
    cout<<"File analyzed: "<<file1<<endl;
    cout<<"Fit range for GAIN: fit1Low = "<<fit1Low<<"; fit1High = "<<fit1High<<"; fit2Low = "<<fit2Low<<"; fit2High = "<<fit2High<<";"<<endl;
    cout<<"   GAIN = ("<<gain<<" +- "<<errgain<<") mV"<<endl;
    cout<<"   DCR at 0.5 pe = ("<<DCR_pe_0_5_vect[0]*TMath::Power(10,-6)<<" +- "<<errDCR_pe_0_5_vect[0]*TMath::Power(10,-6)<<") MHz"<<endl;
    
}

//------------------------------------------------------------------------------
void Ana_LED(string file1, int last_event_n){
    //VARIABLES:
    //TRUE:
    average = true;
    find_offset_bool = true;
    find_peak_in_a_selected_window = true;
    DO_NOT_DELETE_HIST_LED = true;
    
    find_charge_window_bool = true;
    
    min_peak_window = minLED;
    max_peak_window = maxLED;
    
    nfile = 0; //I only consider 1 file
    
    //Analysis
    Analysis(file1, last_event_n, false);
    
    new TCanvas();
    ptrHistCharge->Draw();
    ptrHist->Draw();
}

//------------------------------------------------------------------------------
void Ana_Ped(string file1, int last_event_n){
    //VARIABLES:
    //TRUE:
    find_offset_bool = true;
    
    //FALSE
    DLED_bool = false;
    fill_hist = false;
    remove_0_peak_bool = false;
    
    //Charge Pedestal:
    find_charge_window_bool = true;
    min_peak_window = 200; //ns
    max_peak_window = 250; //ns



    nfile = 0; //I only consider 1 file
    
    //Analysis
    Analysis(file1, last_event_n, false);
    
    new TCanvas();
    ptrAllTrace->Draw();
    
    TCanvas *c4 = new TCanvas();
    offset_hist->Draw();
    
    new TCanvas();
    ptrHistCharge->Draw();
    
}





//------------------------------------------------------------------------------
//-------------------------[   SECONDARY FUNCTIONS   ]--------------------------
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
void Analysis(string file, int last_event_n, bool display){
    gROOT->Reset();
    
    if(first_time_main_called){
        if(Agilent_MSO6054A or Digitizer_CAEN) 
            Read_Agilent_CAEN(file, last_event_n, display);
        if(DRS4_Evaluation_Board)
            ReadBin(file, last_event_n, display);
        if(DRS4_Evaluation_Board_Mod)
          ReadRootFile(file, last_event_n, display);
    }else{ //in order to speed up when I want to find peaks at different thresholds: I read the file only one time
        find_peaks_from_vector();
    }
    
    
    
//***** AVERAGE
    if(average){
        for(i=0; i<trace_length; i++){
            trace_AVG[1][i] = trace_AVG[1][i]/n_ev_tot;
        }
        
        TCanvas *cAVG = new TCanvas("AVG","AVG");
        cAVG->SetGrid();
        miny = -100; maxy = 30;
        show_trace(cAVG,trace_AVG[0], trace_AVG[1], trace_length, miny, maxy,true,false);
       
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
            cAllPeaks1->Update();
            
            if(fitHistAllPeaks){
                if(file == "20180221_HD3-2_1_DARK_34_AS_2_01.txt" or file == "20180221_HD3-2_2_DARK_34_AS_2_02.txt" or file == "20180221_HD3-2_3_DARK_34_AS_2_01.txt")
                {fit1Low = 12; fit1High = 26; fit2Low = 28; fit2High = 42;}else{
                if(file == "20180221_HD3-2_1_DARK_35_AS_2_01.txt" or file == "20180221_HD3-2_2_DARK_35_AS_2_02.txt" or file == "20180221_HD3-2_3_DARK_35_AS_2_01.txt")
                {fit1Low = 12; fit1High = 28; fit2Low = 32; fit2High = 46;}else{
                if(file == "20180221_HD3-2_1_DARK_36_AS_2_01.txt" or file == "20180221_HD3-2_2_DARK_36_AS_2_02.txt" or file == "20180221_HD3-2_3_DARK_36_AS_2_01.txt")
                {fit1Low = 12; fit1High = 28; fit2Low = 36; fit2High = 50;}else{
                    cout<<"--------------------"<<endl;
                    cout<<"File not recognized. Please enter:"<<endl;
                    cout<<"fit1Low  "; scanf("%f", &fit1Low);
                    cout<<"fit1High "; scanf("%f", &fit1High);
                    cout<<"fit2Low  "; scanf("%f", &fit2Low);
                    cout<<"fit2High "; scanf("%f", &fit2High);
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
        
        // fit_hist_del
        fit_hist_del(expDelLow_max, expDelHigh_max);
    }
    
    if(showHist_bool){
        TCanvas *cHist = new TCanvas("hist_GAIN","hist_GAIN",w,h);
        cHist->SetGrid();
        cHist->cd();
        if(SetLogyHist) cHist->SetLogy();
        ptrHist->Draw("hist");
    }
    
    
        
    if(!DO_NOT_DELETE_HIST_LED)delete ptrHist;
    
    //since this function is called, I set first_time_main_called to false
    first_time_main_called = false;
    
}

//------------------------------------------------------------------------------
TGraphErrors *DCR_func(string file1, int last_event_n, int tot_files){ 
    
    bool display = false;
    
    
    
    thr_to_find_peaks = min_thr_to_find_peaks;
    
    n_DCR = (int)((max_thr_to_find_peaks - min_thr_to_find_peaks)/gap_between_thr);
    
    if(first_time_DCR_called){
        DCR = new float*[tot_files];
            for(int i = 0; i < tot_files; i++) {
                DCR[i] = new float[n_DCR];
        }
        errDCR = new float*[tot_files]; 
            for(int i = 0; i < tot_files; i++) {
                errDCR[i] = new float[n_DCR];
        }
        thr_to_find_peaks_vect = new float*[tot_files]; 
            for(int i = 0; i < tot_files; i++) {
                thr_to_find_peaks_vect[i] = new float[n_DCR];
        }
        n_smooth = (int)n_DCR/n_mean;
        der_DCR = new float*[tot_files]; 
            for(int i = 0; i < tot_files; i++) {
                der_DCR[i] = new float[n_smooth];
        }
        thr_to_find_peaks_vect_mean = new float*[tot_files]; 
            for(int i = 0; i < tot_files; i++) {
                thr_to_find_peaks_vect_mean[i] = new float[n_smooth];
        }
        DCR_mean = new float*[tot_files]; 
            for(int i = 0; i < tot_files; i++) {
                DCR_mean[i] = new float[n_smooth];
        }
        errDCR_mean = new float*[tot_files]; 
            for(int i = 0; i < tot_files; i++) {
                errDCR_mean[i] = new float[n_smooth];
        }
    }
    first_time_DCR_called = false;
    
    int h = 0;
    float *DCR_thr = new float[n_DCR]; 
    
    while(thr_to_find_peaks <= max_thr_to_find_peaks){ // loop for DCR vs thresholds graph
        Analysis(file1,last_event_n,display); 
        Get_DCR_temp_and_errDCR_temp();
        ptrHistDelays[nfile]->Reset();
        DCR_thr[h] = thr_to_find_peaks; //mV
        DCR[nfile][h] = DCR_temp[nfile];
        errDCR[nfile][h] = errDCR_temp[nfile];
        thr_to_find_peaks_vect[nfile][h] = thr_to_find_peaks;
        thr_to_find_peaks = thr_to_find_peaks + gap_between_thr; //I jump to the next thr in order to evaluate the new DRC
        h++;
    }
    
    find_DCR_0_5_pe_and_1_5_pe();
    
    
    
    TGraphErrors *gDCR = new TGraphErrors(n_DCR, DCR_thr, DCR[nfile],NULL, errDCR[nfile]);
    
     gDCR->SetLineWidth(2);
    
    if(nfile==0){gDCR->SetLineColor(kGreen+1);gDCR->SetFillColorAlpha(kGreen+1, 0.3);}
    if(nfile==1){gDCR->SetLineColor(kBlue);   gDCR->SetFillColorAlpha(kBlue, 0.3);   }
    if(nfile==2){gDCR->SetLineColor(kRed+1);  gDCR->SetFillColorAlpha(kRed+1, 0.3);  }
    
    cout<<endl<<endl;
    
    return gDCR;
    
}
 
//------------------------------------------------------------------------------
void DLED(int trace_length, int dleddt){
    for(ii=0; ii<trace_DLED_length; ii++){
        trace_DLED[0][ii] = trace[0][ii + dleddt];
        trace_DLED[1][ii] = trace[1][ii + dleddt]-trace[1][ii];
    }
}

//------------------------------------------------------------------------------
int find_peak_fix_time(int mintp, int maxtp){
    max_func = -10000;
    index_func=0;
    for( ii=mintp; ii<maxtp; ii++){
        //cout<<mintp<<"\t"<<maxtp<<endl;
        if(trace_DLED[1][ii]>max_func){
            max_func=trace_DLED[1][ii];
            index_func=ii;
        }
    }
    return index_func;
}

//------------------------------------------------------------------------------
void find_peaks(float thr_to_find_peaks, int max_peak_width, int min_peak_width,int blind_gap,  bool DCR_DELAYS_bool){ //I look for every peaks in the trace, only if I'm in DARK mode
    ii=2*GSPS;
    int before_ind = (int)2*GSPS;
    int index_peak;
    int DCR_cnt_temp = 0;
    int index_old = 0;
    int index_new = 0;
    int peak_width = max_peak_width;
    num_peaks=0;
    for(int i=0; i<max_peak_num; i++)index_vect[i]=0;
    while(ii<trace_DLED_length){//I find peaks after the DLED procedure
        if((trace_DLED[1][ii]>thr_to_find_peaks) and (trace_DLED[1][ii-before_ind]<thr_to_find_peaks) and (ii+max_peak_width<trace_DLED_length)){//I only consider points above thr_to_find_peaks on the rising edge
            DCR_cnt_temp++; //I've seen a peak; if I'm in dark mode it's DCR
            
            //Now I want to see the peak amplitude.
            for(int k=ii+min_peak_width; (k<ii+max_peak_width); k++){//I open a window and look for a point below thr
                if(trace_DLED[1][k]<thr_to_find_peaks){
                    peak_width = k-ii;
                    break;
                }
            }
            
            //Now I look for the peak in that window
            if(ii+peak_width<trace_DLED_length)
                index_peak = find_peak_fix_time(ii, ii+peak_width);
            else 
                index_peak = find_peak_fix_time(ii, trace_DLED_length);
            
            if((index_peak-index_old)>blind_gap){
                index_new=index_peak;
                //DCR from the delay (Itzler Mark - Dark Count Rate Measure (pag 5 ss))
                if(DCR_DELAYS_bool){
                    //I fill the hist with delays between 1 pe peaks (or higher):
                    if(index_old>0){
                        ptrHistDelays[nfile] -> Fill(index_new - index_old); 
                        peaks_all_delay[0][ind_peaks_all_delay] = index_peak;
                        peaks_all_delay[1][ind_peaks_all_delay] = trace_DLED[1][index_peak];
                        ind_peaks_all_delay++;
                    }
                }
                index_old = index_new;  
                ptrHistAllPeaks[nfile]->Fill(trace_DLED[1][index_peak]);
                if(display_peaks_now and num_peaks<max_peak_num){
                    index_vect[num_peaks] = index_peak;
                    num_peaks++;
                }
                ii=ii+peak_width;  
            }else{
                ii++;
            }
            }else{
                ii++;
        }
    }
    peaks_all_delay[0][ind_peaks_all_delay] = -1;
    peaks_all_delay[1][ind_peaks_all_delay] = -1;
    ind_peaks_all_delay++;
    
}

//------------------------------------------------------------------------------
void find_peaks_from_vector(){
    int n,i,k;
    int index_old, index_new;
    index_old = index_new = 0;
    for(int i=0;i<ind_peaks_all_delay; i++){
        if(peaks_all_delay[0][i]==-1){
            index_old = 0;
            index_new = 0;
        }else{
        if(peaks_all_delay[1][i] > thr_to_find_peaks){
            index_new = peaks_all_delay[0][i];
            if(index_old>0){
                ptrHistDelays[nfile] -> Fill(index_new - index_old);
            }
            index_old = index_new;
        }
        }
        }
//          cout<<peaks_all_delay[0][i]<<"\t"<<peaks_all_delay[1][i]<<endl;   
    
}

//------------------------------------------------------------------------------
void show_trace(TCanvas* canv, float *x, float *y, int trace_length, float miny, float maxy, bool line_bool, bool delete_bool){
    for(ii=0; ii<trace_length; ii++){
        x[ii] = x[ii];
        y[ii] = y[ii];
        
        if(reverse_bool)
            y[ii] = -y[ii];
    }
    TGraphErrors *graph = new TGraphErrors(trace_length,x,y,0,0);
    graph->SetTitle("");
    graph->GetXaxis()->SetTitle("Time (ns)");
    graph->GetYaxis()->SetTitle("Amplitude (mV)");
    graph->GetYaxis()->SetTitleOffset(1.2);
    graph->GetXaxis()->SetTitleOffset(1.2);
    graph->GetYaxis()-> SetRangeUser(miny,maxy);
    graph->Draw("aplsame");
    if(line_bool){
        TLine *lmin = new TLine (minLED, miny, minLED, maxy);
        TLine *lmax = new TLine (maxLED, miny, maxLED, maxy);
        lmin->SetLineColor(kBlue);
        lmax->SetLineColor(kBlue);
        lmin->Draw("aplsame");
        lmax->Draw("aplsame");
    }
    canv->Update();
    if(delete_bool) {delete graph;}

}

//------------------------------------------------------------------------------
void show_trace2(TCanvas* canv, float *x1, float *y1, float *x2, float *y2, int trace_length1, int trace_length2, float miny1, float maxy1, float miny2, float maxy2, bool line_bool, bool delete_bool){
    
    for(ii=0; ii<trace_length1; ii++){
        x1[ii] = x1[ii];
        if(reverse_bool)
            y1[ii] = -y1[ii];
    }
    for(ii=0; ii<trace_length2; ii++){
        x2[ii] = x2[ii];
    }
    canv->cd(1);
    TGraphErrors *graph1 = new TGraphErrors(trace_length1,x1,y1,0,0);
    graph1->SetName("graph1");
    graph1->SetTitle("graph1");
    graph1->Draw("apl");
    graph1->GetXaxis()->SetTitle("Time (ns)");
    graph1->GetYaxis()->SetTitle("Amplitude (mV)");
    graph1->GetYaxis()->SetTitleOffset(1.2);
    graph1->GetXaxis()->SetTitleOffset(1.2);
    graph1->GetYaxis()-> SetRangeUser(miny1,maxy1);
    graph1->Draw("apl");
    graph1->SetEditable(kFALSE);
   if(line_bool){
        TLine *lmin = new TLine (minLED, miny1, minLED, maxy1);
        TLine *lmax = new TLine (maxLED, miny1, maxLED, maxy1);
        lmin->SetLineColor(kBlue);
        lmax->SetLineColor(kBlue);
        lmin->Draw("plsame");
        lmax->Draw("plsame");
    }
 
    canv->cd(2);
    TGraphErrors *graph2 = new TGraphErrors(trace_length2,x2,y2,0,0);
    graph2->SetName("graph2");
    graph2->SetTitle("graph2");
    graph2->Draw("apl");
    graph2->GetXaxis()->SetTitle("Time (ns)");
    graph2->GetYaxis()->SetTitle("Amplitude (mV)");
    graph2->GetYaxis()->SetTitleOffset(1.2);
    graph2->GetXaxis()->SetTitleOffset(1.2);
    graph2->GetYaxis()-> SetRangeUser(miny2,maxy2);
    graph2->Draw("apl");
    graph2->SetEditable(kFALSE);
if(line_bool){
        TLine *lmin = new TLine (minLED, miny2, minLED, maxy2);
        TLine *lmax = new TLine (maxLED, miny2, maxLED, maxy2);
        lmin->SetLineColor(kBlue);
        lmax->SetLineColor(kBlue);
        lmin->Draw("plsame");
        lmax->Draw("plsame");
    }
    
    TGraphErrors *graphPeaks = nullptr;
    TGraphErrors *graphPeaks_DLED = nullptr;

    if(display_peaks_now){
        float x_peaks[num_peaks], y_peaks[num_peaks], x_peaks_DLED[num_peaks], y_peaks_DLED[num_peaks];
        for(int i=0; i<num_peaks; i++){
            x_peaks[i] = trace[0][index_vect[i]+dleddt];
            y_peaks[i] = trace[1][index_vect[i]+dleddt];
            x_peaks_DLED[i] = trace_DLED[0][index_vect[i]];
            y_peaks_DLED[i] = trace_DLED[1][index_vect[i]];
        }

        //graphPeaks
        c->cd(1);
        graphPeaks = new TGraphErrors(num_peaks,x_peaks,y_peaks,0,0);
        graphPeaks->SetName("graphPeaks");
        graphPeaks->SetTitle("graphPeaks");
        graphPeaks->SetMarkerStyle(20);
        graphPeaks->SetMarkerColor(kRed);
        graphPeaks->Draw("psame");
        graphPeaks->SetEditable(kFALSE);
        

        //graphPeaks_DLED
        c->cd(2);
        graphPeaks_DLED = new TGraphErrors(num_peaks,x_peaks_DLED, y_peaks_DLED,0,0);
        graphPeaks_DLED->SetName("graphPeaks_DLED");
        graphPeaks_DLED->SetTitle("graphPeaks_DLED");
        graphPeaks_DLED->SetMarkerStyle(20);
        graphPeaks_DLED->SetMarkerColor(kGreen+1);
        graphPeaks_DLED->Draw("psame");
        graphPeaks_DLED->SetEditable(kFALSE);
        
    }
    
    canv->Update();

    // let the user to interact with the canvas, like zooming the axis
    // show next trace by pressing any key

    if(!running_graph){
      while(!gSystem->ProcessEvents()) {
         if (c->WaitPrimitive()==0)
         {
            break;
         }
      }
    }

    if(delete_bool) {
        delete graph1; delete graph2; 
        if(display_peaks_now){delete graphPeaks; delete graphPeaks_DLED;}
    }
}

//------------------------------------------------------------------------------
void Get_DCR_temp_and_errDCR_temp(){
    DCR_temp[nfile] = expDel->GetParameter(0)*TMath::Power(10,9);
    errDCR_temp[nfile] = expDel->GetParError(0)*TMath::Power(10,9);
}


//------------------------------------------------------------------------------
void find_DCR_0_5_pe_and_1_5_pe(){
    float min_0_5=100000000;
    float min_1_5=100000000;
    int min_index_0_5 = 0;
    int min_index_1_5 = 0;
    int j,k;
    
    //initialization
    for(int i=0; i<n_smooth; i++){
        der_DCR[nfile][i] = 0;
        thr_to_find_peaks_vect_mean[nfile][i] = 0;
        DCR_mean[nfile][i] = 0;
        errDCR_mean[nfile][i] = 0;
    }
    der_DCR[nfile][0] = min_1_5;
    
    //mean of n_mean elements
    j=k=0;
    while(j<n_DCR){
        for(int i=j; i<j+n_mean; i++){
            //sum
            DCR_mean[nfile][k] = DCR_mean[nfile][k] + DCR[nfile][i];
            thr_to_find_peaks_vect_mean[nfile][k] = thr_to_find_peaks_vect_mean[nfile][k] + thr_to_find_peaks_vect[nfile][i];
            
            //error propagation: sigma^2 = sigma1^2 + sigma2^2 + ...
            errDCR_mean[nfile][k] = errDCR_mean[nfile][k] + errDCR[nfile][i]*errDCR[nfile][i];
        }
        
        //division for n_mean
        DCR_mean[nfile][k] = DCR_mean[nfile][k]/n_mean;
        thr_to_find_peaks_vect_mean[nfile][k] = thr_to_find_peaks_vect_mean[nfile][k]/n_mean;
        
        //error propagation: sigma = sqrt(sigma^2)
        errDCR_mean[nfile][k] = TMath::Sqrt(errDCR_mean[nfile][k]);
        
        j = j+n_mean;
        k = k+1;
        
    }
    
    //abs(1st derivative of DCR mean)
    for(int i=1; i<n_smooth; i++){
        der_DCR[nfile][i] = TMath::Abs(DCR_mean[nfile][i] - DCR_mean[nfile][i-1]);
    }
    
    //find min for 0.5 pe and 1.5 pe
    for(int i=0; i<n_smooth; i++){
        //0.5 pe
        if(thr_to_find_peaks_vect_mean[nfile][i]>min_pe_0_5 and thr_to_find_peaks_vect_mean[nfile][i]<max_pe_0_5){
            if(der_DCR[nfile][i]<min_0_5){
                min_0_5 = der_DCR[nfile][i];
                min_index_0_5 = i;
            }
        }
        
        //1.5 pe
        if(thr_to_find_peaks_vect_mean[nfile][i]>min_pe_1_5 and thr_to_find_peaks_vect_mean[nfile][i]<max_pe_1_5){
            if(der_DCR[nfile][i]<min_1_5){
                min_1_5 = der_DCR[nfile][i];
                min_index_1_5 = i;
            }
        }
    }
    
    //0.5pe
    pe_0_5_vect[nfile] = thr_to_find_peaks_vect_mean[nfile][min_index_0_5];
    DCR_pe_0_5_vect[nfile] = DCR_mean[nfile][min_index_0_5];
    errDCR_pe_0_5_vect[nfile] = errDCR_mean[nfile][min_index_0_5];
   
    //1.5pe
    pe_1_5_vect[nfile] = thr_to_find_peaks_vect_mean[nfile][min_index_1_5];
    DCR_pe_1_5_vect[nfile] = DCR_mean[nfile][min_index_1_5];
    errDCR_pe_1_5_vect[nfile] = errDCR_mean[nfile][min_index_1_5];
    
}

//------------------------------------------------------------------------------
void fit_hist_del(float expDelLow, float expDelHigh){ //fit hists filled with time delays in order to find DCR
    cout<<"Fit hist delays file "<<nfile<<endl;
    expDel-> SetParName(0,"DCR"); expDel-> SetParName(1,"MultCost");
    expDel->SetParameter(0,0.001);
    ptrHistDelays[nfile] -> Fit(expDel, "", "", expDelLow, expDelHigh);
    expDel->Draw("same");
}

//------------------------------------------------------------------------------
void fit_hist_all_peaks(TCanvas *c, TH1D *hist, float fit1Low, float fit1High, float fit2Low, float fit2High){
    
    float mean1, mean2;
    float errmean1, errmean2;
    
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
void find_offset(){
    offset = 0.;
    int last_n, first_n;
    last_n = (int)(trace_length * max_time_offset / trace[0][trace_length- 1]);
    first_n = (int)(trace_length * min_time_offset / trace[0][trace_length- 1]);
    for(int i=first_n; i<last_n; i++){
        offset = offset + trace[1][i];
    }
    offset = offset/(last_n-first_n);
}


//------------------------------------------------------------------------------
void find_offset_mod(){  // USING TSpectrum

  float bin_width = (trace[0][1023] - trace[0][0])/1024;
  float x_low = trace[0][0] - bin_width/2.;
  float x_up = trace[0][1023] + bin_width/2.;
  float n_bin = x_up - x_low;

  // trace hist
  TH1F *trace_histo = new TH1F("histTrace","",n_bin,x_low,x_up);

  // filling trace hist:
  for(ii=x_low; ii<x_up; ii++){
        trace_histo->Fill(trace[0][ii], trace[1][ii]);
  }

 

  TH1* b = new TH1D();
  TSpectrum *s = new TSpectrum();
  
  //find BACKGROUND:
  b = s->Background(trace_histo,100);
  
  //subtract BACKGROUND:
  trace_histo->Add(b,trace_histo,-1,1);

  // update trace
  for(ii=x_low; ii<x_up; ii++){
        trace[0][ii] = trace_histo->GetBinCenter(ii) - bin_width/2.;
        trace[1][ii] = trace_histo->GetBinContent(ii);
        //ptrAllTrace->Fill(trace[1][ii]);
  }

  delete trace_histo;
  delete s;
  delete b;

}

//------------------------------------------------------------------------------
void find_offset_mod_2(){
    
  offset = 0.;

  TH1F *amplitude_hist = new TH1F("amplitude_hist","Amplitude histogram;V [mV];freq",1000,-60.,60.);
  
  for (int i = start_blind_gap; i < trace_length-end_bling_gap; ++i)
  {
    amplitude_hist->Fill(trace[1][i]);
  }

  offset = amplitude_hist->GetBinCenter(amplitude_hist->GetMaximumBin());
  //offset_hist->Fill(offset);

//   if (offset < -4.)
//   {
//     cout << "Offset for this trace is: " << offset << " mV" << endl;
//   }
  //TCanvas *c3 = new TCanvas();
  //amplitude_hist->Draw();

  delete amplitude_hist;

}

//------------------------------------------------------------------------------
void find_offset_mod_3(){

  TH1F *amplitude_hist = new TH1F("amplitude_hist","Amplitude histogram;V [mV];freq",1000,-60.,60.);
  
  for (int i = start_blind_gap; i < trace_length-end_bling_gap; ++i)
  {
    amplitude_hist->Fill(trace[1][i]);
  }

  offset = 0.;
  max_bin = 0;
  
  n_offset = 3;
  
  for(int i=0; i<n_offset; i++){
      max_bin = amplitude_hist->GetMaximumBin();
      offset += amplitude_hist->GetBinCenter(max_bin);
      amplitude_hist->SetBinContent(max_bin, 0);
  }
  
  offset=offset/n_offset;
  
  offset_hist->Fill(offset);

//   if (offset < -4.)
//   {
//     cout << "Offset for this trace is: " << offset << " mV" << endl;
//   }

  delete amplitude_hist;

}

//------------------------------------------------------------------------------
void find_offset_mod_4(){

    offset = 0.;
    for(int i=start_blind_gap; i<trace_length-end_bling_gap; i++){
        offset = offset + trace[1][i];
    }
    offset = offset/(trace_length - start_blind_gap - end_bling_gap);

}



//------------------------------------------------------------------------------
void subtract_offset(){
    for(int i=0; i<trace_length; i++){
        trace[1][i] = trace[1][i] - offset;
    }
}


//------------------------------------------------------------------------------
void remove_peak_0_half(){
    bool flag_0;
    flag_0=false;
    for(int i=0; i<trace_length; i++){
        if(trace[1][i]==0 and flag_0 and i!=0 and i!=(trace_length-1)){
            trace[1][i]=( trace[1][i-1] + trace[1][i+1] )/2;
            flag_0=false;
        }else{
            if(trace[1][i]==0)
                flag_0=true;
        }
        
    }
}


//------------------------------------------------------------------------------
void remove_peak_0_all(){
    
    for(int i = 1; i < trace_length - 1; i++){
        if(trace[1][i] == 0){
            trace[1][i] = (trace[1][i-1] + trace[1][i+1])/2;
        }
    }

}


//------------------------------------------------------------------------------
void find_charge_selected_window(int mintp, int maxtp){
    charge = 0;
    for(int i=mintp; i<maxtp; i++){
        charge = charge + trace[1][i];
    }
    ptrHistCharge->Fill(charge);
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
            trace_length = atoi(temp);
            OpenFile>>temp>>temp;
        }
        else{
            if(Digitizer_CAEN and not Agilent_MSO6054A){
                OpenFile>>temp>>temp; 
                OpenFile>>temp;
                trace_length = atoi(temp);
                if(all_events_same_window){
                    trace_length = atoi(temp)*last_event_n;
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
            cout<<"Number points "<<trace_length<<endl;
        }        
        
//***** CREATE TRACE
        //trace
        trace = new float*[2];
        for(i = 0; i < 2; i++) {
            trace[i] = new float[trace_length];
        }
        //trace_AVG
        if(average==true and n_ev==0){
            trace_AVG = new float*[2];
            for(i = 0; i < 2; i++) {
                trace_AVG[i] = new float[trace_length];
            }
        }


//***** READ TRACE FROM FILE
        if(Agilent_MSO6054A){
                for(i=0; i<trace_length; i++){
                    OpenFile>>temp;
                    trace[0][i] = atof(temp);
                    OpenFile>>temp;
                    trace[1][i]  = -atof(temp); //AdvanSid is an inverting amplifier, so I have negative signals. In order to analyze them, I reverse the signals. For the plot I can re-reverse them
        }
        }
        else{
            if(Digitizer_CAEN){
                for(i=0; i<trace_length; i++){
                    trace[0][i] = i; //1point=1ns
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
        if(DLED_bool){
            DLED(trace_length,dleddt);
            //Now: trace_DLED
        }else{
            for(ii=0; ii<trace_DLED_length; ii++){
                trace_DLED[0][ii] = trace[0][ii];
                trace_DLED[1][ii] = trace[1][ii];
            }
        }
        
//***** AVERAGE
        if(average){
            if(n_ev==0){
                for(i=0; i< trace_length; i++){
                    trace_AVG[0][i]=trace[0][i];
                    trace_AVG[1][i]=0;
                }
            }
            average_func(trace_length);
        }
        
        if(display){
            if(!display_one_ev) display_peaks_now = display_peaks;
            if(display_one_ev and (n_ev==ev_to_display)) display_peaks_now = display_peaks;
        } 
        
        if(find_peak_in_a_selected_window){
            float *peak = new float[2];
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
        if(display){
            miny1 = -10; maxy1 = 100; miny2 = -10; maxy2 = 100;
            if(!display_one_ev){
                if(n_ev==0){
                    c->Divide(1,2);
                    c->SetGrid();
                }
                show_trace2(c, trace[0], trace[1], trace_DLED[0], trace_DLED[1], trace_length, trace_DLED_length, miny1, maxy1, miny2, maxy2, line_bool, true);
            }else{
                if(n_ev==ev_to_display){
                    c->Divide(1,2);
                    c->SetGrid();
                    show_trace2(c, trace[0], trace[1], trace_DLED[0], trace_DLED[1], trace_length, trace_DLED_length, miny1, maxy1, miny2, maxy2, line_bool, false);
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
    
    n_ev_tot = n_ev;
    cout<<"Last event "<<n_ev_tot<<endl;
}

//------------------------------------------------------------------------------
void ReadBin(string filename, int last_event_n, bool display)
{
// from read_binary.cpp, created by Stefan Ritt
// use it to read binary from DRS4 Evaluation Board
// please do not modify until the end of this section (until the first "===" row)
   FHEADER  fh;
   THEADER  th;
   BHEADER  bh;
   EHEADER  eh;
   TCHEADER tch;
   CHEADER  ch;
   
   
   unsigned int scaler;
   unsigned short voltage[1024];
   float waveform[16][4][1024], time[16][4][1024];
   float bin_width[16][4][1024];
   int i, j, b, chn, n, chn_index, n_boards;
   float t1, t2, dt;

   float threshold;
   n_ev = 0;
   
   
   int len = filename.length();
   char file_for_fopen[len];
   
   for(int k=0; k< len; k++){
       file_for_fopen[k] = filename[k];
    }   
   
   // open the binary waveform file
   FILE *f = fopen(file_for_fopen, "rb");
   if (f == NULL) {
      //printf("Cannot find file \'%s\'\n", filename);
//       return 0;
       cout<<"Cannot find file"<<endl;
   }

   // read file header
   fread(&fh, sizeof(fh), 1, f);
   if (fh.tag[0] != 'D' || fh.tag[1] != 'R' || fh.tag[2] != 'S') {
      //printf("Found invalid file header in file \'%s\', aborting.\n", filename);
//       return 0;
       cout<<"Found invalid file header in file"<<endl;
   }
   
   if (fh.version != '2') {
      //printf("Found invalid file version \'%c\' in file \'%s\', should be \'2\', aborting.\n", fh.version, filename);
//       return 0;
      cout<<"Found invalid file version"<<endl;
   }

   // read time header
   fread(&th, sizeof(th), 1, f);
   if (memcmp(th.time_header, "TIME", 4) != 0) {
      //printf("Invalid time header in file \'%s\', aborting.\n", filename);
//       return 0;
      cout<<"Invalid time header in file"<<endl;
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
         // fix for 2048 bin mode: float channel
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
      
      if(n_ev%1000==0)
            cout<<"Read ev\t"<<n_ev<<endl;
      if(n_ev==last_event_n)
          break;
      
      //printf("Found event #%d %d %d\n", eh.event_serial_number, eh.second, eh.millisecond);
      
      // loop over all boards in data file
      for (b=0 ; b<n_boards ; b++) {
         
         // read board header
         fread(&bh, sizeof(bh), 1, f);
         if (memcmp(bh.bn, "B#", 2) != 0) {
            //printf("Invalid board header in file \'%s\', aborting.\n", filename);
//             return 0;
            cout<<"Invalid board header in file"<<endl;
         }
         
         // read trigger cell
         fread(&tch, sizeof(tch), 1, f);
         if (memcmp(tch.tc, "T#", 2) != 0) {
            //printf("Invalid trigger cell header in file \'%s\', aborting.\n", filename);
            cout<<"Invalid trigger cell header in file"<<endl;
//             return 0;
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
         
        
//==============================================================================
// now you can modify
         
         trace_length = 1024;
         
//***** CREATE TRACE
        //trace
        trace = new float*[2];
        for(i = 0; i < 2; i++) {
            trace[i] = new float[trace_length];
        }
        //trace_AVG
        if(average==true and n_ev==0){
            trace_AVG = new float*[2];
            for(i = 0; i < 2; i++) {
                trace_AVG[i] = new float[trace_length];
            }
        }
        trace_DLED_length = trace_length - dleddt;
        //CREATE TRACE DLED
        trace_DLED = new float*[2];
        for(ii = 0; ii < 2; ii++) {
            trace_DLED[ii] = new float[trace_DLED_length];
        }
    
        for(int k=0; k<trace_length; k++){
             trace[0][k] = time[0][0][k];
             if(reverse_bool) trace[1][k] = -waveform[0][0][k]*1000; //to convert in mV
             else trace[1][k] = waveform[0][0][k]*1000;
        }
        
        if(display){
            if(!display_one_ev) display_peaks_now = display_peaks;
            if(display_one_ev and (n_ev==ev_to_display)) display_peaks_now = display_peaks;
        }
        
        // FIND OFFSET
  /* if(find_offset_bool){
      find_offset_mod_3();
      subtract_offset();
    } // end if offset
*/
    
    if(ptrAllTrace_bool){
    for(int i=start_blind_gap; i<trace_length-end_bling_gap; i++){
        ptrAllTrace->Fill(trace[1][i]);
    }   
    }
            
//***** DLED
        if(DLED_bool){
            DLED(trace_length,dleddt);
            //Now: trace_DLED
        }else{
            for(ii=0; ii<trace_DLED_length; ii++){
                trace_DLED[0][ii] = trace[0][ii];
                trace_DLED[1][ii] = trace[1][ii];
            }
        }
        
        
        if(find_peak_in_a_selected_window){
            float *peak = new float[2];
            mintp = (int)(1024*min_peak_window/trace[0][trace_length-1]);
            maxtp = (int)(1024*max_peak_window/trace[0][trace_length-1]);
            index_for_peak = find_peak_fix_time(mintp, maxtp);
            peak[0] = trace_DLED[0][index_for_peak];
            peak[1] = trace_DLED[1][index_for_peak];
            ptrHist->Fill(peak[1]);
        }
        
//***** PEAKS FINDING
        if(find_peaks_bool){
            find_peaks(thr_to_find_peaks,max_peak_width, min_peak_width,blind_gap,DCR_DELAYS_bool);
        }
        
//***** FIND OFFSET
        if(find_offset_bool){
            find_offset();
            subtract_offset();
        }
        
//***** FIND CHARGE
        if(find_charge_window_bool){
            mintp = (int)(1024*min_peak_window/trace[0][trace_length-1]);
            maxtp = (int)(1024*max_peak_window/trace[0][trace_length-1]);
            find_charge_selected_window(mintp, maxtp);
        }
        
//***** AVERAGE
        if(average){
            if(n_ev==0){
                for(i=0; i< trace_length; i++){
                    trace_AVG[0][i]=trace[0][i];
                    trace_AVG[1][i]=0;
                }
            }
            for (int i=0; i<trace_length; i++){
                trace_AVG[1][i] = trace_AVG[1][i] + trace[1][i];
            }
        }        

        
//***** DISPLAY     
        if(display){
            miny1 = -180; maxy1 = 100; miny2 = -70; maxy2 = 150;
            if(!display_one_ev){
                if(n_ev==0){
                    c->Divide(1,2);
                    c->SetGrid();
                }
                show_trace2(c, trace[0], trace[1], trace_DLED[0], trace_DLED[1], trace_length, trace_DLED_length, miny1, maxy1, miny2, maxy2, line_bool, true);
            }else{
                if(n_ev==ev_to_display){
                    c->Divide(1,2);
                    c->SetGrid();
                    show_trace2(c, trace[0], trace[1], trace_DLED[0], trace_DLED[1], trace_length, trace_DLED_length, miny1, maxy1, miny2, maxy2, line_bool, false);
                }
            }
            
        }
        
        
        delete []trace[0];
        delete []trace[1];
        delete []trace_DLED[0];
        delete []trace_DLED[1];
        delete []peak;
        n_ev++;
//================================================================================
// Please do not modify below, until the end of the function

      }//end loop boards
   }//loop events
  n_ev_tot = n_ev;
  cout<<"Last event "<<n_ev_tot<<endl;

}

//------------------------------------------------------------------------------
void ReadRootFile(string filename, int last_event_n, bool display){
        
  TFile *file_root = TFile::Open(filename.c_str());
  TTree *wave      = (TTree *)file_root->Get("Waveforms");

  trace_length = 1024;

  int bad_peak;
  int nch = 0;

 /* if(!wave->FindBranch("time_ch_1") && !wave->FindBranch("time_ch_2") && !wave->FindBranch("volt_ch_1") && !wave->FindBranch("volt_ch_2")){

    cout << "At least two channels are needed (one wit the signal and one open for spike finding and removal!" << endl;
    cout << "Also check that the branches are called time_ch_[chn], where chn is the channel number (1,2,3 or 4)." << endl;
    cout << "After you fix this, run again this macro. Now it will exit." << endl;

    return;

  }
  else
  {
    nch += 2; 
  }
*/

  Int_t nentries = wave->GetEntries();

  if (last_event_n > nentries)
  {
    cout << "Last event is beyond number of entries of the TTree. Setting last event to number of TTree entries." << endl;
    last_event_n = nentries;
  }

  // create traces. they will be filled by the traces inside the root file

  trace = new float*[2];
  for(i = 0; i < 2; i++) {
    trace[i] = new float[trace_length];
  }

  wave->SetBranchAddress("time_ch_2",trace[0]);
  wave->SetBranchAddress("volt_ch_2",trace[1]);
  
  

  // create average trace if needed

  if(average){
    trace_AVG = new float*[2];
    for(i = 0; i < 2; i++) {
      trace_AVG[i] = new float[trace_length];
    }
  } // end if average

  // create DLED trace

  trace_DLED_length = trace_length - dleddt;
  
  trace_DLED = new float*[2];
  for(ii = 0; ii < 2; ii++) {
    trace_DLED[ii] = new float[trace_DLED_length];
  }

  // loop over ttree entries. Each entry is a trace.
  
  for (int entry = 0; entry < last_event_n; ++entry){

    if(entry%1000==0) cout<<"Read event "<<entry<<endl;

    bad_peak = 0;

    wave->GetEntry(entry);
    
    /*for(int i=0; i<trace_length; i++){
      if((trace[1][i]<-0.09) && (trace[1][i]>-0.11)) cout<<"Trace 0.01"<<endl;
      printf("%.5lf\t%.5lf\n",trace[0][i],trace[1][i]);
      cout<<trace[0][i]<<"\t"<<trace[1][i]<<endl;
    } */          
    
    //REMOVE PEAK AT 0
    if(remove_0_peak_bool){
        remove_peak_0_half();
    }
    
    
    
    for (int k = 0; k < trace_length; ++k)
    {
        if (reverse_bool){
        trace[1][k] = -trace[1][k];
        }
        
        //if (k!=0){
          //if (abs(trace[1][k] - trace[1][k-1]) > 10.)
          //{
           // bad_peak += 1;
          //}
        //}


    } // end for loop over waveform

      //cout << "Number of bad peaks: " << bad_peak << endl;

      //if (bad_peak >= 50)
      //{
        //cout << "Removing entry " << entry << endl;
        //continue;
      //}

    if(display){
      if(!display_one_ev) display_peaks_now = display_peaks;
      if(display_one_ev and (entry==ev_to_display)) display_peaks_now = display_peaks;
    } // end if display

    // DLED

    if(DLED_bool){
      DLED(trace_length,dleddt);
            //Now: trace_DLED
    }
    else
    {
      for(ii=0; ii<trace_DLED_length; ii++){
        trace_DLED[0][ii] = trace[0][ii];
        trace_DLED[1][ii] = trace[1][ii];
      }
    } // end if DLED

    if(find_peak_in_a_selected_window){
      float *peak = new float[2];
      index_for_peak = find_peak_fix_time(mintp, maxtp);
      peak[0] = trace_DLED[0][index_for_peak];
      peak[1] = trace_DLED[1][index_for_peak];
      ptrHist->Fill(peak[1]);
    } // end if find peak in selected window
        
    // PEAKS FINDING
    if(find_peaks_bool){
      find_peaks(thr_to_find_peaks,max_peak_width, min_peak_width,blind_gap,DCR_DELAYS_bool);
    } // end if peak finding
        
    // FIND OFFSET
    if(find_offset_bool){
      find_offset_mod_3();
      subtract_offset();
    } // end if offset

    
    for(int i=start_blind_gap; i<trace_length-end_bling_gap; i++){
        ptrAllTrace->Fill(trace[1][i]);
    }

/*    bool is_shitty = false;
    if(fill_hist){
      for (int l = 20; l < trace_length - 100; ++l)
      {
        if(trace[1][l]<-10.){
          //cout << "Entry number with V<-10.0 mV: " << entry << endl;
          is_shitty = true;
          break;
        }
        if(trace[1][l]>10.){
          cout << "Entry number with V<-10.0 mV: " << entry << endl;
          //is_shitty = true;
          //break;
        }
      }

      if (!is_shitty){
        for (int l = 20; l < trace_length - 100; ++l){ 
          ptrAllTrace->Fill(trace[1][l]);
        }
      }
    } // end if fill_hist
*/

    //***** FIND CHARGE for LED
    if(find_charge_window_bool){
      mintp = (int)(1024*minLED/trace[0][trace_length-1]);
      maxtp = (int)(1024*maxLED/trace[0][trace_length-1]);
      find_charge_selected_window(mintp, maxtp);
    } // end if find charge
        
    //***** AVERAGE
    if(average){

      if(DLED_bool){

        if(entry==0){
            for(i=0; i< trace_DLED_length; i++){
                trace_AVG[0][i]=trace_DLED[0][i];
                trace_AVG[1][i]=0;
            }
        }
        for (ii=0; ii<trace_DLED_length; ii++){
            trace_AVG[1][ii] = trace_AVG[1][ii] + trace_DLED[1][ii];
        }
      }
      else{

        if(entry==0){
            for(i=0; i< trace_length; i++){
                trace_AVG[0][i]=trace[0][i];
                trace_AVG[1][i]=0;
            }
        }
        for (ii=0; ii<trace_length; ii++){
            trace_AVG[1][ii] = trace_AVG[1][ii] + trace[1][ii];
        }
      }
    } // end if average     

        
    //***** DISPLAY     
    if(display){
        miny1 = -100; maxy1 = 100; miny2 = -100; maxy2 = 100;
        if(!display_one_ev){
            if(entry==0){
                c->Divide(1,2);
                c->SetGrid();
            }
            show_trace2(c, trace[0], trace[1], trace_DLED[0], trace_DLED[1], trace_length, trace_DLED_length, miny1, maxy1, miny2, maxy2, line_bool, true);
        }else{
            if(entry==ev_to_display){
                c->Divide(1,2);
                c->SetGrid();
                show_trace2(c, trace[0], trace[1], trace_DLED[0], trace_DLED[1], trace_length, trace_DLED_length, miny1, maxy1, miny2, maxy2, line_bool, false);
            }
        }   
    } // end if display

  } // end entry loop

  cout << endl << "Last event " << last_event_n <<endl;

  n_ev_tot = last_event_n;
    
  delete []trace[0];
  delete []trace[1];
  delete []trace_DLED[0];
  delete []trace_DLED[1];
  delete []peak;

  file_root->Close();
  delete file_root;
}


//------------------------------------------------------------------------------
//---------------------------------[   HELP   ]---------------------------------
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
void help(){
    cout<<endl;
    cout<<"Ana_Traces_SiPM.cxx"<<endl;
    cout<<"PREDEFINED FUNCTIONS:"<<endl;
    
    cout<<"\tvoid DCR_CT_1SiPM_1HV(string file1, int last_event_n);"<<endl;
    cout<<"\tvoid DCR_CT_1SiPM_3HVs(string file1, string file2, string file3, int last_event_n);"<<endl;
    cout<<"\tvoid Ana1(string file1, int last_event_n, bool display_one_ev_param);"<<endl;
    cout<<"\tvoid Ana_LED(string file1, int last_event_n);"<<endl;
    cout<<"\tvoid Ana_Ped(string file1, int last_event_n);"<<endl;
    
    cout<<endl;
    cout<<"See Ana_Traces_SiPM_ReadMe.md for more information"<<endl<<endl;
    cout<<"Davide Depaoli 2018"<<endl;
    cout<<endl;
}
