// Ana_Traces_SiPM.cxx

/******************************************************************************
 *  Ana_Traces_SiPM.cxx                                                       *
 *                                                                            *
 *                                                                            *
 *  Davide Depaoli 2018                                                       *
 *  Alessio Berti  2018                                                       *
 *                                                                            *
 ******************************************************************************/


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
#include <algorithm>
#include <functional>
#include <numeric>
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
#include "TTask.h"
#include "TVirtualGraphPainter.h"
#include "TGaxis.h"
#include "TROOT.h"
#include "TString.h"

using namespace std;



#define nfilemax 10
#define max_peak_num 50
#define max_peaks 5000000

#define tau -225

#define n6 TMath::Power(10,-6)


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


// PREDEFINED
void Ana1(string file1, int last_event_n, float threshold, bool display_one_ev_param);
void Ana3(string file1, string file2, string file3, int last_event_n);
void Ana_LED(string file1, int last_event_n);
void Ana_Ped(string file1, int last_event_n);
void DCR_CT_1SiPM_1HV(string file, int last_event_n);
void DCR_CT_1SiPM_nHVs(string filelist, int nfile_in_list, int last_event_n);


// SECONDARY
void Analysis(string file, int last_event_n, bool display, TCanvas *c);
void DLED(int trace_length, int dleddt);
int find_peak_fix_time(int mintp, int maxtp);
void find_peaks(float thr_to_find_peaks, int max_peak_width, int min_peak_width,int blind_gap, bool DCR_DELAYS_bool);
void FindPeaksRisingFalling(double thr, float **t, double length, int max_peak_width, int rising_points, int falling_points);
void FindPeakPositions(float* vector, Bool_t dled_bool, Int_t dt);
void FindPeaksFromPositions();
void fit_hist_peaks_0pe_1pe_2pe(TCanvas *canv, TH1D *hist);
void fit_hist_peaks_gaus_sum_012(TCanvas *canv, TH1D *hist, bool evaluate_cross_talk);
void fit_hist_del(float expDelLow, float expDelHigh);
void fit_hist_peaks(TCanvas *canv, TH1D *hist);
TGraphErrors* DCR_func(string file1, int last_event_n, int tot_files, TCanvas *c);
TGraphErrors* DCR_func_NO_Delays(string file1, int last_event_n, int tot_files, TCanvas *c);
void discriminator(double thr, float **t, double length);
void Get_DCR_temp_and_errDCR_temp();
void show_trace(TCanvas* canv, float *x, float *y, int trace_length, float miny, float maxy, bool line_bool, bool delete_bool);
void show_trace2(TCanvas* canv, float *x1, float *y1, float *x2, float *y2, int trace_length1, int trace_length2, float miny1, float maxy1, float miny2, float maxy2, bool line_bool, bool delete_bool);
void FindDelaysFromVector();
void FindDCRfromVector();
void find_DCR_0_5_pe_and_1_5_pe_auto();
void find_DCR_0_5_pe_and_1_5_pe_manual();
void DivideAreaDCR_DCR_0_5_DCR_1_5();
void find_offset();
void find_offset_mod_2();
void find_offset_mod_3();
void find_offset_mod_4();
void subtract_offset();
void remove_peak_0_half();
void remove_peak_0_all();
void find_charge_selected_window(int mintp, int maxtp);
double find_area_trace(double center, double low, double high);
void show_AVG_trace_window(TCanvas *c, float *tracet, float *tracev, int trace_length, bool delete_bool);
void DLED_offset_remove();
void smooth_trace_step();
void smooth_trace_3();
void smooth_trace_4();
void smooth_trace_5();
void SmoothTraceMarkov(int window);
void SmoothTraceN(int n);
Float_t GetMean(std::vector<Float_t> vec);
Float_t GetStdDev(std::vector<Float_t> vec);
void EvaluateCrossTalk(double DCR_0_5, double errDCR_0_5, double DCR_1_5, double errDCR_1_5, double* CT);
double GetDCRfromDelays();
double GetErrDCRfromDelays();




// READ FILE
void Read_Agilent_CAEN(string file, int last_event_n, bool display, TCanvas *c);
void ReadBin(string filename, int last_event_n, bool display, TCanvas *c);
void ReadRootFile(string filename, int last_event_n, bool display, TCanvas *canv);

// OLD
void DCR_CT_1SiPM_1HV_NO_Delays(string file1, int last_event_n);
void DCR_CT_1SiPM_3HVs_NO_Delays(string file1, string file2, string file3, int last_event_n);
void DCR_CT_1SiPM_5HVs_NO_Delays(string filelist, int last_event_n);
void DCR_discriminator(string file1, int last_event_n, int thr_to_find_peaks);
void DCR_CT_No_Stair(string file1, int last_event_n, float thr_0_5_pe, float thr_1_5_pe);
void DCR_CT_1SiPM_1HV_all_delays(string file1, int last_event_n);
void DCR_CT_1SiPM_3HVs_all_delays(string file1, string file2, string file3, int last_event_n);


// HELP
void help();
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------



//-----------------------------------------------------------------------------
//----------------------[   SETTING GLOBAL VARIABLES   ]-----------------------
//-----------------------------------------------------------------------------

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
bool DLED_offset_remove_bool = false;
bool fill_hist_peaks_when_found = true;

bool automatic_find_thr_1pe_2pe = false;

//-----------------
//-----------------


//---------------
//---[ TRACE ]---
//---------------
int start_blind_gap = 20;
int end_bling_gap = 100;
int n_smooth_trace = 2;

//--------------
//---[ GSPS ]---
//--------------
double GSPS = 1;

//---------------
//---[ PEAKS ]---
//---------------

// DLED and PEAKS FINDING
int dleddt = 5;//8;//5;//9*GSPS; //10ns is approx the rise time used for HD3_2 on AS out 2. Expressed in points: 9 @ 1GSPS
int blind_gap = 2*dleddt; //ns
int max_peak_width = 50; //used for find_peaks
int min_peak_width =  0; //used for find_peaks
int gap_between_peaks = 10;
int rise_time = dleddt;

// ONLY for DCR_CT_1SiPM_1HV and DCR_CT_1SiPM_3HVs:
float min_thr_to_find_peaks = 8;  //first thr value in the DCR vs thr plot (mV)
float max_thr_to_find_peaks = 80; //last thr value in the DCR vs thr plot (mV)
float gap_between_thr = 0.1; //gap between thresholds in the DCR vs thr plot (mV)
float min_pe_0_5 = 8;  //min value for 0.5pe threshold (mV)
float max_pe_0_5 = 15; //max value for 0.5pe threshold (mV)
float min_pe_1_5 = 28; //min value for 1.5pe threshold (mV)
float max_pe_1_5 = 33; //max value for 1.5pe threshold (mV)
int n_mean = 10; //number of points used for smoothing the DCR vs thr plot
float pe_0_5_vect[nfilemax] = {10., 10., 10., 10., 10., 10., 10., 10., 10., 10.};
float pe_1_5_vect[nfilemax] = {23., 25., 28., 30., 30., 30., 30., 30., 30., 30.};


// AREA
double Area = 36.00;

// ONLY for LED measures
int minLED_amp = 168;//290;//115;  // window: min time for peak (ns) for LED
int maxLED_amp = 176;//305;//125;  // window: max time for peak (ns) for LED
double time_area_low = 30;  // time for the area before the LED peak (ns)
double time_area_high = 200; // time for the area after the LED peak (ns)
int dcr_mintp  = minLED_amp + 200;
int dcr_maxtp  = maxLED_amp + 200;

// threshold
float thr_to_find_peaks = 8; //thr_to_find_peaks, as seen in DLED trace (in V); it should be similar to pe_0_5.


// 0 high, 1 low
float range1_low_low_mV   = 5;//5;  // 5;  //5;
float range1_low_high_mV  = 20;//25; // 25; //15;
// 1 high, 2 low
float range2_low_low_mV   = 20;//25; // 25; //15;
float range2_low_high_mV  = 30;//40; // 30; //25;
// 2 high
float range2_high_low_mV  = 30;//50; // 40; //30;
float range2_high_high_mV = 50;//60; // 50; //40;

int min_time_offset = 20; //min time for offset (ns)
int max_time_offset = 40; //max time for offset (ns)

int minLED_charge = 110;
int maxLED_charge = 180;

const int trace_window_length = maxLED_charge - minLED_charge + 1;

//---------------
//---------------


//---------------
//---[ HISTS ]---
//---------------

float maxyhist = 200;
float maxyHistCharge = 50;
float maxyhistAllPeaks = 200;
float maxyhistDCR = 200;
float maxyhistDelays = 500;
float w = 1000;
float h = 800;
int bins_Volt = 204;
int bins_DCR = 206;
int bins_Delays = 50;
int bins_Charge = 100;

//---------------
//---------------


//-------------
//---[ FIT ]---
//-------------

// DELAYS distribution
float expDelLow_max= 40.;
float expDelHigh_max = 500.;

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
float **AVG_trace_window;
float **DCR;
float **errDCR;
float **DCR_Area;
float **errDCR_Area;
float **DCR_thr;
float *peak;
float *peaks;
float *peak_LED;
float **thr_to_find_peaks_vect;
float **der_DCR;
float **thr_to_find_peaks_vect_mean;
float **DCR_mean;
float **errDCR_mean;

std::vector<Int_t> peak_pos;

TGraphErrors *gDCR[nfilemax];
TGraphErrors *gDCR_temp;

TGraphErrors *gDCR_1;
TGraphErrors *gDCR_2;
TGraphErrors *gDCR_3;
TGraphErrors *gDCR_4;
TGraphErrors *gDCR_5;


double DCR_temp[nfilemax];
double errDCR_temp[nfilemax];

double DCR_pe_0_5_vect[nfilemax];
double DCR_pe_1_5_vect[nfilemax];
double errDCR_pe_0_5_vect[nfilemax];
double errDCR_pe_1_5_vect[nfilemax];
double DCR_pe_0_5_delays_vect[nfilemax];
double DCR_pe_1_5_delays_vect[nfilemax];
double errDCR_pe_0_5_delays_vect[nfilemax];
double errDCR_pe_1_5_delays_vect[nfilemax];

double DCR_pe_0_5_Area_vect[nfilemax];
double DCR_pe_1_5_Area_vect[nfilemax];
double errDCR_pe_0_5_Area_vect[nfilemax];
double errDCR_pe_1_5_Area_vect[nfilemax];
double DCR_pe_0_5_Area_delays_vect[nfilemax];
double DCR_pe_1_5_Area_delays_vect[nfilemax];
double errDCR_pe_0_5_Area_delays_vect[nfilemax];
double errDCR_pe_1_5_Area_delays_vect[nfilemax];


int trace_DLED_length; int ii=0; int i=0; int index_func = 0; int nfile = 0; int n_DCR = 0; int DCR_cnt = 0; int DCR_cnt_temp = 0; int discriminator_cnt = 0; int trace_length = 0; int n_ev, index_for_peak; int one_window=0;
int nfiletot = 1; int n_smooth = 0;

float miny=0; float maxy=0; float miny1=0; float maxy1=0; float miny2=0; float maxy2=0; float gain, errgain; float DCR_time = 0.; float DCR_from_discriminator = 0.; float max_func;

double DCR_from_cnt = 0.;
double errDCR_from_cnt = 0.;

bool first_time_main_called = true; bool reading = true; bool last_event_flag = false;
bool first_time_DCR_called = true;
char temp[20];

int num_peaks=0;
int index_vect[max_peak_num] = {0};
int mintp = 0; int maxtp = 0;
int index_peak_LED = 0;

float peaks_all_delay[2][max_peaks];

int ind_peaks_all_delay = 0;
float peaks_all[max_peaks];
float peaks_all_time[max_peaks];
int ind_peaks_all = 0;
int n_ev_tot = 0;


float offset = 0.; float charge = 0.;

int max_bin, n_offset;

int min_line = 0;
int max_line = 0;

double min_peak_window[2], max_peak_window[2];

float fit0Low, fit0High, fit1Low, fit1High, fit2Low, fit2High;

int range1_low_low_bin = 0;
int range1_low_high_bin = 0;
int range2_low_low_bin = 0;
int range2_low_high_bin = 0;
int range2_high_low_bin = 0;
int range2_high_high_bin = 0;

double trace_time = 0;
double trace_time_raw = 0;

auto color_file_1 = kBlue+1;
auto color_file_2 = kBlue+1;
auto color_file_3 = kBlue+1;
auto color_file_4 = kBlue+1;
auto color_file_5 = kBlue+1;

int color_file[] = {kBlue,kBlue,kBlue,kBlue,kBlue,kBlue,kBlue,kBlue,kBlue,kBlue};

double opacity = 0.3;

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
bool fit_hist_del_bool = false;

bool drawHistAllPeaks = false; // to draw hist of all peaks in traces
bool fitHistAllPeaks = false; // fit hist of all peaks -> for GAIN
// bool drawHistAllPeaksAll = false; // to draw hist of all peaks in traces for the 3 files (superimpose)
bool show_hists_DCR_DELAYS  = false;
bool showHist_bool = false;
bool SetLogyHist = false;
bool running_graph = false; // to see traces in an osc mode (display must be true)
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
bool show_peak_LED_bool = false;
bool peak_rejected = false;
bool led_and_dcr_0pe = false;
bool find_area_trace_bool = false;
bool find_1phe_bool = false;
bool find_peaks_discriminator_bool = false;
bool DCR_from_cnt_bool = false;

bool smooth_trace_bool = false;


//------------------------------------------------------------------------------
//------------------------------------------------------------------------------




//------------------------------------------------------------------------------
//--------------------[   GLOBAL HISTs, CANVs and FUNCs   ]---------------------
//------------------------------------------------------------------------------

int histLED_low = -50;
int histLED_high = 700;
float histLED_binw = 0.3; //mV
int histLED_nbins = (int)(((float)histLED_high-(float)histLED_low)/histLED_binw);//750;
TH1D *ptrHistLED = new TH1D("histLED","",histLED_nbins, histLED_low, histLED_high);
TH1D *ptrHistArea = new TH1D("histArea","",700, 0, 7000);
TH1D *ptrHistDCR_window = new TH1D("HistDCR_window","",histLED_nbins, histLED_low, histLED_high);


TH1F *ptrAll = new TH1F("histAll","",500,-100,100);
TH1F *ptrAllTrace = new TH1F("histAllT","",86,-10.0,20.0);

TH1F *offset_hist = new TH1F("offset_hist","Offset histogram;V [mV];freq",1000,-20.,20.);

TH1D *ptrHistCharge = new TH1D("HistCharge","",bins_Charge,-maxyHistCharge,maxyHistCharge);

TH1D *ptrHistAllPeaks[nfilemax];
TH1D *ptrHistDCRthr[nfilemax];
TH1D *ptrHistDelays[nfilemax];

TF1 *expDel = new TF1("expDel","[1]*TMath::Exp(-[0]*x)",expDelLow_max,expDelHigh_max);
TF1 *gausFit0 = new TF1("gausFit0","gaus",-100,100);
TF1 *gausFit1 = new TF1("gausFit1","gaus",-100,100);
TF1 *gausFit2 = new TF1("gausFit2","gaus",-100,100);

TF1 *gaus_sum_012 = new TF1("gaus_sum_012", "[0]*TMath::Exp( - (x-[7])*(x-[7])/( 2*[5]*[5] ) ) + [1]*TMath::Exp( - (x-[7]-[4]-[8])*(x-[7]-[4]-[8])/( 2*([5]*[5] + [6]*[6] )) ) + [2]*TMath::Exp( - (x-[7]-2*[4])*(x-[7]-2*[4])/( 2*([5]*[5] + 4*[6]*[6] ) ) ) + [3]*TMath::Exp( - (x-[7]-3*[4])*(x-[7]-3*[4])/( 2*([5]*[5] + 9*[6]*[6]) ) )" ,-100,100);
//[H0]*TMath::Exp( - (x-[V0])*(x-[V0])/( 2*[s0]*[s0] ) ) + [H1]*TMath::Exp( - (x-[V0]-[gain])*(x-[V0]-[gain])/( 2*([s0]*[s0] + [sadd]*[sadd] )) ) + [H2]*TMath::Exp( - (x-[V0]-2*[gain])*(x-[V0]-2*[gain])/( 2*([s0]*[s0] + 4*[sadd]*[sadd] ) ) ) + [H3]*TMath::Exp( - (x-[V0]-3*[gain])*(x-[V0]-3*[gain])/( 2*([s0]*[s0] + 9*[sadd]*[sadd]) ) )



//I want to create ONLY one time the Canvas below... this is not the smartest way but it should work...
// TCanvas *c = new TCanvas("Trace","Trace",w,h);
// TCanvas *cDCR = new TCanvas("hist_DCR","hist_DCR",w,h);
// TCanvas *cAllPeaks = new TCanvas("AllPeaks","AllPeaks",w,h);

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------


// TCanvas *c1 = new TCanvas("c1","c1",w,h);


//------------------------------------------------------------------------------
//-------------------------[   PREDEFINED FUNCTIONS   ]-------------------------
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
void Ana1(string file1, int last_event_n, float thr, bool display_one_ev_param){
    //VARIABLES:
    //TRUE:
    find_peaks_bool = true;
    drawHistAllPeaks = true; // to draw hist of all peaks in traces
    fitHistAllPeaks = false; // fit hist of all peaks -> for GAIN
    DCR_DELAYS_bool = true; //DCR from delays
    fit_hist_del_bool = true;
    show_hists_DCR_DELAYS  = true;
    display_one_ev = display_one_ev_param;
    DO_NOT_DELETE_HIST_LED = true;
    display_peaks = true;
    find_offset_bool = true;
    DCR_from_cnt_bool = true;

    smooth_trace_bool = false;


    //Charge:
    line_bool = false;
    find_charge_window_bool = true;

    nfile = 0; //I only consider 1 file

    thr_to_find_peaks = thr;

    ptrHistAllPeaks[0]  = new TH1D("histAllPeaks","",bins_DCR,0,maxyhistAllPeaks);
    ptrHistDelays[0]    = new TH1D("histDelays","",bins_Delays,0,maxyhistDelays);
    ptrHistDCRthr[0]    = new TH1D("histDCRthr","",bins_DCR,0,maxyhistDCR);

    TCanvas *c = new TCanvas("Trace","Trace",w,h);

    //Analysis
    Analysis(file1, last_event_n, true, c);

    TCanvas *chargehist = new TCanvas("chargehist", "chargehist", 5);
    chargehist->cd();
    ptrHistCharge->Draw();
    ptrHistLED->Draw();


    //Get DCR (only @ threshold, set in the 'SETTING GLOBAL VARIABLES' section)
    Get_DCR_temp_and_errDCR_temp();
    DCR_pe_0_5_vect[nfile] = DCR_temp[nfile];
    errDCR_pe_0_5_vect[nfile] = errDCR_temp[nfile];


    cout<<endl<<endl;
    cout<<"-------------------------"<<endl;
    cout<<"-------[ RESULTS ]-------"<<endl;
    cout<<"-------------------------"<<endl<<endl;

    cout<<"File analyzed: "<<file1<<endl;
    // cout<<"Fit range for GAIN: fit1Low = "<<fit1Low<<"; fit1High = "<<fit1High<<"; fit2Low = "<<fit2Low<<"; fit2High = "<<fit2High<<";"<<endl;
    // cout<<"   GAIN = ("<<gain<<" +- "<<errgain<<") mV"<<endl;
    // cout<<"   DCR at 0.5 pe = "<<endl;
    cout<<"   DCR at "<<thr<<" mV = "<<endl;
    cout<<"         ("<<DCR_pe_0_5_vect[0]*TMath::Power(10,-6)<<" +- "<<errDCR_pe_0_5_vect[0]*TMath::Power(10,-6)<<") MHz, from exp fit"<<endl;
    cout<<"         ("<<DCR_from_cnt*TMath::Power(10,-6)<<" +- "<<errDCR_from_cnt*TMath::Power(10,-6)<<") MHz, from cnt"<<endl;

    //////////////
    /// TEST Z ///
    //////////////

    // double test_z_DCR;
    // test_z_DCR = TMath::Abs(DCR_pe_0_5_vect[0] - DCR_from_cnt) / TMath::Sqrt( errDCR_pe_0_5_vect[0]*errDCR_pe_0_5_vect[0] + errDCR_from_cnt*errDCR_from_cnt );
    //
    // if(test_z_DCR<=1.96){
    //     cout<<"DCR are COMPATIBLE, z = "<<test_z_DCR<<endl;
    // }else{
    //     cout<<"DCR are *NOT* COMPATIBLE, z = "<<test_z_DCR<<endl;
    // }

    cout<<"*****************************"<<endl;

}


//------------------------------------------------------------------------------
void Ana3(string file1, string file2, string file3, int last_event_n){

    //TRUE:
    find_peaks_bool = true;
    DCR_DELAYS_bool = true; //DCR from delays
    fit_hist_del_bool = true;
    CROSS_TALK_bool = true; //DCR must be true
    DO_NOT_DELETE_HIST_LED = true;

    smooth_trace_bool = true;


    //I have 3 files
    nfiletot = 3;

    TCanvas *c = new TCanvas("Trace","Trace",w,h);


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
    first_time_main_called=true;
    Analysis(file1, last_event_n, false, c);


    //file2:
    nfile = 1;
    first_time_main_called=true;
    Analysis(file2, last_event_n, false, c);

    //file3:
    nfile = 2;
    first_time_main_called=true;
    Analysis(file3, last_event_n, false, c);


    // Draw peak distribution for 3 different HVs
    TCanvas *cAllPeaks3 = new TCanvas("hist_All3","hist_All3",w,h);
    cAllPeaks3->cd();
    cAllPeaks3->SetGrid();
    ptrHistAllPeaks[1]->SetLineColor(kGreen+1);
    ptrHistAllPeaks[2]->SetLineColor(kRed+1);

    for(int i=0; i<nfiletot; i++){
        ptrHistAllPeaks[i]->GetXaxis()->SetTitle("mV");
        ptrHistAllPeaks[i]->GetYaxis()->SetTitle("Counts");
        ptrHistAllPeaks[i]->Draw("histsame");
    }

    auto legend = new TLegend(0.7,0.7,0.9,0.9);

    legend->AddEntry(ptrHistAllPeaks[0],"HV = 34.00 V","l");
    legend->AddEntry(ptrHistAllPeaks[1],"HV = 35.00 V","l");
    legend->AddEntry(ptrHistAllPeaks[2],"HV = 36.00 V","l");

    legend->Draw();

    // cout number of entries
    cout<<"Number of entries hist 34 V: "<<ptrHistAllPeaks[0]->Integral()<<endl;
    cout<<"Number of entries hist 35 V: "<<ptrHistAllPeaks[1]->Integral()<<endl;
    cout<<"Number of entries hist 36 V: "<<ptrHistAllPeaks[2]->Integral()<<endl;
}

//------------------------------------------------------------------------------
void Ana_LED(string file1, int last_event_n){
    //VARIABLES:
    //TRUE:
    average = true;
    find_offset_bool = true;
    find_peak_in_a_selected_window = true;
    find_peaks_bool = true;
    DO_NOT_DELETE_HIST_LED = true;

    find_charge_window_bool = true;

    led_and_dcr_0pe = true;

    smooth_trace_bool = true;

    min_peak_window[0] = minLED_amp;
    max_peak_window[0] = maxLED_amp;

    min_peak_window[1] = dcr_mintp;
    max_peak_window[1] = dcr_maxtp;

    find_area_trace_bool = true;



    // DISPLAY
    bool display_trace_LED = false;
    display_peaks_now = true;
    display_peaks = true;
    show_peak_LED_bool = true;
    line_bool = true;

    nfile = 0; //I only consider 1 file
    ptrHistAllPeaks[0]  = new TH1D("histAllPeaks","",bins_DCR,0,maxyhistAllPeaks);

    TCanvas *c = new TCanvas("Trace","Trace",w,h);


    //Analysis
    Analysis(file1, last_event_n, display_trace_LED, c);

    // LED
    TCanvas *canvLED = new TCanvas("canvLED", "canvLED", w,h);
    canvLED->cd();
    canvLED-> SetGrid();
    ptrHistLED->GetXaxis()->SetTitle("mV");
    ptrHistLED->GetYaxis()->SetTitle("Counts");
    ptrHistLED->Draw();
    canvLED->Update();

    // fit hist LED
    bool evaluate_cross_talk = true;
    // fit_hist_peaks_0pe_1pe_2pe(canvLED, ptrHistLED);
    cout<<endl;
    cout<<"-------------"<<endl;
    cout<<"---[ LED ]---"<<endl;
    cout<<"-------------"<<endl;
    fit_hist_peaks_gaus_sum_012(canvLED, ptrHistLED, evaluate_cross_talk);
    cout<<"---------------------------------------------------------------------"<<endl;

    if(led_and_dcr_0pe){
      // DCR from WINDOW
      TCanvas *canvDCR_window = new TCanvas("canvDCR_window", "canvDCR_window", w,h);
      canvDCR_window->cd();
      canvDCR_window-> SetGrid();
      ptrHistDCR_window->GetXaxis()->SetTitle("mV");
      ptrHistDCR_window->GetYaxis()->SetTitle("Counts");
      ptrHistDCR_window->Draw();
      canvDCR_window->Update();

      // fit hist ptrHistDCR_window
      evaluate_cross_talk = false;
      cout<<endl;
      cout<<"--------------"<<endl;
      cout<<"---[ DARK ]---"<<endl;
      cout<<"--------------"<<endl;
      fit_hist_peaks_gaus_sum_012(canvDCR_window, ptrHistDCR_window, evaluate_cross_talk);
      cout<<"---------------------------------------------------------------------"<<endl;

    }


    // Area LED
    TCanvas *canvArea = new TCanvas("canvArea", "canvArea", w,h);
    ptrHistArea->Draw();

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

    smooth_trace_bool = true;

    //Charge Pedestal:
    find_charge_window_bool = true;
    min_peak_window[0] = 200; //ns
    max_peak_window[0] = 250; //ns



    nfile = 0; //I only consider 1 file

    TCanvas *c = new TCanvas("Trace","Trace",w,h);

    //Analysis
    Analysis(file1, last_event_n, false, c);

    new TCanvas();
    ptrAllTrace->Draw();

    new TCanvas();
    offset_hist->Draw();

    new TCanvas();
    ptrHistCharge->Draw();

}



//------------------------------------------------------------------------------
void DCR_CT_1SiPM_1HV(string file1, int last_event_n){
    //TRUE:
    find_peaks_bool = true;
    DCR_from_cnt_bool = true;
    DO_NOT_DELETE_HIST_LED = true;
    DCR_DELAYS_bool = true;
    fit_hist_del_bool = false;

    smooth_trace_bool = true;

    // Print the name of the file:
    cout<<"Analysing file:"<<endl;
    cout<<"file = "<<file1<<endl;

    nfile = 0;

    TCanvas *c = new TCanvas("Trace","Trace",w,h);

    char h1[20], h2[20], h3[20];
    char k_temp[2] = "1";
    sprintf(h1, "histDelays");
    sprintf(h2, "histAllPeaks");
    sprintf(h3, "histDCRthr");

    //new hists:
    ptrHistDelays[0]   = new TH1D(strcat(h1,k_temp),"",bins_Delays,0,maxyhistDelays);
    ptrHistAllPeaks[0] = new TH1D(strcat(h2,k_temp),"",bins_DCR,0,maxyhistAllPeaks);
    ptrHistDCRthr[0]   = new TH1D(strcat(h3,k_temp),"",bins_DCR,0,maxyhistDCR);

    // colors:
    color_file[0] = kGreen+1;

    opacity = 0.3;


    // LOOP ON FILES:

    ptrHistAllPeaks[0]->Reset();
    ptrHistDelays[0]->Reset();
    DCR_cnt = 0;
    DCR_from_cnt = 0;
    trace_time = 0;
    n_ev_tot = 0;

    first_time_main_called = true; //will be set to false after the Analysis function is called
    ind_peaks_all = 0;
    for(int j=0; j<max_peaks; j++){
        peaks_all[j] = 0;
    }
    // DCR_func
    gDCR[0] = DCR_func_NO_Delays(file1,last_event_n, 1, c);
    delete[] gDCR_temp;



    TMultiGraph *DCR_mg = new TMultiGraph("DCR_mg", ";THR (mV); $\\frac{DCR}{mm^2}$ (Hz)");
    DCR_mg->Add(gDCR[0]);


    TCanvas *cDCR_loop = new TCanvas("cDCR_loop", "cDCR_loop");

    cDCR_loop->SetGrid();
    cDCR_loop->SetLogy();
    DCR_mg->Draw("A3L");

    // 0.5 pe
    cout<<"     DCR at 0.5 pe from cnt    = ("<<DCR_pe_0_5_Area_vect[0]*n6*Area<<" +- "<<errDCR_pe_0_5_Area_vect[0]*n6*Area<<") MHz"<<endl;
    cout<<"     DCR at 0.5 pe from delays = ("<<DCR_pe_0_5_Area_delays_vect[0]*n6*Area<<" +- "<<errDCR_pe_0_5_Area_delays_vect[0]*n6*Area<<") MHz"<<endl;

    // 1.5 pe
    cout<<"     DCR at 1.5 pe from cnt = ("<<DCR_pe_1_5_Area_vect[0]*n6*Area<<" +- "<<errDCR_pe_1_5_Area_vect[0]*n6*Area<<") MHz"<<endl;
    cout<<"     DCR at 1.5 pe from delays = ("<<DCR_pe_1_5_Area_delays_vect[0]*n6*Area<<" +- "<<errDCR_pe_1_5_Area_delays_vect[0]*n6*Area<<") MHz"<<endl;

    cout<<endl;


}



//------------------------------------------------------------------------------
void DCR_CT_1SiPM_nHVs(string filelist, int nfile_in_list, int last_event_n){
    //TRUE:
    find_peaks_bool = true;
    DCR_from_cnt_bool = true;
    DO_NOT_DELETE_HIST_LED = true;
    DCR_DELAYS_bool = true;
    fit_hist_del_bool = false;

    smooth_trace_bool = true;

    // OPEN FILE AND READ FILE LIST
    ifstream OpenFile (filelist.c_str());
    string file[nfilemax], legend_entry[nfilemax];
    nfiletot=0;

    while(nfiletot<nfile_in_list and nfiletot<nfilemax){
        OpenFile>>file[nfiletot];
        OpenFile>>legend_entry[nfiletot];
        OpenFile>>legend_entry[nfiletot];
        OpenFile>>legend_entry[nfiletot];
        nfiletot++;
    }
    OpenFile.close();

    // Print the name of the files:
    cout<<"Analysing files:"<<endl;
    for(int i=0; i<nfiletot; i++){
        cout<<"file "<<i+1<<" = "<<file[i]<<endl;
    }

    TCanvas *c = new TCanvas("Trace","Trace",w,h);

    for(int k=0; k<nfiletot; k++){
        //In order to set n different titles:
        char h1[20], h2[20], h3[20];
        char k_temp[2];
        sprintf(h1, "histDelays");
        sprintf(h2, "histAllPeaks");
        sprintf(h3, "histDCRthr");
        sprintf(k_temp, "%d", k);

        //new hists:
        ptrHistDelays[k]   = new TH1D(strcat(h1,k_temp),"",bins_Delays,0,maxyhistDelays);
        ptrHistAllPeaks[k] = new TH1D(strcat(h2,k_temp),"",bins_DCR,0,maxyhistAllPeaks);
        ptrHistDCRthr[k]   = new TH1D(strcat(h3,k_temp),"",bins_DCR,0,maxyhistDCR);
    }

    // colors:
    color_file[0] = kBlack;
    color_file[1] = kOrange+5;
    color_file[2] = kRed;
    color_file[3] = kOrange;
    color_file[4] = kGreen+1;
    color_file[5] = kCyan;
    color_file[6] = kBlue;
    color_file[7] = kViolet;
    color_file[8] = kGray;
    color_file[9] = kGray+3;

    opacity = 0.3;

    nfile = 0;

    // LOOP ON FILES:
    for(int i=0; i<nfiletot; i++){

        // reset
        ptrHistAllPeaks[0]->Reset();
        ptrHistDelays[0]->Reset();
        DCR_cnt = 0;
        DCR_from_cnt = 0;
        trace_time = 0;
        n_ev_tot = 0;

        first_time_main_called = true; //will be set to false after the Analysis function is called
        ind_peaks_all = 0;
        for(int j=0; j<max_peaks; j++){
            peaks_all[j] = 0;
        }
        // DCR_func
        gDCR[i] = DCR_func_NO_Delays(file[i],last_event_n, nfiletot, c);
        delete[] gDCR_temp;

        nfile++;
    }


    TMultiGraph *DCR_mg = new TMultiGraph("DCR_mg", ";THR (mV); $\\frac{DCR}{mm^2}$ (Hz)");
    for(int i=0; i<nfiletot; i++){
        DCR_mg->Add(gDCR[i]);
    }

    TCanvas *cDCR_loop = new TCanvas("cDCR_loop", "cDCR_loop");

    cDCR_loop->SetGrid();
    cDCR_loop->SetLogy();
    DCR_mg->Draw("A3L");

    auto legendDCR_loop = new TLegend(0.75,0.75,0.9,0.9);
    for(int i=0; i<nfiletot; i++){
        legendDCR_loop->AddEntry(gDCR[i],legend_entry[i].c_str(),"l");
    }
    legendDCR_loop->SetNColumns(3);
    legendDCR_loop->Draw();



    for(int i=0; i<nfiletot; i++){
        // nfile
        cout<<legend_entry[i]<<endl;

        // 0.5 pe
        cout<<"     DCR at 0.5 pe from cnt    = ("<<DCR_pe_0_5_Area_vect[i]*n6*Area<<" +- "<<errDCR_pe_0_5_Area_vect[i]*n6*Area<<") MHz"<<endl;
        cout<<"     DCR at 0.5 pe from delays = ("<<DCR_pe_0_5_Area_delays_vect[i]*n6*Area<<" +- "<<errDCR_pe_0_5_Area_delays_vect[i]*n6*Area<<") MHz"<<endl;

        // 1.5 pe
        cout<<"     DCR at 1.5 pe from cnt = ("<<DCR_pe_1_5_Area_vect[i]*n6*Area<<" +- "<<errDCR_pe_1_5_Area_vect[i]*n6*Area<<") MHz"<<endl;
        cout<<"     DCR at 1.5 pe from delays = ("<<DCR_pe_1_5_Area_delays_vect[i]*n6*Area<<" +- "<<errDCR_pe_1_5_Area_delays_vect[i]*n6*Area<<") MHz"<<endl;

        cout<<endl;
    }

}



//------------------------------------------------------------------------------
//-------------------------[   SECONDARY FUNCTIONS   ]--------------------------
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
void Analysis(string file, int last_event_n, bool display, TCanvas *c){
    gROOT->Reset();

    if(first_time_main_called){
        if(Agilent_MSO6054A or Digitizer_CAEN)
            Read_Agilent_CAEN(file, last_event_n, display, c);
        if(DRS4_Evaluation_Board)
            ReadBin(file, last_event_n, display, c);
        if(DRS4_Evaluation_Board_Mod)
          ReadRootFile(file, last_event_n, display, c);
    }else{ //in order to speed up when I want to find peaks at different thresholds: I read the file only one time
        FindDelaysFromVector();
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

//***** AVERAGE INTEGRAL
    if(average){
        double sum=0;
        for( int j=0; j<trace_length; j++){
            sum+=trace_AVG[1][j];
        }

        cout << "Integral of the average is: " << sum << endl;
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
                fit_hist_peaks(cAllPeaks1, ptrHistAllPeaks[nfile]);
            }

        }
        // if(drawHistAllPeaksAll){
        //     cAllPeaks->cd();
        //
        //     if(nfile==1) ptrHistAllPeaks[nfile]->SetLineColor(kGreen+1);
        //     if(nfile==2) ptrHistAllPeaks[nfile]->SetLineColor(kRed+1);
        //     ptrHistAllPeaks[nfile]->Draw("histsame");
        //
        //     auto legend = new TLegend(0.7,0.7,0.9,0.9);
        //     for(int k=0; k<nfilemax; k++){
        //         if(nfile==0) legend->AddEntry(ptrHistDCRthr[nfile],"HV = 34.00 V","l");
        //         if(nfile==1) legend->AddEntry(ptrHistDCRthr[nfile],"HV = 35.00 V","l");
        //         if(nfile==2) legend->AddEntry(ptrHistDCRthr[nfile],"HV = 36.00 V","l");
        //     }
        //     legend->Draw();
        // }
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
        if(fit_hist_del_bool){
            fit_hist_del(expDelLow_max, expDelHigh_max);
        }

    }

    if(showHist_bool){
        TCanvas *cHist = new TCanvas("hist_GAIN","hist_GAIN",w,h);
        cHist->SetGrid();
        cHist->cd();
        if(SetLogyHist) cHist->SetLogy();
        ptrHistLED->Draw("hist");
    }



    if(!DO_NOT_DELETE_HIST_LED)delete ptrHistLED;

    //since this function is called, I set first_time_main_called to false
    first_time_main_called = false;

}

//------------------------------------------------------------------------------
TGraphErrors *DCR_func(string file1, int last_event_n, int tot_files, TCanvas *c){

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
        DCR_Area = new float*[tot_files];
            for(int i = 0; i < tot_files; i++) {
                DCR[i] = new float[n_DCR];
        }
        errDCR_Area = new float*[tot_files];
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
        Analysis(file1,last_event_n,display, c);
        Get_DCR_temp_and_errDCR_temp();
        ptrHistDelays[nfile]->Reset();
        DCR_thr[h] = thr_to_find_peaks; //mV
        DCR[nfile][h] = DCR_temp[nfile];
        errDCR[nfile][h] = errDCR_temp[nfile];
        thr_to_find_peaks_vect[nfile][h] = thr_to_find_peaks;
        thr_to_find_peaks = thr_to_find_peaks + gap_between_thr; //I jump to the next thr in order to evaluate the new DRC
        h++;
    }

    if(automatic_find_thr_1pe_2pe)  find_DCR_0_5_pe_and_1_5_pe_auto();
    else                            find_DCR_0_5_pe_and_1_5_pe_manual();


    // find DCR / Area
    DivideAreaDCR_DCR_0_5_DCR_1_5();


    TGraphErrors *gDCR = new TGraphErrors(n_DCR, DCR_thr, DCR_Area[nfile],NULL, errDCR_Area[nfile]);

     gDCR->SetLineWidth(2);


    if(nfile==0){gDCR->SetLineColor(color_file_1);  gDCR->SetFillColorAlpha(color_file_1, opacity);}
    if(nfile==1){gDCR->SetLineColor(color_file_2);  gDCR->SetFillColorAlpha(color_file_2, opacity);   }
    if(nfile==2){gDCR->SetLineColor(color_file_3);  gDCR->SetFillColorAlpha(color_file_3, opacity);  }

    cout<<endl<<endl;

    return gDCR;

}


//------------------------------------------------------------------------------
TGraphErrors *DCR_func_NO_Delays(string file1, int last_event_n, int tot_files, TCanvas *c){
    bool display = false;
    DCR_DELAYS_bool = true;

    thr_to_find_peaks = min_thr_to_find_peaks;

    n_DCR = (int)((max_thr_to_find_peaks - min_thr_to_find_peaks)/gap_between_thr);

    if(first_time_DCR_called){ // first_time_DCR_called
        DCR = new float*[tot_files];
            for(int i = 0; i < tot_files; i++) {
                DCR[i] = new float[n_DCR];
        }
        errDCR = new float*[tot_files];
            for(int i = 0; i < tot_files; i++) {
                errDCR[i] = new float[n_DCR];
        }
        DCR_Area = new float*[tot_files];
            for(int i = 0; i < tot_files; i++) {
                DCR_Area[i] = new float[n_DCR];
        }
        errDCR_Area = new float*[tot_files];
            for(int i = 0; i < tot_files; i++) {
                errDCR_Area[i] = new float[n_DCR];
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

    } // end first_time_DCR_called

    first_time_DCR_called = false;

    float *DCR_thr = new float[n_DCR];
    double delta = (pe_1_5_vect[nfile] - pe_0_5_vect[nfile])/n_DCR;

    /////////////////////////////////////////////
    ///     RESET
    /////////////////////////////////////////////
    ind_peaks_all = 0;
    ind_peaks_all_delay = 0;

    for(int i=0; i<max_peaks; i++){
        peaks_all_delay[0][i] = 0;
        peaks_all_delay[1][i] = 0;

        peaks_all[i] = 0;
        peaks_all_time[i] = 0;
    }


    int h = 0;

    while(thr_to_find_peaks <= max_thr_to_find_peaks){ // loop for DCR vs thresholds graph
        // reset
        DCR_cnt = 0;
        DCR_from_cnt = 0;
        trace_time = 0;


        if(h==0){
            Analysis(file1,last_event_n,display, c);
        }else{
            FindDCRfromVector();
        }


        DCR_thr[h] = thr_to_find_peaks; //mV
        DCR[nfile][h] = DCR_from_cnt;
        errDCR[nfile][h] = errDCR_from_cnt;


        // 0.5 pe
        if((thr_to_find_peaks > pe_0_5_vect[nfile] - delta) and (thr_to_find_peaks < pe_0_5_vect[nfile] + delta)){
            // reset
            ptrHistAllPeaks[nfile]->Reset();
            ptrHistDelays[nfile]->Reset();

            // From delays
            FindDelaysFromVector();
            fit_hist_del(expDelLow_max, expDelHigh_max);

            DCR_pe_0_5_delays_vect[nfile] = GetDCRfromDelays();
            errDCR_pe_0_5_delays_vect[nfile] = GetErrDCRfromDelays();

            // From cnt
            DCR_pe_0_5_vect[nfile] = DCR[nfile][h];
            errDCR_pe_0_5_vect[nfile] = errDCR[nfile][h];

            // cout<<"DCR_thr[h] "<<DCR_thr[h]<<endl;
            // cout<<"DCR_pe_0_5_delays_vect[nfile]"<<DCR_pe_0_5_delays_vect[nfile]<<endl;
        }

        // 1.5 pe
        if((thr_to_find_peaks > pe_1_5_vect[nfile] - delta) and (thr_to_find_peaks < pe_1_5_vect[nfile] + delta)){
            // reset
            ptrHistAllPeaks[nfile]->Reset();
            ptrHistDelays[nfile]->Reset();

            // From delays
            FindDelaysFromVector();
            fit_hist_del(expDelLow_max, expDelHigh_max);

            DCR_pe_1_5_delays_vect[nfile] = GetDCRfromDelays();
            errDCR_pe_1_5_delays_vect[nfile] = GetErrDCRfromDelays();

            // From cnt
            DCR_pe_1_5_vect[nfile] = DCR[nfile][h];
            errDCR_pe_1_5_vect[nfile] = errDCR[nfile][h];

            // cout<<"DCR_thr[h] "<<DCR_thr[h]<<endl;
            // cout<<"DCR[nfile][h] "<<DCR[nfile][h]<<endl;
            // cout<<"DCR_pe_1_5_vect[nfile] "<<DCR_pe_1_5_vect[nfile]<<endl;
        }

        thr_to_find_peaks_vect[nfile][h] = thr_to_find_peaks;
        thr_to_find_peaks += gap_between_thr; //I jump to the next thr in order to evaluate the new DRC
        h++;
    }

    // find DCR / Area
    DivideAreaDCR_DCR_0_5_DCR_1_5();



    TGraphErrors *gDCR_temp = new TGraphErrors(n_DCR, DCR_thr, DCR_Area[nfile],NULL, errDCR_Area[nfile]);

    gDCR_temp->SetLineWidth(2);

    gDCR_temp->SetLineColor(color_file[nfile]);
    gDCR_temp->SetFillColorAlpha(color_file[nfile], opacity);


    cout<<endl<<endl;

    delete[] DCR_thr;

    return gDCR_temp;

}



//------------------------------------------------------------------------------
void discriminator(double thr, float **t, double length){
    int i = 0;
    bool find_rising_edge = true;
    discriminator_cnt = 0;
    while(i<length-4){ // loop on trace t

        if(find_rising_edge){ // find a new peak
            if(t[1][i] > thr){ // loop > thr
                if((t[1][i]<t[1][i+1]) and (t[1][i+1]<t[1][i+2])){ // rising edge
                    discriminator_cnt ++;
                    i+=2;
                    find_rising_edge = false;
                }
                else{
                    i++;
                } // end rising edge

            }else{
                i++;
            } // end loop > thr
        } // end find a new peak
        else{ // wait until the peak ends
            if((t[1][i]>t[1][i+1]) and (t[1][i+1]>t[1][i+2])){ // falling edge
                find_rising_edge = true;
                i++;
            } // end falling edge
            else{
                i++;
            }
        } //end wait until the peak ends



    } // end loop on trace t
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
    index_func=-1;
    for( ii=mintp; ii<maxtp; ii++){
        //cout<<mintp<<"\t"<<maxtp<<endl;
        if(trace_DLED[1][ii]>max_func){
            max_func=trace_DLED[1][ii];
            index_func=ii;
        }
    }

    return index_func;

    // if(trace_DLED[1][index_func]<thr_to_find_peaks){
    //       return index_func;
    // }
    // else{
    //     if(trace_DLED[1][index_func-1]<max_func &&
    //        trace_DLED[1][index_func-2]<max_func &&
    //        trace_DLED[1][index_func-3]<max_func &&
    //        trace_DLED[1][index_func+1]<max_func &&
    //        trace_DLED[1][index_func+2]<max_func &&
    //        trace_DLED[1][index_func+3]<max_func
    //       )
    //     {
    //        return index_func;
    //     }
    //     else{
    //         return -1;
    //     }
    // }
}

//------------------------------------------------------------------------------
void find_peaks(float thr_to_find_peaks, int max_peak_width, int min_peak_width,int blind_gap,  bool DCR_DELAYS_bool){ //I look for every peaks in the trace, only if I'm in DARK mode
    ii=2*GSPS;
    int before_ind = (int)2*GSPS;
    int index_peak;
    DCR_cnt_temp = 0;
    int index_old = 0;
    int index_new = 0;
    int peak_width = max_peak_width;
    num_peaks=0;
    for(int i=0; i<max_peak_num; i++)index_vect[i]=0;
    ind_peaks_all_delay = 0;

    while(ii<trace_DLED_length){//I find peaks after the DLED procedure

        //I only consider points above thr_to_find_peaks on the rising edge
        if((trace_DLED[1][ii]>thr_to_find_peaks) and (trace_DLED[1][ii-before_ind]<thr_to_find_peaks) and (ii+max_peak_width<trace_DLED_length)){ // loop > thr and rising edge
            DCR_cnt_temp++; //I've seen a peak; if I'm in dark mode it's DCR

            //Now I want to see the peak width
            for(int k=ii+min_peak_width; (k<ii+max_peak_width); k++){//I open a window and look for a point below thr
                if(trace_DLED[1][k]<thr_to_find_peaks){
                    peak_width = k-ii;
                    break;
                }
            }

            //Now I look for the peak in that window
            if(ii+peak_width<trace_DLED_length)
                index_new = find_peak_fix_time(ii, ii+peak_width);
            else
                index_new = find_peak_fix_time(ii, trace_DLED_length);

            if((index_new-index_old)>blind_gap){ // start loop blind_gap
                // I fill the hist of all the peaks
                if(fill_hist_peaks_when_found){
                    ptrHistAllPeaks[nfile]->Fill(trace_DLED[1][index_new]);
                }

                // DCR from the delay (Itzler Mark - Dark Count Rate Measure (pag 5 ss))
                if(DCR_DELAYS_bool){
                    // I fill the hist with delays between 1 pe peaks (or higher):
                    if(index_old>0){
                        ptrHistDelays[nfile] -> Fill(index_new - index_old);
                        peaks_all_delay[0][ind_peaks_all_delay] = index_new;
                        peaks_all_delay[1][ind_peaks_all_delay] = trace_DLED[1][index_new];
                        ind_peaks_all_delay++;
                    }
                }

                // now index_new will be index_old, for the following loop cycle
                index_old = index_new;

                if(num_peaks<max_peak_num){
                    index_vect[num_peaks] = index_new;
                    num_peaks++;
                }

                ii=ii+peak_width;
            }else{ // I'm in the blind_gap
                index_old = index_new;
                ii++;
            } // end loop blind_gap
        }else{ // point NOT > thr and rising edge
                ii++;
        } // end loop > thr and rising edge
    }
    peaks_all_delay[0][ind_peaks_all_delay] = -1;
    peaks_all_delay[1][ind_peaks_all_delay] = -1;
    ind_peaks_all_delay++;
}

//------------------------------------------------------------------------------
void FindPeaksRisingFalling(double thr, float **t, double length, int max_peak_width, int rising_points, int falling_points){
    int i=1+rising_points;
    int index_old = -1;
    int index_new = -1;
    int peak_start = 0;
    int peak_end = 0;
    int peak_width = max_peak_width;
    double time_delay = 0;
    num_peaks=0;
    bool rising = true;
    bool falling = true;
    int half_rising_points = 0;
    int half_falling_points = 0;
    int diff = 0;

    for(int i=0; i<max_peak_num; i++){
            index_vect[i]=0;
    }

    bool find_rising_edge = true;

    DCR_cnt_temp = 0;

    while(i<length-1){ // loop on trace t

        if(find_rising_edge){ // find rising edge
            if(t[1][i] > thr){ // loop > thr

                // find if rising edge or not
                rising = true;
                half_rising_points = (int)(rising_points/2);
                diff = rising_points%2;
                // cout<<endl<<"rising_points = "<<rising_points<<", i = "<<i<<endl;
                for(int j=i-half_rising_points; j<i+half_rising_points+diff; j++){ // rising?
                    // cout<<j<<", ";
                    if( t[1][j] >= t[1][j+1] ){
                        rising = false;
                    }
                } // end rising?

                // getchar();

                if(rising){ // rising edge
                    peak_start = i-1;
                    i++;
                    find_rising_edge = false;
                }
                else{
                    i++;
                } // end rising edge

            }else{
                i++;
            } // end loop > thr
        } // end find rising edge
        else{ // wait until the peak ends

            // find if falling edge or not
            falling = true;
            half_falling_points = (int)(falling_points/2);
            diff = falling_points%2;
            // cout<<endl<<"falling_points = "<<falling_points<<", i = "<<i<<endl;
            for(int j=i-half_falling_points; j<i+half_falling_points+diff; j++){ // falling?
                // cout<<j<<", ";
                if( t[1][j] <= t[1][j+1] ){
                    falling = false;
                }
            } // end falling?

            if(falling){ // falling edge
                peak_end = i+1;
                peak_width = peak_end-peak_start;

                // I have found a falling edge. The following time I will look for a new peak:
                find_rising_edge = true;


                // I have to check if this can be a peak:
                if(peak_width < max_peak_width){ // peak_width < max_peak_width
                    // I have found a new peak
                    DCR_cnt++;
                    DCR_cnt_temp++;

                    //Now I look for the maximum value in that window
                    index_new = find_peak_fix_time(peak_start, peak_end);

                    // I fill the hist of all the peaks
                    if(fill_hist_peaks_when_found){
                        ptrHistAllPeaks[nfile]->Fill(trace_DLED[1][index_new]);
                    }

                    // I fill also the vector:
                    peaks_all_time[ind_peaks_all] = t[0][index_new];
                    peaks_all[ind_peaks_all] = t[1][index_new];
                    ind_peaks_all++;

                    // DCR from the delay (Itzler Mark - Dark Count Rate Measure (pag 5 ss))
                    if(DCR_DELAYS_bool){ // DCR from delays
                        // I fill the hist with delays between 1 pe peaks (or higher):
                        if(index_old>0){
                            time_delay = t[0][index_new] - t[0][index_old];
                            // if( (time_delay>expDelLow_max) and (time_delay<expDelHigh_max) and (t[0][index_new]>2*expDelLow_max) and (index_new<900) ){
                            if( (time_delay>expDelLow_max) and (time_delay<expDelHigh_max) ){
                                ptrHistDelays[nfile] -> Fill(time_delay);
                                peaks_all_delay[0][ind_peaks_all_delay] = t[0][index_new];
                                peaks_all_delay[1][ind_peaks_all_delay] = t[1][index_new];
                                ind_peaks_all_delay++;
                            }
                        }
                    } // end DCR from delays

                    // now index_new will be index_old, for the following loop cycle
                    index_old = index_new;

                    if(num_peaks<max_peak_num){
                        index_vect[num_peaks] = index_new;
                        num_peaks++;
                    }

                    i+=gap_between_peaks;
                } // end peak_width < max_peak_width
                else{ // this is not a peak, too wide in time
                    i++;
                    // cout<<"*******   PEAK TOO WIDE   *******"<<endl;
                }
            } // end falling edge
            else{ // I have not found a falling edge
                i++;
            }
        } //end wait until the peak ends

    } // end loop on trace t

    // in peaks_all_delay I put -1 at the end of the trace
    peaks_all_delay[0][ind_peaks_all_delay] = -1;
    peaks_all_delay[1][ind_peaks_all_delay] = -1;

    peaks_all[ind_peaks_all] = -1;
    peaks_all[ind_peaks_all] = -1;

    ind_peaks_all_delay++;

}

//------------------------------------------------------------------------------
void FindPeakPositions(float* vector, Bool_t dled_bool, Int_t dt)
{
    std::vector<int> vec(vector, vector + sizeof vector / sizeof vector[0]);
   // Routine to find peaks. Code adapted from:
   // https://stackoverflow.com/questions/22583391/peak-signal-detection-in-realtime-timeseries-data/46956908#46956908
   // (original Q/A: https://stackoverflow.com/questions/22583391/peak-signal-detection-in-realtime-timeseries-data)

   //lag 5 for the smoothing functions
   UInt_t lag = 20;
   //3.5 standard deviations for signal
   Float_t threshold = 1;
   //between 0 and 1, where 1 is normal influence, 0.5 is half
   Float_t influence = 1.;

   if (vec.size() <= lag + 2)
   {
      std::vector<Int_t> emptyVec;
      peak_pos = emptyVec;
   }

   // DLED stuff

   UInt_t blind_gap = 0;

   if(dled_bool){
      blind_gap = 2*dt;
   }

   //Initialise variables
   std::vector<Int_t> signals(vec.size(), 0.0);
   std::vector<Float_t> filteredY(vec.size(), 0.0);
   std::vector<Float_t> avgFilter(vec.size(), 0.0);
   std::vector<Float_t> stdFilter(vec.size(), 0.0);
   std::vector<Float_t> subVecStart(vec.begin(), vec.begin() + lag);

   avgFilter[lag] = GetMean(subVecStart);
   stdFilter[lag] = GetStdDev(subVecStart);

   UInt_t old_peak = 0;

   for (size_t i = lag + 1; i < vec.size()-2; i++)
   {
      if (std::abs(vec[i] - avgFilter[i - 1]) > threshold * stdFilter[i - 1])
      {
         if (vec[i] > avgFilter[i - 1])
         {
            if ((vec[i-2] < vec[i]) && (vec[i-1] < vec[i]) && (vec[i+1] < vec[i]) && (vec[i+2] < vec[i]))
            {
               if (i - old_peak > blind_gap)
               {
                  signals[i] = 1; //# Positive signal
                  old_peak = i;
               }
               else
               {
                  signals[i] = 0;
               }

            }
            else
            {
               signals[i] = 0;
            }
         }
         else
         {
            if ((vec[i-2] > vec[i]) && (vec[i-1] > vec[i]) && (vec[i+1] > vec[i]) && (vec[i+2] > vec[i]))
            {
               signals[i] = -1; //# Negative signal
            }
            else
            {
               signals[i] = 0;
            }
         }
         //Make influence lower
         filteredY[i] = influence* vec[i] + (1 - influence) * filteredY[i - 1];
      }
      else
      {
         signals[i] = 0; //# No signal
         filteredY[i] = vec[i];
      }

      //Adjust the filters
      std::vector<Float_t> subVec(filteredY.begin() + i - lag, filteredY.begin() + i);
      avgFilter[i] = GetMean(subVec);
      stdFilter[i] = GetStdDev(subVec);
   }

   peak_pos = signals;
}


//-----------------------------------------------------------------------------
void FindPeaksFromPositions(){

   // get only positive peaks and create two vector containing the
   // times and amplitudes of the peaks



   int index_old = -1;
   int index_new = -1;
   double time_delay = 0;

   for(UInt_t i = 0; i < peak_pos.size(); i++){
      if(peak_pos[i] != 1){
         continue;
      }
      else // I have a peak at the i-th position
      {
          index_new = i;

          // I fill the hist of all the peaks
          if(fill_hist_peaks_when_found){
              ptrHistAllPeaks[nfile]->Fill(trace_DLED[1][index_new]);
          }

          // DCR from the delay (Itzler Mark - Dark Count Rate Measure (pag 5 ss))
          if(DCR_DELAYS_bool){
              // I fill the hist with delays between 1 pe peaks (or higher):
              if(index_old>0){
                  time_delay = trace_DLED[0][index_new] -trace_DLED[1][index_old];
                  ptrHistDelays[nfile] -> Fill(time_delay);
                  peaks_all_delay[0][ind_peaks_all_delay] = trace_DLED[0][index_new];
                  peaks_all_delay[1][ind_peaks_all_delay] = trace_DLED[1][index_new];
                  ind_peaks_all_delay++;
              }
          }

          // now index_new will be index_old, for the following loop cycle
          index_old = index_new;

          if(num_peaks<max_peak_num){
              index_vect[num_peaks] = index_new;
              num_peaks++;
          }

      }
   }

}

//------------------------------------------------------------------------------
void FindDCRfromVector(){

    for(int i=0; i<ind_peaks_all; i++){
        if(peaks_all[i] > thr_to_find_peaks){
            DCR_cnt++;
        }
    }

    // trace_time = trace_time_raw - DCR_cnt * 2 * dleddt - n_ev_tot * rise_time;
    trace_time = trace_time_raw - DCR_cnt * 2 * (dleddt) - n_ev_tot * rise_time;


    double err_DCR_cnt = 0.;
    double err_trace_time = 0.;

    // DCR from cnt
    trace_time *= TMath::Power(10,-9); // trace time is in ns
    trace_time /= n_ev_tot;
    DCR_from_cnt = (double)DCR_cnt / (trace_time * n_ev_tot);

    // errors
    err_trace_time = TMath::Sqrt( (5*n_ev_tot) * TMath::Power(10,-9) * TMath::Power(10,-9) );
    err_DCR_cnt = TMath::Sqrt((double)DCR_cnt);
    errDCR_from_cnt = TMath::Sqrt( TMath::Power( err_DCR_cnt / (trace_time * n_ev_tot), 2) + TMath::Power( (double)DCR_cnt*err_trace_time / ( trace_time * n_ev_tot * trace_time * n_ev_tot ), 2) );

}

//------------------------------------------------------------------------------
void FindDelaysFromVector(){
    //int n,i,k;
    double time_new = -1., time_old = -1.;
    double time_delay = 0.;

    for(int i=0;i<ind_peaks_all; i++){ // loop on peaks_all_delay
        // cout<<time_old<<endl;
        // cout<<peaks_all_delay[0][i]<<endl;
        if(peaks_all_time[i]==-1){ // new trace
            time_old = -1;
            time_new = -1;
        }
        else{ // trace
            if(peaks_all[i] > thr_to_find_peaks){ // > thr
                time_new = peaks_all_time[i];
                if(time_old>0){ // not the first peak trace
                    time_delay = time_new - time_old;
                    // if( (time_delay>expDelLow_max) and (time_delay<expDelHigh_max) and (t[0][time_new]>2*expDelLow_max) and (time_new<900) ){
                    if( (time_delay>expDelLow_max) and (time_delay<expDelHigh_max) ){
                        ptrHistDelays[nfile] -> Fill(time_delay);
                    }
                } // end not the first peak in trace
                time_old = time_new;

            } // end > thr
        } // end trace

    } // end loop on peaks_all_delay

}

//------------------------------------------------------------------------------
void DivideAreaDCR_DCR_0_5_DCR_1_5(){
    for(int i=0; i<n_DCR; i++){
        DCR_Area[nfile][i]    = DCR[nfile][i] / Area;
        errDCR_Area[nfile][i] = errDCR[nfile][i] / Area;
    }

    DCR_pe_0_5_Area_vect[nfile]    = DCR_pe_0_5_vect[nfile] / Area;
    errDCR_pe_0_5_Area_vect[nfile] = errDCR_pe_0_5_vect[nfile] / Area;
    DCR_pe_1_5_Area_vect[nfile]    = DCR_pe_1_5_vect[nfile] / Area;
    errDCR_pe_1_5_Area_vect[nfile] = errDCR_pe_1_5_vect[nfile] / Area;

    DCR_pe_0_5_Area_delays_vect[nfile] = DCR_pe_0_5_delays_vect[nfile] / Area;
    errDCR_pe_0_5_Area_delays_vect[nfile] = errDCR_pe_0_5_delays_vect[nfile] / Area;
    DCR_pe_1_5_Area_delays_vect[nfile] = DCR_pe_1_5_delays_vect[nfile] / Area;
    errDCR_pe_1_5_Area_delays_vect[nfile] = errDCR_pe_1_5_delays_vect[nfile] / Area;
}

//------------------------------------------------------------------------------
void show_trace(TCanvas* canv, float *x, float *y, int trace_length, float miny, float maxy, bool line_bool, bool delete_bool){

    if(reverse_bool){
      for(ii=0; ii<trace_length; ii++){
        y[ii] = -y[ii];
      }
    }

    TGraphErrors *graph = new TGraphErrors(trace_length,x,y,0,0);
    graph->SetTitle("Trace");
    graph->Draw("apl");
    graph->GetXaxis()->SetTitle("Time (ns)");
    graph->GetYaxis()->SetTitle("Amplitude (mV)");
    graph->GetYaxis()->SetTitleOffset(1.2);
    graph->GetXaxis()->SetTitleOffset(1.2);
    graph->GetYaxis()-> SetRangeUser(miny,maxy);
    graph->Draw("apl");
    if(line_bool){
        TLine *lmin = new TLine (minLED_amp, miny, minLED_amp, maxy);
        TLine *lmax = new TLine (maxLED_amp, miny, maxLED_amp, maxy);
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
    if(reverse_bool){
      for(ii=0; ii<trace_length1; ii++){
        y1[ii] = -y1[ii];
      }
    }
    canv->cd(1);
    TGraphErrors *graph1 = new TGraphErrors(trace_length1,x1,y1,0,0);
    graph1->SetName("graph1");
    graph1->SetTitle("Original trace");
    graph1->Draw("apl");
    graph1->GetXaxis()->SetTitle("Time (ns)");
    graph1->GetYaxis()->SetTitle("Amplitude (mV)");
    graph1->GetYaxis()->SetTitleOffset(1.2);
    graph1->GetXaxis()->SetTitleOffset(1.2);
    graph1->GetYaxis()-> SetRangeUser(miny1,maxy1);
    graph1->Draw("apl");
    //graph1->SetEditable(kFALSE);
   if(line_bool){
        TLine *lmin = new TLine (min_line, miny1, min_line, maxy1);
        TLine *lmax = new TLine (max_line, miny1, max_line, maxy1);
        lmin->SetLineColor(kBlue);
        lmax->SetLineColor(kBlue);
        lmin->Draw("plsame");
        lmax->Draw("plsame");
    }

    canv->cd(2);
    TGraphErrors *graph2 = new TGraphErrors(trace_length2,x2,y2,0,0);
    graph2->SetName("graph2");
    graph2->SetTitle("Trace after DLED");
    graph2->Draw("apl");
    graph2->GetXaxis()->SetTitle("Time (ns)");
    graph2->GetYaxis()->SetTitle("Amplitude (mV)");
    graph2->GetYaxis()->SetTitleOffset(1.2);
    graph2->GetXaxis()->SetTitleOffset(1.2);
    graph2->GetYaxis()-> SetRangeUser(miny2,maxy2);
    graph2->Draw("apl");
    //graph2->SetEditable(kFALSE);
if(line_bool){
        TLine *lmin = new TLine (min_line, miny2, min_line, maxy2);
        TLine *lmax = new TLine (max_line, miny2, max_line, maxy2);
        lmin->SetLineColor(kBlue);
        lmax->SetLineColor(kBlue);
        lmin->Draw("plsame");
        lmax->Draw("plsame");
    }

    TGraphErrors *graphPeaks = nullptr;
    TGraphErrors *graphPeaks_DLED = nullptr;
    TGraphErrors *graphPeak_LED = nullptr;
    TGraphErrors *graphPeak_LED_DLED = nullptr;

    if(display_peaks_now){
        float x_peaks[num_peaks], y_peaks[num_peaks], x_peaks_DLED[num_peaks], y_peaks_DLED[num_peaks];
        for(int i=0; i<num_peaks; i++){
            x_peaks[i] = trace[0][index_vect[i]+dleddt];
            y_peaks[i] = trace[1][index_vect[i]+dleddt];
            x_peaks_DLED[i] = trace_DLED[0][index_vect[i]];
            y_peaks_DLED[i] = trace_DLED[1][index_vect[i]];
        }

        //graphPeaks
        canv->cd(1);
        graphPeaks = new TGraphErrors(num_peaks,x_peaks,y_peaks,0,0);
        graphPeaks->SetName("graphPeaks");
        graphPeaks->SetTitle("graphPeaks");
        graphPeaks->SetMarkerStyle(20);
        graphPeaks->SetMarkerColor(kRed);
        graphPeaks->Draw("psame");
        //graphPeaks->SetEditable(kFALSE);


        //graphPeaks_DLED
        canv->cd(2);
        graphPeaks_DLED = new TGraphErrors(num_peaks,x_peaks_DLED, y_peaks_DLED,0,0);
        graphPeaks_DLED->SetName("graphPeaks_DLED");
        graphPeaks_DLED->SetTitle("graphPeaks_DLED");
        graphPeaks_DLED->SetMarkerStyle(20);
        graphPeaks_DLED->SetMarkerColor(kGreen+1);
        graphPeaks_DLED->Draw("psame");
        //graphPeaks_DLED->SetEditable(kFALSE);

    }

    if(show_peak_LED_bool && !peak_rejected){
      // On trace
      float x_g[1], y_g[1];
      x_g[0] = trace[0][index_peak_LED+dleddt];
      y_g[0] = trace[1][index_peak_LED+dleddt];
      canv->cd(1);
      graphPeak_LED = new TGraphErrors(1,x_g,y_g,0,0);
      graphPeak_LED->SetMarkerStyle(20);
      graphPeak_LED->SetMarkerColor(kBlue);
      graphPeak_LED->Draw("psame");

      // On trace_DLED
      float x_g_DLED[1], y_g_DLED[1];
      x_g_DLED[0] = trace_DLED[0][index_peak_LED];
      y_g_DLED[0] = trace_DLED[1][index_peak_LED];
      canv->cd(2);
      graphPeak_LED_DLED = new TGraphErrors(1,x_g_DLED,y_g_DLED,0,0);
      graphPeak_LED_DLED->SetMarkerStyle(20);
      graphPeak_LED_DLED->SetMarkerColor(kBlue);
      graphPeak_LED_DLED->Draw("psame");
    }

    canv->Update();

    // let the user to interact with the canvas, like zooming the axis
    // show next trace by pressing any key

    if(!running_graph and delete_bool){
      while(!gSystem->ProcessEvents()) {
         if (canv->WaitPrimitive()==0)
         {
            break;
         }
      }
    }

    if(delete_bool) {
        delete graph1; delete graph2;
        if(display_peaks_now){delete graphPeaks; delete graphPeaks_DLED;}
        if(show_peak_LED_bool){delete graphPeak_LED;}
    }
}


//------------------------------------------------------------------------------
void show_AVG_trace_window(TCanvas *canv, float *tracet, float *tracev, int trace_length, bool delete_bool){

	if(reverse_bool){
          for(int xx=0; xx<trace_length; xx++){
            tracev[xx] = -tracev[xx];
          }
        }

	canv->cd();
	TGraph *trace = new TGraph(trace_length, tracet, tracev);
	trace->Draw("AL");
	//if(delete_bool) delete trace;

}


//------------------------------------------------------------------------------
void Get_DCR_temp_and_errDCR_temp(){
    DCR_temp[nfile] = expDel->GetParameter(0)*TMath::Power(10,9);
    errDCR_temp[nfile] = expDel->GetParError(0)*TMath::Power(10,9);
}

//-----------------------------------------------------------------------------
double GetDCRfromDelays(){
    return expDel->GetParameter(0)*TMath::Power(10,9);
}

//-----------------------------------------------------------------------------
double GetErrDCRfromDelays(){
    return expDel->GetParError(0)*TMath::Power(10,9);
}

//------------------------------------------------------------------------------
void find_DCR_0_5_pe_and_1_5_pe_manual(){

    double delta = (pe_1_5_vect[nfile] - pe_0_5_vect[nfile])/n_DCR;
    int index_0_5=0, index_1_5=0;

    for(int i=0; i<n_DCR; i++){
        if((thr_to_find_peaks_vect[nfile][i] > pe_0_5_vect[nfile] - delta) and (thr_to_find_peaks_vect[nfile][i] < pe_0_5_vect[nfile] + delta))
            index_0_5 = i;
        if((thr_to_find_peaks_vect[nfile][i] > pe_1_5_vect[nfile] - delta) and (thr_to_find_peaks_vect[nfile][i] < pe_1_5_vect[nfile] + delta))
            index_1_5 = i;
    }

    //0.5pe
    DCR_pe_0_5_vect[nfile] = DCR[nfile][index_0_5];
    errDCR_pe_0_5_vect[nfile] = errDCR[nfile][index_0_5];

    //1.5pe
    DCR_pe_1_5_vect[nfile] = DCR[nfile][index_1_5];
    errDCR_pe_1_5_vect[nfile] = errDCR[nfile][index_1_5];
}


//------------------------------------------------------------------------------
void find_DCR_0_5_pe_and_1_5_pe_auto(){
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
    expDel-> SetParName(0,"DCR");
    expDel-> SetParName(1,"MultCost");
    expDel->SetParameter(0,0.001);
    ptrHistDelays[nfile] -> Fit(expDel, "", "", expDelLow, expDelHigh);
    expDel->Draw("same");
}

//------------------------------------------------------------------------------
void fit_hist_peaks(TCanvas *canv, TH1D *hist){

    float mean1, mean2;
    float errmean1, errmean2;

    float fit1Low, fit1High, fit2Low, fit2High;

    cout<<"--------------------"<<endl;
    cout<<"Please enter:"<<endl;
    cout<<"fit1Low  "; if (scanf("%f", &fit1Low) == 1){;}
    cout<<"fit1High "; if (scanf("%f", &fit1High) == 1){;}
    cout<<"fit2Low  "; if (scanf("%f", &fit2Low) == 1){;}
    cout<<"fit2High "; if (scanf("%f", &fit2High) == 1){;}

    canv->cd();

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
void fit_hist_peaks_0pe_1pe_2pe(TCanvas *canv, TH1D *hist){

  float Mean_peak_0, Mean_peak_1, Mean_peak_2;
  float errMean_peak_0, errMean_peak_1, errMean_peak_2;

  float Sigma_peak_0, Sigma_peak_1, Sigma_peak_2;
  float errSigma_peak_0, errSigma_peak_1, errSigma_peak_2;

  float H_peak_0, H_peak_1, H_peak_2;
  float errH_peak_0, errH_peak_1, errH_peak_2;

  float Mean_hg, errMean_hg, Std_hg, errStd_hg;

  int Entries;

  // find range values to find fit fit ranges
  bool range1_low_low_found = false;
  bool range1_low_high_found = false;
  bool range2_low_low_found = false;
  bool range2_low_high_found = false;
  bool range2_high_low_found = false;
  bool range2_high_high_found = false;


  for(i=0; i<histLED_nbins; i++){
    if((hist->GetBinCenter(i) > range1_low_low_mV) && !range1_low_low_found){
      range1_low_low_bin = i;
      range1_low_low_found = true;
    }
    if((hist->GetBinCenter(i) > range1_low_high_mV) && !range1_low_high_found){
      range1_low_high_bin = i;
      range1_low_high_found = true;
    }
    if((hist->GetBinCenter(i) > range2_low_low_mV) && !range2_low_low_found){
      range2_low_low_bin = i;
      range2_low_low_found = true;
    }
    if((hist->GetBinCenter(i) > range2_low_high_mV) && !range2_low_high_found){
      range2_low_high_bin = i;
      range2_low_high_found = true;
    }
    if((hist->GetBinCenter(i) > range2_high_low_mV) && !range2_high_low_found){
      range2_high_low_bin = i;
      range2_high_low_found = true;
    }
    if((hist->GetBinCenter(i) > range2_high_high_mV) && !range2_high_high_found){
      range2_high_high_bin = i;
      range2_high_high_found = true;
    }
  }
  cout<<"0 high, 1 low\t"<<ptrHistLED->GetBinCenter(range1_low_low_bin)<<endl;
  cout<<"0 high, 1 low\t"<<ptrHistLED->GetBinCenter(range1_low_high_bin)<<endl;
  cout<<"1 high, 2 low\t"<<ptrHistLED->GetBinCenter(range2_low_low_bin)<<endl;
  cout<<"1 high, 2 low\t"<<ptrHistLED->GetBinCenter(range2_low_high_bin)<<endl;
  cout<<"2 high       \t"<<ptrHistLED->GetBinCenter(range2_high_low_bin)<<endl;
  cout<<"2 high       \t"<<ptrHistLED->GetBinCenter(range2_high_high_bin)<<endl;

  // find fit ranges
  fit0Low   = -20;
  fit0High  = 0;
  fit1Low   = 0;
  fit1High  = 0;
  fit2Low   = 0;
  fit2High  = 0;

  int min;

  // find fit0High and fit1Low
  min = 100000;
  for(int i=range1_low_low_bin; i<range1_low_high_bin; i++){
    if(hist->GetBinContent(i) < min){
      min = hist->GetBinContent(i);
      fit0High = (int)(hist->GetBinCenter(i));
    }
  }
  fit1Low = fit0High;

  // find fit1High and fit2Low
  min = 100000;
  for(int i=range2_low_low_bin; i<range2_low_high_bin; i++){
    if(hist->GetBinContent(i) < min){
      min = hist->GetBinContent(i);
      fit1High = (int)(hist->GetBinCenter(i));
    }
  }
  fit2Low = fit1High;

  // find fit2High
  min = 100000;
  for(int i=range2_high_low_bin; i<range2_high_high_bin; i++){
    if(hist->GetBinContent(i) < min){
      min = hist->GetBinContent(i);
      fit2High = (int)(hist->GetBinCenter(i));
    }
  }

  cout<<endl;
  cout<<"fit0High\t"<<fit0High<<endl;
  cout<<"fit1Low\t\t"<<fit1Low<<endl;
  cout<<"fit1High\t"<<fit1High<<endl;
  cout<<"fit2Low\t\t"<<fit2Low<<endl;
  cout<<"fit2High\t"<<fit2High<<endl;
  cout<<endl;

  canv->cd();

  //-------------
  //---[ 0pe ]---
  //-------------
  hist -> Fit("gausFit0", "+", "", fit0Low, fit0High); //gaus fit of the 1st peak
  H_peak_0   = gausFit0->GetParameter(0);
  errH_peak_0= gausFit0->GetParError(0);
  Mean_peak_0     = gausFit0->GetParameter(1);
  errMean_peak_0  = gausFit0->GetParError(1);
  Sigma_peak_0     = gausFit0->GetParameter(2);
  errSigma_peak_0  = gausFit0->GetParError(2);
  gausFit0->Draw("same");

  //-------------
  //---[ 1pe ]---
  //-------------
  hist -> Fit("gausFit1", "+", "", fit1Low, fit1High); //gaus fit of the 1st peak
  H_peak_1   = gausFit1->GetParameter(0);
  errH_peak_1= gausFit1->GetParError(0);
  Mean_peak_1     = gausFit1->GetParameter(1);
  errMean_peak_1  = gausFit1->GetParError(1);
  Sigma_peak_1     = gausFit1->GetParameter(2);
  errSigma_peak_1  = gausFit1->GetParError(2);
  gausFit1->Draw("same");

  //-------------
  //---[ 2pe ]---
  //-------------
  hist -> Fit("gausFit2", "+", "", fit2Low, fit2High); //gaus fit of the 2st peak
  H_peak_2    = gausFit2->GetParameter(0);
  errH_peak_2 = gausFit2->GetParError(0);
  Mean_peak_2      = gausFit2->GetParameter(1);
  errMean_peak_2   = gausFit2->GetParError(1);
  Sigma_peak_2     = gausFit2->GetParameter(2);
  errSigma_peak_2  = gausFit2->GetParError(2);
  gausFit2->Draw("same");

  gain = Mean_peak_2 - Mean_peak_1; //gain is the difference between 1pe peak and 2pe peak
  errgain = TMath::Sqrt( errMean_peak_1*errMean_peak_1 + errMean_peak_2*errMean_peak_2 ); //error propagation

  //---------------------
  //---[ GLOBAL HIST ]---
  //---------------------
  Entries = hist->GetEntries();
  Mean_hg = hist->GetMean();
  errMean_hg = hist->GetMeanError();
  Std_hg  = hist->GetStdDev();
  errStd_hg = hist->GetStdDevError();


  // CROSS TALK from LED
  double Area0 = H_peak_0*Sigma_peak_0*TMath::Power(2*TMath::Pi(),0.5)/1;
  double Prob_0pe = Area0/Entries;
  double Mu = -TMath::Log(Prob_0pe);
  double Prob_1pe = Mu*TMath::Exp(-Mu);
  double Area1 = H_peak_1*Sigma_peak_1*TMath::Power(2*TMath::Pi(),0.5)/1;
  double Prob_1peS = Area1/Entries;
  double Prob_Cross_Talk = 1-(Prob_1peS/Prob_1pe);



  //-----------------
  //---[ RESULTS ]---
  //-----------------
  int index = 0;
  cout<<endl;
  // cout<<"Enter vector index:"<<endl;
  // cin>>index;

  cout<<endl;
  cout<<"// Window for LED peak: "<<"("<<minLED_amp<<", "<<maxLED_amp<<") ns"<<"  (minLED_amp, maxLED_amp)"<<endl;

  // peak 0
  cout<<"H_peak_0["<<index<<"]\t\t= "<<H_peak_0<<";"<<endl;
  cout<<"errH_peak_0["<<index<<"]\t\t= "<<errH_peak_0<<";"<<endl;
  cout<<"Sigma_peak_0["<<index<<"]\t\t= "<<Sigma_peak_0<<";"<<endl;
  cout<<"errSigma_peak_0["<<index<<"]\t= "<<errSigma_peak_0<<";"<<endl;

  // peak 1
  cout<<"H_peak_1["<<index<<"]\t\t= "<<H_peak_1<<";"<<endl;
  cout<<"errH_peak_1["<<index<<"]\t\t= "<<errH_peak_1<<";"<<endl;
  cout<<"Sigma_peak_1["<<index<<"]\t\t= "<<Sigma_peak_1<<";"<<endl;
  cout<<"errSigma_peak_1["<<index<<"]\t= "<<errSigma_peak_1<<";"<<endl;
  cout<<"Mean_peak_1["<<index<<"]\t\t= "<<Mean_peak_1<<";"<<endl;
  cout<<"errMean_peak_1["<<index<<"]\t= "<<errMean_peak_1<<";"<<endl;

  // peak 2
  cout<<"Mean_peak_2["<<index<<"]\t\t= "<<Mean_peak_2<<";"<<endl;
  cout<<"errMean_peak_2["<<index<<"]\t= "<<errMean_peak_2<<";"<<endl;

  // GLOBAL HIST
  cout<<"Mean_hg["<<index+2<<"]\t\t= "<<Mean_hg<<";"<<endl;
  cout<<"errMean_hg["<<index+2<<"]\t\t= "<<errMean_hg<<";"<<endl;
  cout<<"Std_hg["<<index+2<<"]\t\t= "<<Std_hg<<";"<<endl;
  cout<<"errStd_hg["<<index+2<<"]\t\t= "<<errStd_hg<<";"<<endl;
  cout<<"Entries["<<index<<"]\t\t= "<<Entries<<";"<<endl;


  // Cross Talk
  cout<<endl;
  cout << "Cross Talk\t\t= " << Prob_Cross_Talk*100 <<" \%" << endl;

}



//------------------------------------------------------------------------------
void fit_hist_peaks_gaus_sum_012(TCanvas *canv, TH1D *hist, bool evaluate_cross_talk){

  float Mean_peak_0, Mean_peak_1, Mean_peak_2;
  float errMean_peak_0, errMean_peak_1, errMean_peak_2;

  float Sigma_peak_0, Sigma_peak_1, Sigma_peak_2;
  float errSigma_peak_0, errSigma_peak_1, errSigma_peak_2;

  float H_peak_0, H_peak_1, H_peak_2;
  float errH_peak_0, errH_peak_1, errH_peak_2;

  float Gain, errGain, Sigma_add, errSigma_add, Diff_Gain_0_1, errDiff_Gain_0_1;

  float Mean_hg, errMean_hg, Std_hg, errStd_hg;

  int Entries;


  canv->cd();

  // Set Name
  gaus_sum_012->SetParName(0, "H0");
  gaus_sum_012->SetParName(1, "H1");
  gaus_sum_012->SetParName(2, "H2");
  gaus_sum_012->SetParName(3, "H3");
  gaus_sum_012->SetParName(4, "gain");
  gaus_sum_012->SetParName(5, "s0");
  gaus_sum_012->SetParName(6, "sadd");
  gaus_sum_012->SetParName(7, "V0");
  gaus_sum_012->SetParName(8, "dg01");

  // Initialize Parameters
  gaus_sum_012->SetParameter(0, 800);  // H0
  gaus_sum_012->SetParameter(1, 800);  // H1
  gaus_sum_012->SetParameter(2, 800);  // H2
  gaus_sum_012->SetParameter(3, 800);  // H3
  gaus_sum_012->SetParameter(4, 15);   // g
  gaus_sum_012->SetParameter(5, 2);    // s0
  gaus_sum_012->SetParameter(6, 1);    // sadd
  gaus_sum_012->SetParameter(7, 1);    // V0
  gaus_sum_012->SetParameter(8, 0.1);  // dg01

  // Fit Range
  float fit_low  = -10;
  float fit_high = 55;


  // FIT
  gaus_sum_012->SetRange(fit_low, fit_high);
  hist -> Fit("gaus_sum_012", "+", "", fit_low, fit_high);

  // gain
  Gain              = gaus_sum_012->GetParameter(4);
  errGain           = gaus_sum_012->GetParError(4);
  Diff_Gain_0_1     = gaus_sum_012->GetParameter(8);
  errDiff_Gain_0_1  = gaus_sum_012->GetParError(8);

  // sigma
  Sigma_add        = gaus_sum_012->GetParameter(6);
  errSigma_add     = gaus_sum_012->GetParError(6);


  //--------------
  //---[ 0pe ]---
  //-------------
  H_peak_0          = gaus_sum_012->GetParameter(0);
  errH_peak_0       = gaus_sum_012->GetParError(0);
  Mean_peak_0       = gaus_sum_012->GetParameter(7);
  errMean_peak_0    = gaus_sum_012->GetParError(7);
  Sigma_peak_0      = gaus_sum_012->GetParameter(5);
  errSigma_peak_0   = gaus_sum_012->GetParError(5);


  //--------------
  //---[ 1pe ]---
  //-------------
  H_peak_1          = gaus_sum_012->GetParameter(1);
  errH_peak_1       = gaus_sum_012->GetParError(1);
  Mean_peak_1       = Mean_peak_0 + Gain + Diff_Gain_0_1;
  errMean_peak_1    = TMath::Sqrt( errMean_peak_0*errMean_peak_0 + errGain*errGain + errDiff_Gain_0_1*errDiff_Gain_0_1);
  Sigma_peak_1      = TMath::Sqrt( Sigma_peak_0*Sigma_peak_0 + Sigma_add*Sigma_add );
  errSigma_peak_1   = TMath::Sqrt( (Sigma_peak_0*Sigma_peak_0*errSigma_peak_0*errSigma_peak_0 + Sigma_add*Sigma_add*errSigma_add*errSigma_add) / ( Sigma_peak_0*Sigma_peak_0 + Sigma_add*Sigma_add ) );


  //--------------
  //---[ 2pe ]---
  //-------------
  H_peak_2          = gaus_sum_012->GetParameter(2);
  errH_peak_2       = gaus_sum_012->GetParError(2);
  Mean_peak_2       = Mean_peak_1 + Gain;
  errMean_peak_2    = TMath::Sqrt( errMean_peak_1*errMean_peak_1 + errGain*errGain );
  Sigma_peak_2      = TMath::Sqrt( Sigma_peak_0*Sigma_peak_0 + 4*Sigma_add*Sigma_add );
  errSigma_peak_2   = TMath::Sqrt( (Sigma_peak_0*Sigma_peak_0*errSigma_peak_0*errSigma_peak_0 + 4*4*Sigma_add*Sigma_add*errSigma_add*errSigma_add) / ( Sigma_peak_0*Sigma_peak_0 + 4*Sigma_add*Sigma_add ) );


  //gaus_sum_012->Draw("same");

  //---------------------
  //---[ GLOBAL HIST ]---
  //---------------------
  Entries = hist->GetEntries();
  Mean_hg = hist->GetMean();
  errMean_hg = hist->GetMeanError();
  Std_hg  = hist->GetStdDev();
  errStd_hg = hist->GetStdDevError();

  //-----------------------------
  //---[ CROSS TALK from LED ]---
  //-----------------------------
  double Area0 = H_peak_0*Sigma_peak_0*TMath::Power(2*TMath::Pi(),0.5)/1;
  double Prob_0pe = Area0/Entries;
  double Mu = -TMath::Log(Prob_0pe);
  double Prob_1pe = Mu*TMath::Exp(-Mu);
  double Area1 = H_peak_1*Sigma_peak_1*TMath::Power(2*TMath::Pi(),0.5)/1;
  double Prob_1peS = Area1/Entries;
  double Prob_Cross_Talk = 1-(Prob_1peS/Prob_1pe);


  //-----------------
  //---[ RESULTS ]---
  //-----------------
  int index = 0;
  cout<<endl;
  // cout<<"Enter vector index:"<<endl;
  // cin>>index;

  cout<<endl;
  cout<<"// Window for LED peak: "<<"("<<minLED_amp<<", "<<maxLED_amp<<") ns"<<"  (minLED_amp, maxLED_amp)"<<endl;

  // peak 0
  cout<<"H_peak_0["<<index<<"]          = "<<H_peak_0<<";"<<endl;
  cout<<"errH_peak_0["<<index<<"]       = "<<errH_peak_0<<";"<<endl;
  cout<<"Sigma_peak_0["<<index<<"]      = "<<Sigma_peak_0<<";"<<endl;
  cout<<"errSigma_peak_0["<<index<<"]   = "<<errSigma_peak_0<<";"<<endl;

  // peak 1
  cout<<"H_peak_1["<<index<<"]          = "<<H_peak_1<<";"<<endl;
  cout<<"errH_peak_1["<<index<<"]       = "<<errH_peak_1<<";"<<endl;
  cout<<"Sigma_peak_1["<<index<<"]      = "<<Sigma_peak_1<<";"<<endl;
  cout<<"errSigma_peak_1["<<index<<"]   = "<<errSigma_peak_1<<";"<<endl;
  cout<<"Mean_peak_1["<<index<<"]       = "<<Mean_peak_1<<";"<<endl;
  cout<<"errMean_peak_1["<<index<<"]    = "<<errMean_peak_1<<";"<<endl;

  // peak 2
  cout<<"Mean_peak_2["<<index<<"]       = "<<Mean_peak_2<<";"<<endl;
  cout<<"errMean_peak_2["<<index<<"]    = "<<errMean_peak_2<<";"<<endl;

  // GAIN
  cout<<"GAIN["<<index<<"]              = "<<Gain<<";"<<endl;
  cout<<"errGAIN["<<index<<"]           = "<<errGain<<";"<<endl;

  // GLOBAL HIST
  cout<<"Mean_hg["<<index+2<<"]           = "<<Mean_hg<<";"<<endl;
  cout<<"errMean_hg["<<index+2<<"]        = "<<errMean_hg<<";"<<endl;
  cout<<"Std_hg["<<index+2<<"]            = "<<Std_hg<<";"<<endl;
  cout<<"errStd_hg["<<index+2<<"]         = "<<errStd_hg<<";"<<endl;
  cout<<"Entries["<<index<<"]           = "<<Entries<<";"<<endl;

  // Info
  cout<<endl;
  cout<<"// Fit based on the function reported in: "<<endl;
  cout<<"//       Mallamaci Manuela, Mariotti Mosé - Report SiPM tests"<<endl;

  // Cross Talk
  if(evaluate_cross_talk){
    cout<<endl;
    cout << "Cross Talk           = " << Prob_Cross_Talk*100 <<" \%" << endl;
  }




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
// void find_offset_mod(){  // USING TSpectrum
//
//   float bin_width = (trace[0][1023] - trace[0][0])/1024;
//   float x_low = trace[0][0] - bin_width/2.;
//   float x_up = trace[0][1023] + bin_width/2.;
//   float n_bin = x_up - x_low;
//
//   // trace hist
//   TH1F *trace_histo = new TH1F("histTrace","",n_bin,x_low,x_up);
//
//   // filling trace hist:
//   for(ii=x_low; ii<x_up; ii++){
//         trace_histo->Fill(trace[0][ii], trace[1][ii]);
//   }
//
//   TSpectrum *s = new TSpectrum();
//
//   //find BACKGROUND:
//   TH1 *bckg = s->Background(trace_histo,100);
//
//   //subtract BACKGROUND:
//   trace_histo->Add(bckg,trace_histo,-1,1);
//
//   // update trace
//   for(ii=x_low; ii<x_up; ii++){
//         trace[0][ii] = trace_histo->GetBinCenter(ii) - bin_width/2.;
//         trace[0][ii] = trace_histo->GetBinCenter(ii) - bin_width/2.;
//         trace[1][ii] = trace_histo->GetBinContent(ii);
//         //ptrAllTrace->Fill(trace[1][ii]);
//   }
//
//   delete trace_histo;
//   delete s;
//   delete bckg;
//
// }

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
double find_area_trace(double center, double low, double high){
    double area = 0;
    double center_index = 0, low_index = 0, high_index = 0;

    bool found_center = false;
    bool found_low = false;
    bool found_high = false;

    // find center
    for(int i=0; i<trace_length; i++){ // loop on the trace (NOT DLED)

      // the first time that time > low
      if((trace[0][i]>low) && !found_low){
          low_index = i; // index corresponding to low
          found_low = true;
      }

      // the first time that time > center
      if((trace[0][i]>center) && !found_center){
          center_index = i; // index corresponding to center
          found_center = true;
      }

      // the first time that time > high
      if((trace[0][i]>high) && !found_high){
          high_index = i; // index corresponding to high
          found_high = true;
      }
    }

    // sum all the voltage values times delta_temp
    for(int i=low_index; i<high_index; i++){
      area += trace[1][i]*(trace[0][i+1]-trace[0][i]); // area of the small rectangle (amplitude x delta_temp)
    }

    return area;
}

//------------------------------------------------------------------------------
void DLED_offset_remove(){
    int index_peak_trace;
    for(int i = 0; i<num_peaks; i++){ //loop over peaks
        index_peak_trace = index_vect[i] + dleddt;

        for(int j=index_vect[i]; j<trace_DLED_length; j++){ //loop on trace_DLED

            trace_DLED[1][j] += trace_DLED[1][index_peak_trace] * ( TMath::Exp( ( j - dleddt - index_peak_trace) / tau ) - TMath::Exp( ( j - index_peak_trace) / tau ) );

        }
    }
}

//------------------------------------------------------------------------------
void SmoothTraceN(int n){
    int max_k = n - n%2;
    double den = n + 1 - n%2;

    // cout<<"SMOOTHING TRACE "<<n<<endl;

    // cout<<max_k<<" "<<den<<endl;
    // getchar();

    for(int i=0; i<trace_length; i+=n){

        // sum on n
        for(int j=1; j<n; j++){
          trace[1][i] += trace[1][i+j];
        }

        // divide for n
        trace[1][i]   /= n;
        for(int j=1; j<n; j++){
          trace[1][i+j] = trace[1][i];
        }

        if(i!=0){
            int diff = -(int)n/2;
            double gap = (double)trace[1][i] - (double)trace[1][i-1];
            double mult = 0;
            // cout<<diff<<endl;
            for(int k=0; k<max_k; k++){
                if(diff<0){
                    mult+=1;
                    // cout<<gap<<" "<<den<<" "<<mult<<" "<<gap/den*mult<<endl;
                    trace[1][i+diff] += gap/den*mult;
                }else{
                    trace[1][i+diff] -= gap/den*mult;
                    mult-=1;
                }
                diff+=1;
            }
            // cout<<diff<<" "<<den<<" "<<mult<<endl;
            // cout<<i<<endl;
        }
    }
    // getchar();
}

//------------------------------------------------------------------------------
Float_t GetMean(std::vector<Float_t> vec)
{
   // Compute mean of a vector.
   // From: https://stackoverflow.com/questions/7616511/calculate-mean-and-standard-deviation-from-a-vector-of-samples-in-c-using-boos/12405793#12405793

   Float_t sum = std::accumulate(vec.begin(), vec.end(), 0.0);

   return (sum / vec.size());
}


//-----------------------------------------------------------------------------
Float_t GetStdDev(std::vector<Float_t> vec)
{
   // Compute standard deviation of a vector.
   // From: https://stackoverflow.com/questions/7616511/calculate-mean-and-standard-deviation-from-a-vector-of-samples-in-c-using-boos/12405793#12405793

   Float_t accum = 0.0;
   Float_t m = GetMean(vec);
   std::for_each (vec.begin(), vec.end(), [&](const Float_t d) {
    accum += (d - m) * (d - m);
   });

   return (sqrt(accum / (vec.size()-1)));
}

//-----------------------------------------------------------------------------
void EvaluateCrossTalk(double DCR_0_5, double errDCR_0_5, double DCR_1_5, double errDCR_1_5, double* CT){
    // CT[0] = CrossTalk, CT[1] = errCrossTalk

    // CrossTalk = (DCR @ 1.5pe) / (DCR @ 0.5pe)
    CT[0] = DCR_1_5/DCR_0_5;
    CT[1] = CT[0] * TMath::Sqrt( (errDCR_0_5/DCR_0_5)*(errDCR_0_5/DCR_0_5) + (errDCR_1_5/DCR_1_5)*(errDCR_1_5/DCR_1_5) );
}


//------------------------------------------------------------------------------
//------------------------------[   READ FILES   ]------------------------------
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
void Read_Agilent_CAEN(string file, int last_event_n, bool display, TCanvas *c){

//     ifstream OpenFile (file.c_str());
//     reading = true;
//
//     n_ev=0;
//
// //***** READ FILE
//     while(!OpenFile.eof() and (reading)){
//
//         if(n_ev%1000==0)
//             cout<<"Read ev\t"<<n_ev<<endl;
//
// //***** READ HEADER FROM FILE
//         //select device
//         if(Agilent_MSO6054A and not Digitizer_CAEN){
//             OpenFile>>temp;
//             OpenFile>>temp;
//             trace_length = atoi(temp);
//             OpenFile>>temp>>temp;
//         }
//         else{
//             if(Digitizer_CAEN and not Agilent_MSO6054A){
//                 OpenFile>>temp>>temp;
//                 OpenFile>>temp;
//                 trace_length = atoi(temp);
//                 if(all_events_same_window){
//                     trace_length = atoi(temp)*last_event_n;
//                     one_window = atoi(temp);
//                 }
//                 OpenFile>>temp>>temp>>temp>>temp>>temp>>temp>>temp>>temp>>temp>>temp>>temp>>temp>>temp>>temp>>temp>>temp>>temp;
//         }else{
//             cout<<"ERROR: check acquisition device"<<endl;
//         }
//         }
//         if(n_ev==0){
//             cout<<"Acquisition device: ";
//             if(Agilent_MSO6054A) cout<<"Agilent_MSO6054A"<<endl;
//             if(Digitizer_CAEN)   cout<<"Digitizer_CAEN"<<endl;
//             cout<<"Number points "<<trace_length<<endl;
//         }
//
// //***** CREATE TRACE
//         //trace
//         trace = new float*[2];
//         for(i = 0; i < 2; i++) {
//             trace[i] = new float[trace_length];
//         }
//         //trace_AVG
//         if(average==true and n_ev==0){
//             trace_AVG = new float*[2];
//             for(i = 0; i < 2; i++) {
//                 trace_AVG[i] = new float[trace_length];
//             }
//         }
//
//
// //***** READ TRACE FROM FILE
//         if(Agilent_MSO6054A){
//                 for(i=0; i<trace_length; i++){
//                     OpenFile>>temp;
//                     trace[0][i] = atof(temp);
//                     OpenFile>>temp;
//                     trace[1][i]  = -atof(temp); //AdvanSid is an inverting amplifier, so I have negative signals. In order to analyze them, I reverse the signals. For the plot I can re-reverse them
//         }
//         }
//         else{
//             if(Digitizer_CAEN){
//                 for(i=0; i<trace_length; i++){
//                     trace[0][i] = i; //1point=1ns
//                     OpenFile>>temp;
//                     trace[1][i]  = -atof(temp)/1024 * 1000; //1024 channels from 0 V to 1 V, expressed in mV
//
//                     if(i==(one_window-1)){
//                         OpenFile>>temp>>temp;
//                         OpenFile>>temp;
//                         OpenFile>>temp>>temp>>temp>>temp>>temp>>temp>>temp>>temp>>temp>>temp>>temp>>temp>>temp>>temp>>temp>>temp>>temp;
//                     }
//
//                 }
//         }
//         }
//
// //***** DLED
//         if(DLED_bool){
//             DLED(trace_length,dleddt);
//             //Now: trace_DLED
//         }else{
//             for(ii=0; ii<trace_DLED_length; ii++){
//                 trace_DLED[0][ii] = trace[0][ii];
//                 trace_DLED[1][ii] = trace[1][ii];
//             }
//         }
//
// //***** AVERAGE
//         if(average){
//             if(n_ev==0){
//                 for(i=0; i< trace_length; i++){
//                     trace_AVG[0][i]=trace[0][i];
//                     trace_AVG[1][i]=0;
//                 }
//             }
//         }
//
//         if(display){
//             if(!display_one_ev) display_peaks_now = display_peaks;
//             if(display_one_ev and (n_ev==ev_to_display)) display_peaks_now = display_peaks;
//         }
//
//         if(find_peak_in_a_selected_window){
//
//             float *peak = new float[2];
//             index_for_peak = find_peak_fix_time(mintp, maxtp);
//             peak[0] = trace_DLED[0][index_for_peak];
//             peak[1] = trace_DLED[1][index_for_peak];
//             ptrHistLED->Fill(peak[1]);
//                 delete[] peak;
//             }
//
// //***** DCR
//         if(DCR_DELAYS_bool){
//            find_peaks(thr_to_find_peaks,max_peak_width, min_peak_width,blind_gap,DCR_DELAYS_bool);
//         }
//         else{
//             if(drawHistAllPeaks){//all peaks but not DCR
//                 thr_to_find_peaks = 10; //mV
//                 find_peaks(thr_to_find_peaks,max_peak_width, min_peak_width,blind_gap,DCR_DELAYS_bool);
//             }
//         }
//
//
//
// //***** DISPLAY
//         if(display){
//             miny1 = -10; maxy1 = 100; miny2 = -10; maxy2 = 100;
//             if(!display_one_ev){
//                 if(n_ev==0){
//                     c->Divide(1,2);
//                     c->SetGrid();
//                 }
//                 show_trace2(c, trace[0], trace[1], trace_DLED[0], trace_DLED[1], trace_length, trace_DLED_length, miny1, maxy1, miny2, maxy2, line_bool, true);
//             }else{
//                 if(n_ev==ev_to_display){
//                     c->Divide(1,2);
//                     c->SetGrid();
//                     show_trace2(c, trace[0], trace[1], trace_DLED[0], trace_DLED[1], trace_length, trace_DLED_length, miny1, maxy1, miny2, maxy2, line_bool, false);
//                 }
//             }
//
//         }
//
//
//         if(n_ev==last_event_n-1)
//             reading=false;
//
//         delete []trace[0];
//         delete []trace[1];
//         delete []trace_DLED[0];
//         delete []trace_DLED[1];
//         delete []peak;
//         n_ev++;
//     }//file is closed
//
//     n_ev_tot = n_ev;
//     cout<<"Last event "<<n_ev_tot<<endl;
}

//------------------------------------------------------------------------------
void ReadBin(string filename, int last_event_n, bool display, TCanvas *c){
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
   int count_peak_window = 0;

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
       return;
   }

   // read file header
   if (fread(&fh, sizeof(fh), 1, f) != 1){
    cout << "Problem in reading binary file." << endl;
    fclose(f);
    return;
   }
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
   if (fread(&th, sizeof(th), 1, f) != 1){
    cout << "Problem in reading binary file." << endl;
    fclose(f);
    return;
   }
   if (memcmp(th.time_header, "TIME", 4) != 0) {
      //printf("Invalid time header in file \'%s\', aborting.\n", filename);
//       return 0;
      cout<<"Invalid time header in file"<<endl;
   }

   for (b = 0 ; ; b++) {
      // read board header
      if (fread(&bh, sizeof(bh), 1, f) != 1){
        cout << "Problem in reading binary file." << endl;
        fclose(f);
        return;
      }
      if (memcmp(bh.bn, "B#", 2) != 0) {
         // probably event header found
         fseek(f, -4, SEEK_CUR);
         break;
      }

      //printf("Found data for board #%d\n", bh.board_serial_number);

      // read time bin widths
      memset(bin_width[b], sizeof(bin_width[0]), 0);
      for (chn=0 ; chn<5 ; chn++) {
         if (fread(&ch, sizeof(ch), 1, f) != 1){
          cout << "Problem in reading binary file." << endl;
          fclose(f);
          return;
         }
         if (ch.c[0] != 'C') {
            // event header found
            fseek(f, -4, SEEK_CUR);
            break;
         }
         i = ch.cn[2] - '0' - 1;
         //printf("Found timing calibration for channel #%d\n", i+1);
         if (fread(&bin_width[b][i][0], sizeof(float), 1024, f) != 1024){
          cout << "Problem in reading binary file." << endl;
          fclose(f);
          return;
          }
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

       n_ev = n;

       if(n_ev==last_event_n-1)
              break;

      // read event header
      i = (int)fread(&eh, sizeof(eh), 1, f);
      if (i < 1)
         break;



      if(n_ev%1000==0)
            cout<<"Read ev\t"<<n_ev<<endl;


      //printf("Found event #%d %d %d\n", eh.event_serial_number, eh.second, eh.millisecond);

      // loop over all boards in data file
      for (b=0 ; b<n_boards ; b++) {

         // read board header
         if (fread(&bh, sizeof(bh), 1, f) != 1){
          cout << "Problem in reading binary file." << endl;
          fclose(f);
          return;
         }
         if (memcmp(bh.bn, "B#", 2) != 0) {
            //printf("Invalid board header in file \'%s\', aborting.\n", filename);
//             return 0;
            cout<<"Invalid board header in file"<<endl;
         }

         // read trigger cell
         if (fread(&tch, sizeof(tch), 1, f) != 1){
          cout << "Problem in reading binary file." << endl;
          fclose(f);
          return;
         }
         if (memcmp(tch.tc, "T#", 2) != 0) {
            //printf("Invalid trigger cell header in file \'%s\', aborting.\n", filename);
            cout<<"Invalid trigger cell header in file"<<endl;
//             return 0;
         }



         // reach channel data
         for (chn=0 ; chn<4 ; chn++) {
            // read channel header
            if (fread(&ch, sizeof(ch), 1, f) != 1){
              cout << "Problem in reading binary file." << endl;
              fclose(f);
              return;
            }
            if (ch.c[0] != 'C') {
               // event header found
               fseek(f, -4, SEEK_CUR);
               break;
            }
            chn_index = ch.cn[2] - '0' - 1;
            if (fread(&scaler, sizeof(int), 1, f) != 1){
              cout << "Problem in reading binary file." << endl;
              fclose(f);
              return;
            }

            if (fread(voltage, sizeof(short), 1024, f) != 1024){
              cout << "Problem in reading binary file." << endl;
              fclose(f);
              return;
            }

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

        peak_LED = new float[2];

        // cout<<endl<<"trace[0][k]"<<endl;
        for(int k=0; k<trace_length; k++){
             trace[0][k] = time[0][0][k];
             if(reverse_bool) trace[1][k] = -waveform[0][0][k]*1000; //to convert in mV
             else trace[1][k] = waveform[0][0][k]*1000;

             // if(k<10)cout<<trace[0][k]<<"\t";
        }
        // cout<<endl;

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

      if(smooth_trace_bool){
          int n_SmootTraceN = 2;
          SmoothTraceN(n_SmootTraceN);
          cout<<"Smooth trace "<<n_SmootTraceN<<" points"<<endl;

        // smooth_trace_step();
        // smooth_trace_3();
        // smooth_trace_4();
        // smooth_trace_5();
        // SmoothTraceMarkov(4);

        // //////////////////////////
        // ///   D_E_B_U_G
        // //////////////////////////
        // float temp[1024], temp1[1024];
        // for(int k=0; k<1024; k++) temp[k] = trace[1][k];  // dummmy, the original trace
        // SmoothTraceN(3);
        // for(int k=0; k<1024; k++) temp1[k] = trace[1][k]; // temp1 is trace from smooth_trace_4
        // for(int k=0; k<1024; k++) trace[1][k] = temp[k];  // reset trace
        // smooth_trace_3();
        // for(int k=0; k<1024; k++){
        //     cout<<-trace[1][k]+temp1[k]<<endl;
        // }
        // getchar();
      }

//***** DLED
        if(DLED_bool){
            DLED(trace_length,dleddt);
            //Now: trace_DLED
        }else{
            for(ii=0; ii<trace_DLED_length; ii++){
                trace_DLED[0][ii] = trace[0][ii+dleddt];
                trace_DLED[1][ii] = trace[1][ii+dleddt];
            }
        }

        //--------------------------------------
        //---[ FIND PEAKS DISCRIMINATOR WAY ]---
        //--------------------------------------
        if(find_peaks_discriminator_bool and !find_peaks_bool){
            discriminator(thr_to_find_peaks, trace_DLED, trace_DLED_length);
            DCR_from_discriminator += (double)discriminator_cnt / (trace[0][trace_length-1]*TMath::Power(10,-9)); // time in trace is in ns
        }


        // FIND PEAKS LED

        if(find_peak_in_a_selected_window){

          int led_mintp;
          int led_maxtp;
          int num_windows = 1;

          if(led_and_dcr_0pe) {
            num_windows = 2;
          }
          else {
            num_windows = 1;
          }

            for(int jj = 0; jj < num_windows; jj++){ // LOOP ON NUM WINDOWS (ONE FOR LED, TWO FOR LED AND DCR)

                    // min and max value search
                    bool found_min = false;
                    bool found_max = false;
                    // mintp = (int)(1024*min_peak_window/trace[0][trace_length-1] - dleddt);
                    // maxtp = (int)(1024*max_peak_window/trace[0][trace_length-1] - dleddt);
                    for(int i=dleddt; i<max_peak_window[jj]*2; i++){
                      if((trace[0][i]>min_peak_window[jj]) && !found_min){
                        mintp = i-dleddt;
                        found_min = true;
                      }
                      if((trace[0][i]>max_peak_window[jj]) && !found_max){
                        maxtp = i-dleddt;
                        found_max = true;
                      }
                    } // end for
                    if(!found_min) {
                      mintp = min_peak_window[jj] - dleddt;
                      cerr << "################################################################################"<<endl;
                      cerr << "ERROR: min_peak_window out of value" << endl;
                      cerr << "################################################################################"<<endl;
                    }
                    if(!found_max){
                      mintp = max_peak_window[jj] - dleddt;
                      cerr << "################################################################################"<<endl;
                      cerr << "ERROR: max_peak_window out of value" << endl;
                      cerr << "################################################################################"<<endl;
                    }

                    if(jj==0){
                      led_mintp = mintp;
                      led_maxtp = maxtp;
                    }


                    // cout << "trace_DLED[0][mintp]\t" << trace_DLED[0][mintp] << endl;
                    // cout << "trace_DLED[0][maxtp]\t" << trace_DLED[0][maxtp] << endl;

                    if(jj==0){
                      min_line = trace_DLED[0][mintp];
                      max_line = trace_DLED[0][maxtp];
                    }
                    // cout<<"mintp "<<mintp<<"\t"<<"maxtp "<<maxtp<<"\t";
                    // cout<<"min_line "<<min_line<<"\t"<<"max_line "<<max_line<<"\t";
                    // cout<<"trace_DLED[0][0] "<<trace_DLED[0][0]<<endl;
                    index_for_peak = find_peak_fix_time(mintp, maxtp);
                    if(jj==0){
                      index_peak_LED = index_for_peak;
                    }
                    peak_rejected = true;
                    if(index_for_peak>-1){ // peak rejected or not
                      if(trace_DLED[1][index_for_peak]<thr_to_find_peaks){ // 0pe
                          peak_rejected = false;
                      } // end 0pe
                      else{ // 1pe or more
                          if(  trace_DLED[1][index_for_peak]>trace_DLED[1][mintp] && trace_DLED[1][index_for_peak]>trace_DLED[1][mintp-1]){ // check mintp
                            if(trace_DLED[1][index_for_peak]>trace_DLED[1][maxtp] && trace_DLED[1][index_for_peak]>trace_DLED[1][maxtp+1]){ // check maxtp
                              peak_rejected = false;
                            } // end check maxtp
                          } // end check mintp
                      } // end 1pe or more
                    }// end peak rejected or not

                    // peak OK
                    if(!peak_rejected){
                      peak_LED[0] = trace_DLED[0][index_for_peak];
                      peak_LED[1] = trace_DLED[1][index_for_peak];
                      if(jj==0) {
                        // cout<<mintp<<"\t"<<maxtp<<endl;
                        ptrHistLED->Fill(peak_LED[1]);
                      }
                      if(jj==1){
                        // cout<<mintp<<"\t"<<maxtp<<endl;
                        ptrHistDCR_window->Fill(peak_LED[1]);
                      }
                    } // end peak OK


                } // END LOOP ON NUM WINDOWS

                // reuse the old values: (for LED)
                mintp = led_mintp;
                maxtp = led_maxtp;


        } // end FIND PEAKS LED


        /////////////////////////
        ///   FIND AREA LED   ///
        /////////////////////////
        if(find_area_trace_bool){
          ptrHistArea->Fill(find_area_trace(peak_LED[0], peak_LED[0]-time_area_low, peak_LED[0]+time_area_high));
          // cout<<find_area_trace(peak_LED[0], peak_LED[0]-time_area_low, peak_LED[0]+time_area_high)<<endl;
        }



        /////////////////////
        /// PEAKS FINDING ///
        /////////////////////

        // if(n_ev==0)  {
        //     trace_time = ( trace_DLED[0][trace_DLED_length-1] - trace_DLED[0][0] ); // time in trace is in ns
        // }
        fill_hist_peaks_when_found = false;
        if(find_peaks_bool){
            // find_peaks(thr_to_find_peaks,max_peak_width, min_peak_width,blind_gap,DCR_DELAYS_bool);

            FindPeaksRisingFalling(thr_to_find_peaks, trace_DLED, trace_DLED_length, max_peak_width,2,4);

            if(n_ev==0){
                trace_time = 0;
                trace_time_raw = 0;
            }

            trace_time_raw += trace_DLED[0][trace_DLED_length-1] - trace_DLED[0][0];

            trace_time += ( trace_DLED[0][trace_DLED_length-1] - trace_DLED[0][0] - DCR_cnt_temp * 2 * dleddt - rise_time ); // time in trace is in ns


            // FindPeakPositions(trace_DLED[1], DLED_bool, dleddt);
            // FindPeaksFromPositions();
        }



//***** REMOVE OFFSET DLED
        if(DLED_offset_remove_bool){
            DLED_offset_remove();
        }

        if(!fill_hist_peaks_when_found){
            for(int i=0; i<num_peaks; i++){
                ptrHistAllPeaks[nfile]->Fill(trace_DLED[1][index_vect[i]]);
            }
        }

	//***** FIND CHARGE
        if(find_charge_window_bool){
            mintp = (int)(1024*minLED_charge/trace[0][trace_length-1]);
            maxtp = (int)(1024*maxLED_charge/trace[0][trace_length-1]);
            find_charge_selected_window(mintp, maxtp);
        }

//if(index_vect[1] > maxtp && ((trace[0][index_vect[1]] - trace[0][index_vect[0]]) > 20) && ((trace[1][index_vect[0]]-trace[1][minLED_charge]) > -22) && ((trace[1][index_vect[0]]-trace[1][minLED_charge]) < -18) && count_peak_window < 20)
//&& (abs(trace[1][minLED_charge] - trace[1][index_vect[0]]) > 18) && (abs(trace[1][minLED_charge] - trace[1][index_vect[0]]) < 22)

//cout << n_ev << endl;

//     if(find_1phe_bool){ // if find_1phe_bool
// 	if(index_vect[1] > maxtp && index_vect[0] > mintp && index_vect[0] < maxtp  && count_peak_window < 100 && ((trace[0][index_vect[1]+dleddt] - trace[0][index_vect[0]+dleddt]) > 20) && (abs(trace[1][mintp] - trace[1][index_vect[0]+dleddt]) > 18) && (abs(trace[1][mintp] - trace[1][index_vect[0]+dleddt]) < 22) ){
//
// 	 if(count_peak_window==0){
//
// 	      AVG_trace_window = new float*[2];
//               for(int ii = 0; ii < 2; ii++) {
//                 AVG_trace_window[ii] = new float[trace_window_length];
//               } // end for loop for AVG trace creation
//
//
//           for(int jj=0; jj < trace_window_length; jj++){
//              AVG_trace_window[0][jj]=trace[0][jj+mintp];
//              AVG_trace_window[1][jj]=0;
//           } // end for loop for AVG trace initialization
//          } // end if n_ev == 0
//
// // 	cout << "n_ev\t"<< n_ev << endl;
// // 	cout << " " << endl;
// // 	cout << "index_vect[0]+dleddt\t" <<index_vect[0]+dleddt << endl;
// // 	cout << "mintp\t"<<mintp << endl;
// // 	cout << "trace[0][index_vect[0]+dleddt]\t"<<trace[0][index_vect[0]+dleddt] << endl;
// // 	cout << "trace[1][mintp]\t"<<trace[1][mintp] << endl;
// // 	cout << "trace[1][index_vect[0]+dleddt]\t"<<trace[1][index_vect[0]+dleddt] << endl;
// // 	cout << "(trace[1][mintp+dleddt] - trace[1][index_vect[0]+dleddt])\t"<<(trace[1][mintp+dleddt] - trace[1][index_vect[0]+dleddt]) << endl;
// // 	cout << "(trace[1][mintp] - trace[1][index_vect[0]])\t"<<(trace[1][mintp] - trace[1][index_vect[0]]) << endl;
// // 	//cout <<  << endl;
//
//
// 	count_peak_window++;
//
// 	 for (int kk=0; kk < trace_window_length; kk++){
//                 AVG_trace_window[1][kk] += (trace[1][kk+index_vect[0]-dleddt]-trace[1][index_vect[0]-dleddt]);
//             } // end for loop a
//
// 	if(count_peak_window == 100){
//
// 	  for (int ww=0; ww < trace_window_length; ww++){
//              AVG_trace_window[1][ww] /= count_peak_window;
//           } // end for loop average
//       TCanvas *AVG_trace_canvas = new TCanvas("AVG_trace_canvas", "AVG_trace_canvas", 5);
// 	  show_AVG_trace_window(AVG_trace_canvas, AVG_trace_window[0], AVG_trace_window[1], trace_window_length, true);
// 	} // end if count == 20
//
//
// 	} //end if loop check 1phe waveform
// } // end if find_1phe_bool


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
        delete []peak_LED;

//================================================================================
// Please do not modify below, until the end of the function

      }//end loop boards


   }//loop events
  n_ev_tot = n_ev + 1;
  cout<<"Last event "<<n_ev_tot<<endl;
  fclose(f);

  double err_DCR_cnt = 0.;
  double err_trace_time = 0.;

  DCR_from_discriminator = DCR_from_discriminator / (double)n_ev_tot;

  if(find_peaks_bool and DCR_from_cnt_bool){
      // DCR from cnt
      trace_time *= TMath::Power(10,-9); // trace time is in ns
      trace_time /= n_ev_tot;
      DCR_from_cnt = (double)DCR_cnt / (trace_time * n_ev_tot);

      // errors
      err_trace_time = TMath::Sqrt( (5*n_ev_tot) * TMath::Power(10,-9) * TMath::Power(10,-9) );
      err_DCR_cnt = TMath::Sqrt((double)DCR_cnt);
      errDCR_from_cnt = TMath::Sqrt( TMath::Power( err_DCR_cnt / (trace_time * n_ev_tot), 2) + TMath::Power( (double)DCR_cnt*err_trace_time / ( trace_time * n_ev_tot * trace_time * n_ev_tot ), 2) );
      // errDCR_from_cnt /= trace_time;
      // errDCR_from_cnt /= n_ev_tot;

  }

//   if(find_1phe_bool){ // if find_1phe_bool
//   cout << "Number of traces with one peak within window: " << count_peak_window << endl;
//   cout << index_vect[1] << " " << mintp << " " << maxtp << endl;
//
//   if(count_peak_window < 100 && count_peak_window != 0){
//
// 	cout << "Qui!" << endl;
//
// 	  for (int zz=0; zz < trace_window_length; zz++){
//              AVG_trace_window[1][zz] /= count_peak_window;
//           }
//
// 	  show_AVG_trace_window(AVG_trace_canvas, AVG_trace_window[0], AVG_trace_window[1], trace_window_length, true);
// 	}
// } // end if find_1phe_bool

}

//------------------------------------------------------------------------------
void ReadRootFile(string filename, int last_event_n, bool display, TCanvas *c){

  // TFile *file_root = TFile::Open(filename.c_str());
  // TTree *wave      = (TTree *)file_root->Get("Waveforms");
  //
  // trace_length = 1024;
  //
  // Int_t nentries = wave->GetEntries();
  //
  // if (last_event_n > nentries)
  // {
  //   cout << "Last event is beyond number of entries of the TTree. Setting last event to number of TTree entries." << endl;
  //   last_event_n = nentries;
  // }
  //
  // // create traces. they will be filled by the traces inside the root file
  //
  // trace = new float*[2];
  // for(i = 0; i < 2; i++) {
  //   trace[i] = new float[trace_length];
  // }
  //
  // wave->SetBranchAddress("time_ch_2",trace[0]);
  // wave->SetBranchAddress("volt_ch_2",trace[1]);
  //
  //
  //
  // // create average trace if needed
  //
  // if(average){
  //   trace_AVG = new float*[2];
  //   for(i = 0; i < 2; i++) {
  //     trace_AVG[i] = new float[trace_length];
  //   }
  // } // end if average
  //
  // // create DLED trace
  //
  // trace_DLED_length = trace_length - dleddt;
  //
  // trace_DLED = new float*[2];
  // for(ii = 0; ii < 2; ii++) {
  //   trace_DLED[ii] = new float[trace_DLED_length];
  // }
  //
  // // loop over ttree entries. Each entry is a trace.
  //
  // for (int entry = 0; entry < last_event_n; ++entry){
  //
  //   if(entry%1000==0) cout<<"Read event "<<entry<<endl;
  //
  //   wave->GetEntry(entry);
  //
  //   /*for(int i=0; i<trace_length; i++){
  //     if((trace[1][i]<-0.09) && (trace[1][i]>-0.11)) cout<<"Trace 0.01"<<endl;
  //     printf("%.5lf\t%.5lf\n",trace[0][i],trace[1][i]);
  //     cout<<trace[0][i]<<"\t"<<trace[1][i]<<endl;
  //   } */
  //
  //   //REMOVE PEAK AT 0
  //   if(remove_0_peak_bool){
  //       remove_peak_0_half();
  //   }
  //
  //
  //
  //   for (int k = 0; k < trace_length; ++k)
  //   {
  //       if (reverse_bool){
  //       trace[1][k] = -trace[1][k];
  //       }
  //
  //   } // end for loop over waveform
  //
  //   if(display){
  //     if(!display_one_ev) display_peaks_now = display_peaks;
  //     if(display_one_ev and (entry==ev_to_display)) display_peaks_now = display_peaks;
  //   } // end if display
  //
  //   // DLED
  //
  //   if(DLED_bool){
  //     DLED(trace_length,dleddt);
  //           //Now: trace_DLED
  //   }
  //   else
  //   {
  //     for(ii=0; ii<trace_DLED_length; ii++){
  //       trace_DLED[0][ii] = trace[0][ii];
  //       trace_DLED[1][ii] = trace[1][ii];
  //     }
  //   } // end if DLED
  //
  //   if(find_peak_in_a_selected_window){
  //     float *peak = new float[2];
  //     index_for_peak = find_peak_fix_time(mintp, maxtp);
  //     peak[0] = trace_DLED[0][index_for_peak];
  //     peak[1] = trace_DLED[1][index_for_peak];
  //     ptrHistLED->Fill(peak[1]);
  //     delete[] peak;
  //   } // end if find peak in selected window
  //
  //   // PEAKS FINDING
  //   if(find_peaks_bool){
  //     find_peaks(thr_to_find_peaks,max_peak_width, min_peak_width,blind_gap,DCR_DELAYS_bool);
  //   } // end if peak finding
  //
  //   // FIND OFFSET
  //   if(find_offset_bool){
  //     find_offset_mod_3();
  //     subtract_offset();
  //   } // end if offset
  //
  //
  //   for(int i=start_blind_gap; i<trace_length-end_bling_gap; i++){
  //       ptrAllTrace->Fill(trace[1][i]);
  //   }
  //
  //   //***** FIND CHARGE for LED
  //   if(find_charge_window_bool){
  //     mintp = (int)(1024*minLED_charge/trace[0][trace_length-1]);
  //     maxtp = (int)(1024*maxLED_charge/trace[0][trace_length-1]);
  //     find_charge_selected_window(mintp, maxtp);
  //   } // end if find charge
  //
  //   //***** AVERAGE
  //   if(average){
  //
  //     if(DLED_bool){
  //
  //       if(entry==0){
  //           for(i=0; i< trace_DLED_length; i++){
  //               trace_AVG[0][i]=trace_DLED[0][i];
  //               trace_AVG[1][i]=0;
  //           }
  //       }
  //       for (ii=0; ii<trace_DLED_length; ii++){
  //           trace_AVG[1][ii] = trace_AVG[1][ii] + trace_DLED[1][ii];
  //       }
  //     }
  //     else{
  //
  //       if(entry==0){
  //           for(i=0; i< trace_length; i++){
  //               trace_AVG[0][i]=trace[0][i];
  //               trace_AVG[1][i]=0;
  //           }
  //       }
  //       for (ii=0; ii<trace_length; ii++){
  //           trace_AVG[1][ii] = trace_AVG[1][ii] + trace[1][ii];
  //       }
  //     }
  //   } // end if average
  //
  //
  //   //***** DISPLAY
  //   if(display){
  //       miny1 = -100; maxy1 = 100; miny2 = -100; maxy2 = 100;
  //       if(!display_one_ev){
  //           if(entry==0){
  //               c->Divide(1,2);
  //               c->SetGrid();
  //           }
  //           show_trace2(c, trace[0], trace[1], trace_DLED[0], trace_DLED[1], trace_length, trace_DLED_length, miny1, maxy1, miny2, maxy2, line_bool, true);
  //       }else{
  //           if(entry==ev_to_display){
  //               c->Divide(1,2);
  //               c->SetGrid();
  //               show_trace2(c, trace[0], trace[1], trace_DLED[0], trace_DLED[1], trace_length, trace_DLED_length, miny1, maxy1, miny2, maxy2, line_bool, false);
  //           }
  //       }
  //   } // end if display
  //
  // } // end entry loop
  //
  // cout << endl << "Last event " << last_event_n <<endl;
  //
  // n_ev_tot = last_event_n;
  //
  // delete []trace[0];
  // delete []trace[1];
  // delete []trace_DLED[0];
  // delete []trace_DLED[1];
  // delete []peak;
  //
  // file_root->Close();
  // delete file_root;
}


//------------------------------------------------------------------------------
//---------------------------------[   HELP   ]---------------------------------
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
void help(){
    cout<<endl;
    cout<<"Ana_Traces_SiPM.cxx"<<endl;
    cout<<"PREDEFINED FUNCTIONS:"<<endl;

    cout<<"\tvoid DCR_CT_1SiPM_1HV_all_delays(string file1, int last_event_n);"<<endl;
    cout<<"\tvoid DCR_CT_1SiPM_3HVs_all_delays(string file1, string file2, string file3, int last_event_n);"<<endl;
    cout<<"\tvoid Ana1(string file1, int last_event_n, bool display_one_ev_param);"<<endl;
    cout<<"\tvoid Ana_LED(string file1, int last_event_n);"<<endl;
    cout<<"\tvoid Ana_Ped(string file1, int last_event_n);"<<endl;

    cout<<endl;
    cout<<"See Ana_Traces_SiPM_ReadMe.md for more information"<<endl<<endl;
    cout<<"Davide Depaoli 2018"<<endl;
    cout<<endl;
}




//------------------------------------------------------------------------------
//----------------------------[   OLD FUNCTIONS   ]-----------------------------
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
void DCR_CT_1SiPM_1HV_NO_Delays(string file1, int last_event_n){
    //TRUE:
    find_peaks_bool = true;
    DCR_from_cnt_bool = true;
    DO_NOT_DELETE_HIST_LED = true;

    smooth_trace_bool = true;

    nfile = 0; //I only consider 1 file

    TCanvas *c = new TCanvas("Trace","Trace",w,h);

    ptrHistAllPeaks[0]  = new TH1D("histAllPeaks","",bins_DCR,0,maxyhistAllPeaks);
    ptrHistDelays[0]    = new TH1D("histDelays","",bins_Delays,0,maxyhistDelays);
    ptrHistDCRthr[0]    = new TH1D("histDCRthr","",bins_DCR,0,maxyhistDCR);


    // DCR_func
    gDCR_1 = DCR_func_NO_Delays(file1,last_event_n, 1, c);


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

}


//------------------------------------------------------------------------------
void DCR_CT_1SiPM_3HVs_NO_Delays(string file1, string file2, string file3, int last_event_n){
    //TRUE:
    find_peaks_bool = true;
    DCR_from_cnt_bool = true;
    DO_NOT_DELETE_HIST_LED = true;

    smooth_trace_bool = true;

    //I have 3 files
    nfiletot = 3;

    color_file_1 = kGreen+1;
    color_file_2 = kBlue;
    color_file_3 = kRed;

    TCanvas *c = new TCanvas("Trace","Trace",w,h);

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
    ind_peaks_all = 0;
    for(int i=0; i<max_peaks; i++){peaks_all[i] = 0;}
    min_thr_to_find_peaks = min_thr_to_find_peaks;
    // DCR_func
    gDCR_1 = DCR_func_NO_Delays(file1,last_event_n, nfiletot, c);

    //file2:
    nfile = 1;
    first_time_main_called = true; //will be set to false after the Analysis function is called
    ind_peaks_all = 0;
    for(int i=0; i<max_peaks; i++){peaks_all[i] = 0;}
    min_thr_to_find_peaks = min_thr_to_find_peaks;
    // DCR_func
    gDCR_2 = DCR_func_NO_Delays(file2,last_event_n, nfiletot, c);

    //file3:
    nfile = 2;
    first_time_main_called = true; //will be set to false after the Analysis function is called
    ind_peaks_all = 0;
    for(int i=0; i<max_peaks; i++){peaks_all[i] = 0;}
    min_thr_to_find_peaks = min_thr_to_find_peaks;
    // DCR_func
    gDCR_3 = DCR_func_NO_Delays(file3,last_event_n, nfiletot, c);

    TMultiGraph *DCR_mg = new TMultiGraph("DCR_mg", ";THR (mV); DCR (Hz)");
    DCR_mg->Add(gDCR_1);
    DCR_mg->Add(gDCR_2);
    DCR_mg->Add(gDCR_3);

    TCanvas *cDCR_loop = new TCanvas("cDCR_loop", "cDCR_loop");

    cDCR_loop->SetGrid();
    cDCR_loop->SetLogy();
    DCR_mg->Draw("A3L");

    auto legendDCR_loop = new TLegend(0.75,0.75,0.9,0.9);
    legendDCR_loop->AddEntry(gDCR_1,"DCR file 1","l");
    legendDCR_loop->AddEntry(gDCR_2,"DCR file 2","l");
    legendDCR_loop->AddEntry(gDCR_3,"DCR file 3","l");
    // legendDCR_loop->AddEntry(gDCR_1,"HV = 34.00 V","l");
    // legendDCR_loop->AddEntry(gDCR_2,"HV = 35.00 V","l");
    // legendDCR_loop->AddEntry(gDCR_3,"HV = 36.00 V","l");
    legendDCR_loop->Draw();

    cout<<endl<<endl;
    cout<<"-------------------------"<<endl;
    cout<<"-------[ RESULTS ]-------"<<endl;
    cout<<"-------------------------"<<endl<<endl;


    // CROSS TALK
    double CT[3][2];   // CT[nfile][0] = CrossTalk, CT[nfile][1] = errCrossTalk
    for(int i=0; i<nfiletot; i++){
        EvaluateCrossTalk(DCR_pe_0_5_vect[i], errDCR_pe_0_5_vect[i], DCR_pe_1_5_vect[i], errDCR_pe_1_5_vect[i], CT[i]);
    }

    cout<<"// "<<"Files analyzed:"<<endl;
    cout<<"// "<<file1<<"   pe_0_5 = "<<pe_0_5_vect[0]<<" mV; pe_1_5 = "<<pe_1_5_vect[0]<<" mV"<<endl;
    cout<<"// "<<file2<<"   pe_0_5 = "<<pe_0_5_vect[1]<<" mV; pe_1_5 = "<<pe_1_5_vect[1]<<" mV"<<endl;
    cout<<"// "<<file3<<"   pe_0_5 = "<<pe_0_5_vect[2]<<" mV; pe_1_5 = "<<pe_1_5_vect[2]<<" mV"<<endl;
    cout<<"double DCR[] =          {"<<DCR_pe_0_5_vect[0]*TMath::Power(10,-6)<<", "<<DCR_pe_0_5_vect[1]*TMath::Power(10,-6)<<", "<<DCR_pe_0_5_vect[2]*TMath::Power(10,-6)<<"};"<<endl;

    cout<<"double errDCR[] =       {"<<errDCR_pe_0_5_vect[0]*TMath::Power(10,-6)<<", "<<errDCR_pe_0_5_vect[1]*TMath::Power(10,-6)<<", "<<errDCR_pe_0_5_vect[2]*TMath::Power(10,-6)<<"};"<<endl;

    cout<<"double CrossTalk[] =    {"<<CT[0][0]<<", "<<CT[1][0]<<", "<<CT[2][0]<<"};"<<endl;

    cout<<"double errCrossTalk[] = {"<<CT[0][1]<<", "<<CT[1][1]<<", "<<CT[2][1]<<"};"<<endl;

}


//------------------------------------------------------------------------------
void DCR_CT_1SiPM_5HVs_NO_Delays(string filelist, int last_event_n){
    //TRUE:
    find_peaks_bool = true;
    DCR_from_cnt_bool = true;
    DO_NOT_DELETE_HIST_LED = true;

    smooth_trace_bool = true;


    ifstream OpenFile (filelist.c_str());
    // char file1[300], file2[300], file3[300], file4[300], file5[300];
    string file1, file2, file3, file4, file5;
    OpenFile>>file1;
    OpenFile>>file2;
    OpenFile>>file3;
    OpenFile>>file4;
    OpenFile>>file5;

    cout<<"Analysing files:"<<endl;
    cout<<"file1 = "<<file1<<endl;
    cout<<"file2 = "<<file2<<endl;
    cout<<"file3 = "<<file3<<endl;
    cout<<"file4 = "<<file4<<endl;
    cout<<"file5 = "<<file5<<endl;

    OpenFile.close();


    //I have 3 files
    nfiletot = 5;

    TCanvas *c = new TCanvas("Trace","Trace",w,h);

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

    // colors:
    color_file_1 = kGreen;
    color_file_2 = kGreen+1;
    color_file_3 = kBlue;
    color_file_4 = kRed;
    color_file_5 = kRed+1;

    opacity = 0.3;

    //file1:
    nfile = 0;
    first_time_main_called = true; //will be set to false after the Analysis function is called
    ind_peaks_all = 0;
    for(int i=0; i<max_peaks; i++){peaks_all[i] = 0;}
    // DCR_func
    gDCR_1 = DCR_func_NO_Delays(file1,last_event_n, nfiletot, c);

    //file2:
    nfile = 1;
    first_time_main_called = true; //will be set to false after the Analysis function is called
    ind_peaks_all = 0;
    for(int i=0; i<max_peaks; i++){peaks_all[i] = 0;}
    // DCR_func
    gDCR_2 = DCR_func_NO_Delays(file2,last_event_n, nfiletot, c);

    //file3:
    nfile = 2;
    first_time_main_called = true; //will be set to false after the Analysis function is called
    ind_peaks_all = 0;
    for(int i=0; i<max_peaks; i++){peaks_all[i] = 0;}
    // DCR_func
    gDCR_3 = DCR_func_NO_Delays(file3,last_event_n, nfiletot, c);

    //file4:
    nfile = 3;
    first_time_main_called = true; //will be set to false after the Analysis function is called
    ind_peaks_all = 0;
    for(int i=0; i<max_peaks; i++){peaks_all[i] = 0;}
    // DCR_func
    gDCR_4 = DCR_func_NO_Delays(file4,last_event_n, nfiletot, c);

    //file5:
    nfile = 4;
    first_time_main_called = true; //will be set to false after the Analysis function is called
    ind_peaks_all = 0;
    for(int i=0; i<max_peaks; i++){peaks_all[i] = 0;}
    // DCR_func
    gDCR_5 = DCR_func_NO_Delays(file5,last_event_n, nfiletot, c);

    TMultiGraph *DCR_mg = new TMultiGraph("DCR_mg", ";THR (mV); DCR (Hz)");
    DCR_mg->Add(gDCR_1);
    DCR_mg->Add(gDCR_2);
    DCR_mg->Add(gDCR_3);
    DCR_mg->Add(gDCR_4);
    DCR_mg->Add(gDCR_5);

    TCanvas *cDCR_loop = new TCanvas("cDCR_loop", "cDCR_loop");

    cDCR_loop->SetGrid();
    cDCR_loop->SetLogy();
    DCR_mg->Draw("A3L");

    auto legendDCR_loop = new TLegend(0.75,0.75,0.9,0.9);
    legendDCR_loop->AddEntry(gDCR_1,"DCR file 1","l");
    legendDCR_loop->AddEntry(gDCR_2,"DCR file 2","l");
    legendDCR_loop->AddEntry(gDCR_3,"DCR file 3","l");
    legendDCR_loop->AddEntry(gDCR_4,"DCR file 4","l");
    legendDCR_loop->AddEntry(gDCR_5,"DCR file 5","l");

    legendDCR_loop->Draw();

    // cout<<endl<<endl;
    // cout<<"-------------------------"<<endl;
    // cout<<"-------[ RESULTS ]-------"<<endl;
    // cout<<"-------------------------"<<endl<<endl;
    //
    //
    // // CROSS TALK
    // double CT[3][2];   // CT[nfile][0] = CrossTalk, CT[nfile][1] = errCrossTalk
    // for(int i=0; i<nfiletot; i++){
    //     EvaluateCrossTalk(DCR_pe_0_5_vect[i], errDCR_pe_0_5_vect[i], DCR_pe_1_5_vect[i], errDCR_pe_1_5_vect[i], CT[i]);
    // }
    //
    // cout<<"// "<<"Files analyzed:"<<endl;
    // cout<<"// "<<file1<<"   pe_0_5 = "<<pe_0_5_vect[0]<<" mV; pe_1_5 = "<<pe_1_5_vect[0]<<" mV"<<endl;
    // cout<<"// "<<file2<<"   pe_0_5 = "<<pe_0_5_vect[1]<<" mV; pe_1_5 = "<<pe_1_5_vect[1]<<" mV"<<endl;
    // cout<<"// "<<file3<<"   pe_0_5 = "<<pe_0_5_vect[2]<<" mV; pe_1_5 = "<<pe_1_5_vect[2]<<" mV"<<endl;
    // cout<<"double DCR[] =          {"<<DCR_pe_0_5_vect[0]*TMath::Power(10,-6)<<", "<<DCR_pe_0_5_vect[1]*TMath::Power(10,-6)<<", "<<DCR_pe_0_5_vect[2]*TMath::Power(10,-6)<<"};"<<endl;
    //
    // cout<<"double errDCR[] =       {"<<errDCR_pe_0_5_vect[0]*TMath::Power(10,-6)<<", "<<errDCR_pe_0_5_vect[1]*TMath::Power(10,-6)<<", "<<errDCR_pe_0_5_vect[2]*TMath::Power(10,-6)<<"};"<<endl;
    //
    // cout<<"double CrossTalk[] =    {"<<CT[0][0]<<", "<<CT[1][0]<<", "<<CT[2][0]<<"};"<<endl;
    //
    // cout<<"double errCrossTalk[] = {"<<CT[0][1]<<", "<<CT[1][1]<<", "<<CT[2][1]<<"};"<<endl;

}


//------------------------------------------------------------------------------
void DCR_discriminator(string file1, int last_event_n, int thr_to_find_peaks){
    // TRUE:
    find_peaks_discriminator_bool = true;

    smooth_trace_bool = true;

    TCanvas *c = new TCanvas("Trace","Trace",w,h);

    //Analysis
    Analysis(file1, last_event_n, false, c);

    delete c;

    cout<<endl;
    cout<<"// File analyzed: "<<file1<<endl;
    cout<<"// threshold = "<<thr_to_find_peaks<<endl;
    cout<<"DCR = "<< DCR_from_discriminator*TMath::Power(10,-6) <<" MHz"<<endl;
}

//------------------------------------------------------------------------------
void DCR_CT_No_Stair(string file1, int last_event_n, float thr_0_5_pe, float thr_1_5_pe){
    //VARIABLES:
    //TRUE:
    find_peaks_bool = true;
    DCR_DELAYS_bool = true; //DCR from delays
    fit_hist_del_bool = false;
    DO_NOT_DELETE_HIST_LED = true;
    DCR_from_cnt_bool = true;

    smooth_trace_bool = true;


    nfile = 0; //I only consider 1 file

    ptrHistAllPeaks[0]  = new TH1D("histAllPeaks","",bins_DCR,0,maxyhistAllPeaks);
    ptrHistDelays[0]    = new TH1D("histDelays","",bins_Delays,0,maxyhistDelays);
    ptrHistDCRthr[0]    = new TH1D("histDCRthr","",bins_DCR,0,maxyhistDCR);

    TCanvas *c = new TCanvas("Trace","Trace",w,h);

    double DCR_from_cnt_0_5_pe = 0.;
    double errDCR_from_cnt_0_5_pe = 0.;
    double DCR_from_cnt_1_5_pe = 0.;
    double errDCR_from_cnt_1_5_pe = 0.;

    //////////////
    /// 0.5 pe ///
    //////////////

    // Analysis
    first_time_main_called = true;
    thr_to_find_peaks = thr_0_5_pe;
    Analysis(file1, last_event_n, false, c);

    DCR_from_cnt_0_5_pe = DCR_from_cnt*TMath::Power(10,-6);
    errDCR_from_cnt_0_5_pe = errDCR_from_cnt*TMath::Power(10,-6);

    Get_DCR_temp_and_errDCR_temp();
    DCR_pe_0_5_vect[nfile] = DCR_temp[nfile];
    errDCR_pe_0_5_vect[nfile] = errDCR_temp[nfile];


    // reset
    ptrHistAllPeaks[0]->Reset();
    ptrHistDelays[0]->Reset();
    ptrHistDCRthr[0]->Reset();
    DCR_cnt = 0;
    DCR_from_cnt = 0;
    trace_time = 0;
    n_ev_tot = 0;

    //////////////
    /// 1.5 pe ///
    //////////////

    // Analysis for 1.5 pe
    first_time_main_called = true;
    thr_to_find_peaks = thr_1_5_pe;
    Analysis(file1, last_event_n, false, c);

    delete c;

    DCR_from_cnt_1_5_pe = DCR_from_cnt*TMath::Power(10,-6);
    errDCR_from_cnt_1_5_pe = errDCR_from_cnt*TMath::Power(10,-6);

    Get_DCR_temp_and_errDCR_temp();
    DCR_pe_1_5_vect[nfile] = DCR_temp[nfile];
    errDCR_pe_1_5_vect[nfile] = errDCR_temp[nfile];


    // CROSS TALK
    double CT[2];   // CT[0] = CrossTalk, CT[1] = errCrossTalk
    EvaluateCrossTalk(DCR_pe_0_5_vect[0], errDCR_pe_0_5_vect[0], DCR_pe_1_5_vect[0], errDCR_pe_1_5_vect[0], CT);
    cout<<"   Cross Talk    = ("<<CT[0]<<" +- "<<CT[1]<<")"<<endl;



    float CrossTalk = DCR_pe_1_5_vect[0]/DCR_pe_0_5_vect[0];
    float errCrossTalk= CrossTalk * TMath::Sqrt( (errDCR_pe_0_5_vect[0]/DCR_pe_0_5_vect[0])*(errDCR_pe_0_5_vect[0]/DCR_pe_0_5_vect[0]) + (errDCR_pe_1_5_vect[0]/DCR_pe_1_5_vect[0])*(errDCR_pe_1_5_vect[0]/DCR_pe_1_5_vect[0]) );
    cout<<"   Cross Talk    = ("<<CrossTalk<<" +- "<<errCrossTalk<<")"<<endl;



    cout<<endl<<endl;
    cout<<"-------------------------"<<endl;
    cout<<"-------[ RESULTS ]-------"<<endl;
    cout<<"-------------------------"<<endl<<endl;

    cout<<"File analyzed: "<<file1<<endl;

    cout<<"   DCR at "<<thr_0_5_pe<<" mV = "<<endl;
    cout<<"         ("<<DCR_pe_0_5_vect[0]*TMath::Power(10,-6)<<" +- "<<errDCR_pe_0_5_vect[0]*TMath::Power(10,-6)<<") MHz, from exp fit"<<endl;
    cout<<"         ("<<DCR_from_cnt_0_5_pe<<" +- "<<errDCR_from_cnt_0_5_pe<<") MHz, from cnt"<<endl;

    cout<<"   DCR at "<<thr_1_5_pe<<" mV = "<<endl;
    cout<<"         ("<<DCR_pe_1_5_vect[0]*TMath::Power(10,-6)<<" +- "<<errDCR_pe_1_5_vect[0]*TMath::Power(10,-6)<<") MHz, from exp fit"<<endl;
    cout<<"         ("<<DCR_from_cnt_1_5_pe<<" +- "<<errDCR_from_cnt_1_5_pe<<") MHz, from cnt"<<endl;



}


//------------------------------------------------------------------------------
void DCR_CT_1SiPM_1HV_all_delays(string file1, int last_event_n){
    //TRUE:
    find_peaks_bool = true;
    fit_hist_del_bool = true;
    DCR_DELAYS_bool = true; //DCR from delays
    CROSS_TALK_bool = true; //DCR must be true
    DO_NOT_DELETE_HIST_LED = true;

    smooth_trace_bool = true;

    nfile = 0; //I only consider 1 file

    ptrHistDelays[0]    = new TH1D("histDelays","",bins_Delays,0,maxyhistDelays);
    ptrHistAllPeaks[0]  = new TH1D("histAllPeaks","",bins_DCR,0,maxyhistAllPeaks);
    ptrHistDCRthr[0]    = new TH1D("histDCRthr","",bins_DCR,0,maxyhistDCR);

    TCanvas *c = new TCanvas("Trace","Trace",w,h);


    // DCR_func
    gDCR_1 = DCR_func(file1,last_event_n, 1, c);

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
void DCR_CT_1SiPM_3HVs_all_delays(string file1, string file2, string file3, int last_event_n){

    //TRUE:
    find_peaks_bool = true;
    DCR_DELAYS_bool = true; //DCR from delays
    fit_hist_del_bool = true;
    CROSS_TALK_bool = true; //DCR must be true
    DO_NOT_DELETE_HIST_LED = true;

    smooth_trace_bool = true;

    //I have 3 files
    nfiletot = 3;

    color_file_1 = kGreen+1;
    color_file_2 = kBlue;
    color_file_3 = kRed;

    TCanvas *c = new TCanvas("Trace","Trace",w,h);

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
    gDCR_1 = DCR_func(file1,last_event_n, 3, c);


    //file2:
    nfile = 1;
    first_time_main_called = true; //will be set to false after the Analysis function is called
    ind_peaks_all_delay = 0;
    for(int i=0; i<max_peaks; i++){peaks_all_delay[0][i] = 0; peaks_all_delay[1][i] = 0;}
    // DCR_func
    gDCR_2 = DCR_func(file2,last_event_n, 3, c);


    //file3:
    nfile = 2;
    first_time_main_called = true; //will be set to false after the Analysis function is called
    ind_peaks_all_delay = 0;
    for(int i=0; i<max_peaks; i++){peaks_all_delay[0][i] = 0; peaks_all_delay[1][i] = 0;}
    // DCR_func
    gDCR_3 = DCR_func(file3,last_event_n, 3, c);

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

    cout<<"// "<<"Files analyzed:"<<endl;
    cout<<"// "<<file1<<"   pe_0_5 = "<<pe_0_5_vect[0]<<" mV; pe_1_5 = "<<pe_1_5_vect[0]<<" mV"<<endl;
    cout<<"// "<<file2<<"   pe_0_5 = "<<pe_0_5_vect[1]<<" mV; pe_1_5 = "<<pe_1_5_vect[1]<<" mV"<<endl;
    cout<<"// "<<file3<<"   pe_0_5 = "<<pe_0_5_vect[2]<<" mV; pe_1_5 = "<<pe_1_5_vect[2]<<" mV"<<endl;
    cout<<"double DCR[] =         {"<<DCR_pe_0_5_vect[0]*TMath::Power(10,-6)<<", "<<DCR_pe_0_5_vect[1]*TMath::Power(10,-6)<<", "<<DCR_pe_0_5_vect[2]*TMath::Power(10,-6)<<"};"<<endl;

    cout<<"double errDCR[] =      {"<<errDCR_pe_0_5_vect[0]*TMath::Power(10,-6)<<", "<<errDCR_pe_0_5_vect[1]*TMath::Power(10,-6)<<", "<<errDCR_pe_0_5_vect[2]*TMath::Power(10,-6)<<"};"<<endl;

    cout<<"double CrossTalk[] =   {"<<CrossTalk[0]<<", "<<CrossTalk[1]<<", "<<CrossTalk[2]<<"};"<<endl;

    cout<<"double errCrossTalk[] ={"<<errCrossTalk[0]<<", "<<errCrossTalk[1]<<", "<<errCrossTalk[2]<<"};"<<endl;

}

//
//
// //------------------------------------------------------------------------------
// void smooth_trace_step(){ // smooth the trace before DLED
//   for(int i=0; i<trace_length; i+=n_smooth_trace){
//
//     // sum on n_smooth_trace
//     for(int j=i+1; j<i+n_smooth_trace; j++){
//       trace[0][i] += trace[0][j];
//       trace[1][i] += trace[1][j];
//     }
//
//     // divide for n_smooth_trace
//     trace[0][i] /= n_smooth_trace;
//     trace[1][i] /= n_smooth_trace;
//
//     // all the values of the trace from i to i+n_smooth_trace are set equal
//     for(int j=i+1; j<i+n_smooth_trace; j++){
//       trace[0][j] = trace[0][i];
//       trace[1][j] = trace[1][i];
//     }
//   }
// }
//
//
// //------------------------------------------------------------------------------
// void smooth_trace_3(){ // smooth the trace before DLED
//   for(int i=0; i<trace_length; i+=3){
//
//     // sum on 3
//     for(int j=1; j<3; j++){
//       trace[1][i] += trace[1][i+j];
//     }
//
//     // divide for 3
//     trace[1][i]   /= 3;
//     for(int j=1; j<3; j++){
//       trace[1][i+j] = trace[1][i];
//     }
//
//     // now I have something like:
//     //
//     //                    [ i ] [i+1] [i+2]
//     //  [i-3] [i-2] [i-1]
//     //
//
//     // I change the values of [ i ] and [i-1] in order to smooth the trace
//     if(i!=0){
//       double gap = trace[1][i] - trace[1][i-1];
//       trace[1][i-1] += gap/3;
//       trace[1][i]   -= gap/3;
//     }
//     // cout<<i<<endl;
//   }
// }
//
//
// //------------------------------------------------------------------------------
// void smooth_trace_4(){ // smooth the trace before DLED
//   for(int i=0; i<trace_length; i+=4){
//
//     // sum on 4
//     for(int j=1; j<4; j++){
//       trace[1][i] += trace[1][i+j];
//     }
//
//     // divide for 4
//     trace[1][i]   /= 4;
//     for(int j=1; j<4; j++){
//       trace[1][i+j] = trace[1][i];
//     }
//
//     // now I have something like:
//     //
//     //                         [ i ] [i+1] [i+2] [i+3]
//     //  [i-4] [i-3] [i-2] [i-1]
//     //
//
//     // I change the values of [i-2], [i-1], [i] and [i+1] in order to smooth the trace
//     if(i!=0){
//       double gap = (double)trace[1][i+2] - (double)trace[1][i-2];
//       trace[1][i-2] += gap/5;
//       trace[1][i-1] += gap/5*2;
//       trace[1][i]   -= gap/5*2;
//       trace[1][i+1] -= gap/5;
//     }
//   }
// }
//
//
//
// //------------------------------------------------------------------------------
// void smooth_trace_5(){ // smooth the trace before DLED
//   for(int i=0; i<trace_length; i+=5){
//
//     // sum on 5
//     for(int j=1; j<5; j++){
//       trace[1][i] += trace[1][i+j];
//     }
//
//     // divide for 5
//     trace[1][i]   /= 5;
//     for(int j=1; j<5; j++){
//       trace[1][i+j] = trace[1][i];
//     }
//
//     // now I have something like:
//     //
//     //                               [ i ] [i+1] [i+2] [i+3] [i+4]
//     //  [i-5] [i-4] [i-3] [i-2] [i-1]
//     //
//
//     // I change the values of [i-2], [i-1], [i] and [i+1] in order to smooth the trace
//     if(i!=0){
//       double gap = (double)trace[1][i+2] - (double)trace[1][i-2];
//       trace[1][i-2] += gap/5;
//       trace[1][i-1] += gap/5*2;
//       trace[1][i]   -= gap/5*2;
//       trace[1][i+1] -= gap/5;
//     }
//   }
// }
//
//
// //------------------------------------------------------------------------------
// void SmoothTraceMarkov(int window){
//
//     double temp[trace_length];
//     double diff = 1000.0;
//     TSpectrum *s = new TSpectrum();
//     for(int i=0; i<trace_length; i++){
//        temp[i] = (double)trace[1][i] + diff;
//     }
//     s->SmoothMarkov(temp, trace_length, window);
//     for (int i=0; i<trace_length; i++) {
//        trace[1][i] = (float)temp[i] - diff;
//     }
//     delete s;
//
// }






//-----------------------------------------------------------------------------



# ifndef __CINT__
int main()
{
  cout<<"main"<<endl;
  return 0;
}
# endif
