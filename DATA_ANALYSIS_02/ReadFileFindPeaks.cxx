// ReadFileFindPeaks.cxx




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
#define max_peaks 50000000


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


// MAIN
void Analysis(string file, int last_event_n, double thr, bool display);

// SECONDARY
void DLED(int trace_length, int dleddt);
int FindPeaksFixTime(int mintp, int maxtp);
void FindPeaksRisingFalling(double thr, float **t, double length, int max_peak_width, int rising_points, int falling_points);
void show_trace(TCanvas* canv, float *x, float *y, int trace_length, float miny, float maxy, bool line_bool, bool delete_bool);
void show_trace2(TCanvas* canv, float *x1, float *y1, float *x2, float *y2, int trace_length1, int trace_length2, float miny1, float maxy1, float miny2, float maxy2, bool line_bool, bool delete_bool);
void FindDelaysFromVector();
void FindDCRfromVector();
double find_area_trace(double center, double low, double high);
void show_AVG_trace_window(TCanvas *c, float *tracet, float *tracev, int trace_length, bool delete_bool);
void DLED_offset_remove();
void SmoothTraceN(int n);
void WriteIntro();

// READ FILE
void ReadDRS4Bin(string filename, int last_event_n, bool display, TCanvas *c);

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

bool change_title_size_bool = true;

bool fit_trace_bool = false;

bool find_clean_peaks = false;

//-----------------
//-----------------


//---------------
//---[ TRACE ]---
//---------------
int start_blind_gap = 20;
int end_bling_gap = 100;
int n_smooth_trace = 2;

//---------------
//---[ PEAKS ]---
//---------------

// DLED and PEAKS FINDING
int dleddt = 6;//9;//8;//5;//9*GSPS;
    // dleddt = 6; for DCR_CT_1SiPM_nHVs(), 20180725_HD3-2_01_DARK_AgilentE3641A_35.00_AS_2_100000ev_01.dat and similar
    // dleddt = 9; for Ana_LED(), 20180614_HD3-2_1_LASER_PLS_81_PAPER_AGILENT_35_AS_2_50000_01.dat and similar
int blind_gap = 2*dleddt; //ns
int max_peak_width = 50; //used for find_peaks
int min_peak_width =  0; //used for find_peaks
int gap_between_peaks = 10;
int rise_time = dleddt;

// ONLY for LED measures
int minLED_amp = 168;//168;//168;//290;//115;  // window: min time for peak (ns) for LED
int maxLED_amp = 180;//177;//176;//305;//125;  // window: max time for peak (ns) for LED
double time_area_low = 30;  // time for the area before the LED peak (ns)
double time_area_high = 200; // time for the area after the LED peak (ns)

// threshold
float thr_to_find_peaks = 8; //thr_to_find_peaks, as seen in DLED trace (in V); it should be similar to pe_0_5.

//-----------------
//---[ DISPLAY ]---
//-----------------
int ev_to_display = 5;


//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------



//------------------------------------------------------------------------------
//---------------------------[   GLOBAL VARIABLES   ]---------------------------
//------------------------------------------------------------------------------
float **trace;
float **trace_DLED;
float **trace_AVG;
float **AVG_trace_window;
float *peak;
float *peaks;
float *peak_LED;

ofstream FilePeaks;
string FilePeaksName;

int trace_DLED_length; int index_func = 0; int nfile = 0; int n_DCR = 0; int DCR_cnt = 0; int DCR_cnt_temp = 0; int discriminator_cnt = 0; int trace_length = 0; int n_ev, index_for_peak; int one_window=0;
int nfiletot = 1; int n_smooth = 0;

float miny=0; float maxy=0; float miny1=0; float maxy1=0; float miny2=0; float maxy2=0;

double DCR_from_cnt = 0.;
double errDCR_from_cnt = 0.;

int num_peaks=0;
int index_vect[max_peak_num] = {0};
int mintp = 0; int maxtp = 0;
int index_peak_LED = 0;

float peaks_all_delay[2][max_peaks];

int ind_peaks_all_delay = 0;
float peaks_all[max_peaks];
float peaks_all_time[max_peaks];
int peaks_all_index[max_peaks];
int ind_peaks_all = 0;
int ind_peaks_all_old = 0;
int n_ev_tot = 0;

int clean_peaks_num = 0;

float offset = 0.; float charge = 0.;


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
//----------------------------[   MAIN FUNCTION   ]-----------------------------
//------------------------------------------------------------------------------

//______________________________________________________________________________
void Analysis(string file, int last_event_n, double thr, bool display){
    gROOT->Reset();

    thr_to_find_peaks = thr;

    find_peaks_bool = true;
    drawHistAllPeaks = true; // to draw hist of all peaks in traces
    fitHistAllPeaks = false; // fit hist of all peaks -> for GAIN
    DCR_DELAYS_bool = true; //DCR from delays
    fit_hist_del_bool = true;
    show_hists_DCR_DELAYS  = true;
    DO_NOT_DELETE_HIST_LED = true;
    display_one_ev = !display;
    display_peaks = true;
    DCR_from_cnt_bool = true;
    smooth_trace_bool = false;

    char thr_char[50];
    int temp_results;

    temp_results = sprintf(thr_char, "%.2f", thr_to_find_peaks);
    string thr_string(thr_char);

    std::size_t found = file.find_last_of("/\\");

    FilePeaksName = file.substr(0,found) + "/Peaks/" + file.substr(found+1) + "_Peaks_" + thr_string + ".txt";
    cout<<FilePeaksName<<endl;

    FilePeaks.open(FilePeaksName);

    WriteIntro();

    TCanvas *c = new TCanvas("Trace","Trace",800,600);

    ReadDRS4Bin(file, last_event_n, display, c);

    FilePeaks.close();

//***** AVERAGE
    if(average){
        for(int i=0; i<trace_length; i++){
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

}


//------------------------------------------------------------------------------
//-------------------------[   SECONDARY FUNCTIONS   ]--------------------------
//------------------------------------------------------------------------------

//______________________________________________________________________________
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


//______________________________________________________________________________
void DLED(int trace_length, int dleddt){
    for(int i=0; i<trace_DLED_length; i++){
        trace_DLED[0][i] = trace[0][i + dleddt];
        trace_DLED[1][i] = trace[1][i + dleddt]-trace[1][i];
    }
}

//______________________________________________________________________________
int FindPeakLED(int mintp, int maxtp){
    double max_func = -10000;
    index_func=-1;
    for(int i=mintp; i<maxtp; i++){
        if(trace_DLED[1][i]>max_func && trace_DLED[1][i]>trace_DLED[1][i-1] && trace_DLED[1][i]>trace_DLED[1][i+1]){
            max_func=trace_DLED[1][i];
            index_func=i;
        }
    }

    if(index_func==-1) {
        for(int i=mintp; i<maxtp; i++){
            if(trace_DLED[1][i]>max_func && trace_DLED[1][i]>trace_DLED[1][mintp-1] && trace_DLED[1][i]>trace_DLED[1][maxtp+1]){
                max_func=trace_DLED[1][i];
                index_func=i;
            }
        }
    }

    if(index_func==-1){
        for(int i=mintp; i<maxtp; i++){
            if(trace_DLED[1][i]>max_func){
                max_func=trace_DLED[1][i];
                index_func=i;
            }
        }
    }

    return index_func;
}



//______________________________________________________________________________
int FindPeaksFixTime(int mintp, int maxtp){
    double max_func = -10000;
    index_func=-1;
    for(int i=mintp; i<maxtp; i++){
        //cout<<mintp<<"\t"<<maxtp<<endl;
        if(trace_DLED[1][i]>max_func){
            max_func=trace_DLED[1][i];
            index_func=i;
        }
    }

    return index_func;

}


//______________________________________________________________________________
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

    ind_peaks_all_old = ind_peaks_all;

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

                    //Now I look for the maximum value in that window
                    index_new = FindPeaksFixTime(peak_start, peak_end);

                    // cout<<t[0][index_new]<<"\t"<<t[1][index_new]<<endl;
                    FilePeaks<<t[0][index_new]<<"\t"<<t[1][index_new]<<endl;

                    if(num_peaks<max_peak_num){
                        index_vect[num_peaks] = index_new;
                        num_peaks++;
                    }

                    // now index_new will be index_old, for the following loop cycle
                    index_old = index_new;




                    i+=gap_between_peaks;
                    //i++;
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

    FilePeaks<<"N"<<endl;


}


//______________________________________________________________________________
void show_trace(TCanvas* canv, float *x, float *y, int trace_length, float miny, float maxy, bool line_bool, bool delete_bool){

    if(reverse_bool){
      for(int i=0; i<trace_length; i++){
        y[i] = -y[i];
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


//______________________________________________________________________________
void show_trace2(TCanvas* canv, float *x1, float *y1, float *x2, float *y2, int trace_length1, int trace_length2, float miny1, float maxy1, float miny2, float maxy2, bool line_bool, bool delete_bool){
    if(reverse_bool){
      for(int i=0; i<trace_length1; i++){
        y1[i] = -y1[i];
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
    graph1->Draw("al");
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
    graph2->Draw("al");
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

    if(change_title_size_bool){
        gStyle->SetTitleFontSize(0.06);
        graph1->GetXaxis()->SetTitleSize(0.06);
        graph1->GetYaxis()->SetTitleSize(0.06);
        graph1->GetXaxis()->SetLabelSize(0.06);
        graph1->GetYaxis()->SetLabelSize(0.06);
        graph2->GetXaxis()->SetTitleSize(0.06);
        graph2->GetYaxis()->SetTitleSize(0.06);
        graph2->GetXaxis()->SetLabelSize(0.06);
        graph2->GetYaxis()->SetLabelSize(0.06);
        canv->Update();
    }

    canv->Update();

    //--------------------
    // new TCanvas();
    // TGraphErrors *graphtemp = new TGraphErrors(trace_length1,x1,y1,0,0);
    // graphtemp->Draw();
    //--------------------

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


//______________________________________________________________________________
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

//______________________________________________________________________________
double find_area_trace(double center, double low, double high){
    double area = 0;
    int low_index = 0;
    int high_index = 0;

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

//______________________________________________________________________________
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
//------------------------------[   READ FILES   ]------------------------------
//------------------------------------------------------------------------------

//______________________________________________________________________________
void ReadDRS4Bin(string filename, int last_event_n, bool display, TCanvas *c){
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
        for(int k = 0; k < 2; k++) {
            trace_DLED[k] = new float[trace_DLED_length];
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


        if(smooth_trace_bool){
          int n_SmootTraceN = 2;
          SmoothTraceN(n_SmootTraceN);
        }

        // DLED
        if(DLED_bool){
            DLED(trace_length,dleddt);
            //Now: trace_DLED
        }else{
            for(int i=0; i<trace_DLED_length; i++){
                trace_DLED[0][i] = trace[0][i+dleddt];
                trace_DLED[1][i] = trace[1][i+dleddt];
            }
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

                    if(jj==0){
                      min_line = trace_DLED[0][mintp];
                      max_line = trace_DLED[0][maxtp];
                    }


                    index_for_peak = FindPeaksFixTime(mintp, maxtp);

                    if(jj==0){
                      index_peak_LED = index_for_peak;
                    }

                    // peak always ok:
                    peak_rejected = false;

                      peak_LED[0] = trace_DLED[0][index_for_peak];
                      peak_LED[1] = trace_DLED[1][index_for_peak];
                      if(jj==0) {
                        // cout<<mintp<<"\t"<<maxtp<<endl;
                        // ptrHistLED->Fill(peak_LED[1]);
                      }
                      if(jj==1){
                        // cout<<mintp<<"\t"<<maxtp<<endl;
                        // ptrHistDCR_window->Fill(peak_LED[1]);
                      }


                } // END LOOP ON NUM WINDOWS

                // reuse the old values: (for LED)
                mintp = led_mintp;
                maxtp = led_maxtp;


        } // end FIND PEAKS LED


        /////////////////////////
        ///   FIND AREA LED   ///
        /////////////////////////
        if(find_area_trace_bool){
          // ptrHistArea->Fill(find_area_trace(peak_LED[0], peak_LED[0]-time_area_low, peak_LED[0]+time_area_high));
          // cout<<find_area_trace(peak_LED[0], peak_LED[0]-time_area_low, peak_LED[0]+time_area_high)<<endl;
        }



        /////////////////////
        /// PEAKS FINDING ///
        /////////////////////

        fill_hist_peaks_when_found = false;
        if(find_peaks_bool){ // find_peaks_bool
            // find_peaks(thr_to_find_peaks,max_peak_width, min_peak_width,blind_gap,DCR_DELAYS_bool);

            FindPeaksRisingFalling(thr_to_find_peaks, trace_DLED, trace_DLED_length, max_peak_width,2,4);

            if(n_ev==0){
                trace_time = 0;
                trace_time_raw = 0;
            }

            trace_time_raw += trace_DLED[0][trace_DLED_length-1] - trace_DLED[0][0];

            trace_time += ( trace_DLED[0][trace_DLED_length-1] - trace_DLED[0][0] - DCR_cnt_temp * 2 * dleddt - rise_time ); // time in trace is in ns

        } // END find_peaks_bool


        // AVERAGE
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


        // DISPLAY
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


}

//______________________________________________________________________________
void WriteIntro(){
    // write introduction in FilePeaks

    // find name
    std::size_t found = FilePeaksName.find_last_of("/\\");
    FilePeaks<<FilePeaksName.substr(found+1)<<endl;

    // print some important informations
    FilePeaks<<"Peaks found using a threshold = "<<thr_to_find_peaks<<" mV"<<endl;
    if(smooth_trace_bool) FilePeaks<<"YES ";
    else                  FilePeaks<<"NO ";
    FilePeaks<<"trace smoothing"<<endl;
    FilePeaks<<"TimeWindow = "<<1024<<" ns"<<endl;
    FilePeaks<<"DLED related: dleddt = "<<dleddt<<endl;
    FilePeaks<<"              blind_gap = "<<blind_gap<<endl;

    // END of INTRODUCTION
    FilePeaks<<"END_INTRODUCTION"<<endl;

}

//------------------------------------------------------------------------------
//---------------------------------[   HELP   ]---------------------------------
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
void help(){
    cout<<endl;
    cout<<endl;
}


# ifndef __CINT__
int main()
{
  cout<<"main"<<endl;
  return 0;
}
# endif
