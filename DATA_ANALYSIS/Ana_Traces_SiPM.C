//Ana_Traces_SiPM.C

/* Traces acquired by AGILENT MS06054A or by Digitizer CAEN DT 5751
 * 
 * Calculations:
 *      (i)  DCR            (Hamamatsu - MPPC Characterization)
 *      (ii) Cross Talk     (Hamamatsu - MPPC Characterization pag 44)
 * 
 * Tecniques:
 *      (i)  DLED for peak detections
 * 
 */


#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
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

//------------------------------------------------------------------------------
//-------------------------------[   FUNCTIONS   ]------------------------------
//------------------------------------------------------------------------------
int Analysis(string file, int last_event_n, bool display);
int loopAnalysis(string file1, int last_event_n);
int Ana3(string file1, string file2, string file3, int last_event_n);
int loopAna3(string file1, string file2, string file3, int last_event_n);


void DLED(int trace_lenght, int dleddt);
int find_peak_fix_time(int mintp, int maxtp);
int find_peaks(double thr_to_find_peaks, int max_peak_width, int min_peak_width, bool DCR_DELAYS_bool);
void average_func(int trace_lenght);
void show_trace(TCanvas* canv, double *x, double *y, int trace_lenght, double miny, double maxy, int mintp, int maxtp, bool line_bool, bool delete_bool, bool reverse);
void help();
void fit_hist_del(double expDelLow, double expDelHigh);
void ResetHistsDelays();
void Get_DCR_temp_and_errDCR_temp(int nfile);
void fit_hist_all_peaks(TCanvas *c, TH1D *hist, double fit1Low, double fit1High, double fit2Low, double fit2High);
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------


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

double DCR_temp[] = {0., 0., 0.}; //I consider 3 files
double errDCR_temp[] = {0., 0., 0.};
double DCR_pe_0_5_vect[] = {0., 0., 0.};
double DCR_pe_1_5_vect[] = {0., 0., 0.};
double errDCR_pe_0_5_vect[] = {0., 0., 0.};
double errDCR_pe_1_5_vect[] = {0., 0., 0.};

int trace_DLED_lenght;

int ii=0;
double max_func;
int index_func = 0;

double pe_0_5 = 10; // !!! CHECK !!! (see comments in function Ana3(...) for values)
double pe_1_5 = 25; // !!! CHECK !!! (see comments in function Ana3(...) for values)

double thr_to_find_peaks = 10; //thr_to_find_peaks, as seen in DLED trace (in V); it should be similar to pe_0_5

double w = 1000;
double h = 800;

int bins_Volt = 204;
int bins_DCR = 206;
int bins_Delays = 200;

double maxyhistAllPeaks = 200; 
double maxyhistDCR = 200;
double maxyhistDelays = 200;

double delta_pe = 0.005;

double expDelLow_max= 40.;
double expDelHigh_max = 160.;

int nfile = 1;

bool first_time_main_called = true;

//------------------------------------------------------------------------------
//-------------------------------[   OPTIONS   ]--------------------------------
//------------------------------------------------------------------------------

bool find_peak_in_a_selected_window = false; //to find a peak in a selected window (e.g. for LED measures)
bool average = false; //calculate the average of traces (useful for LED measures)

bool DCR_CNT_bool = false; // DCR from counters
bool DCR_DELAYS_bool = false; //DCR from delays
bool CROSS_TALK_bool = false; //DCR must be true

bool drawHistAllPeaks = true; // to draw hist of all peaks in traces
bool fitHistAllPeaks = true;
bool drawHistAllPeaksAll = false; // to draw hist of all peaks in traces for the 3 files (superimpose)
bool show_hists_DCR_DELAYS  = false;
bool showHist_bool = false; 
bool SetLogyHist = false;
bool running_graph = false;// to see traces in an osc mode (display must be true)

bool all_events_same_window = false; //all the events (from 0 to last_event_n) are joined in a single event

/* VALUES:
 * 
 * LED from Digitizer_CAEN: bins_Volt = 204;
 * 
 */
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
//--------------------[   GLOBAL HISTs, CANVs and FUNCs   ]---------------------
//------------------------------------------------------------------------------


TH1D *ptrHistAllPeaks1 = new TH1D("histAllPeaks1","",bins_DCR,0,maxyhistAllPeaks);
TH1D *ptrHistAllPeaks2 = new TH1D("histAllPeaks2","",bins_DCR,0,maxyhistAllPeaks);
TH1D *ptrHistAllPeaks3 = new TH1D("histAllPeaks3","",bins_DCR,0,maxyhistAllPeaks);

TH1D *ptrHistDCRthr1 = new TH1D("histDCRthr1","",bins_DCR,0,maxyhistDCR);
TH1D *ptrHistDCRthr2 = new TH1D("histDCRthr2","",bins_DCR,0,maxyhistDCR);
TH1D *ptrHistDCRthr3 = new TH1D("histDCRthr3","",bins_DCR,0,maxyhistDCR);

TH1D *ptrHistDelays_1 = new TH1D("histDelays_pe_0_5_1","",bins_Delays,0,maxyhistDelays);
TH1D *ptrHistDelays_2 = new TH1D("histDelays_pe_0_5_2","",bins_Delays,0,maxyhistDelays);
TH1D *ptrHistDelays_3 = new TH1D("histDelays_pe_0_5_3","",bins_Delays,0,maxyhistDelays);

TF1 *expDel = new TF1("expDel","TMath::Exp(-[0]*x+[1])",expDelLow_max,expDelHigh_max);
TF1 *gausFit1 = new TF1("gausFit1","gaus",-100,100);
TF1 *gausFit2 = new TF1("gausFit2","gaus",-100,100);

//I want to create ONLY one time the Canvas below... this is not the smartest way but it shuold work...
TCanvas *c = new TCanvas("Trace","Trace",w,h);
TCanvas *cDLED = new TCanvas("DLED","DLED",w,h);
TCanvas *cDCR = new TCanvas("hist_DCR","hist_DCR",w,h);
TCanvas *cAllPeaks = new TCanvas("AllPeaks","AllPeaks",w,h);

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------



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
    int dleddt = 9; //9ns is approx the rise time used for HD3_2 on AS out 2
    int max_peak_width = 12; //used for find_peaks
    int min_peak_width = 0;  //used for find_peaks
    double maxyhist = 200;
    
    double fit1Low = 0;
    double fit1High = 0;
    double fit2Low = 0;
    double fit2High = 0;
    
    //DEVICE:
    bool Agilent_MSO6054A = false; //true if data taken by Agilent_MSO6054A, false otherwise
    bool Digitizer_CAEN = true;  //true if data taken by Digitizer_CAEN, false otherwise
    
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

    if(first_time_main_called){
        if(DCR_CNT_bool and DCR_DELAYS_bool){
            cout<<"WARNING: both DCR_CNT_bool and DCR_DELAYS_bool are true: do you want to evaluate DCR both by counting peaks and fitting the delays distribution? Y/N";
            char opt = getchar();
            if(opt!='Y' and opt!='y')
                return 5;
        }
    }

   
    ifstream OpenFile (file.c_str());
    
    //Local variables
    char temp[20];
    int n_ev, index, i;
    bool reading = true;
    bool last_event_flag = false;
    
    int trace_lenght = 0;
    int one_window=0;
    int DCR_cnt = 0;
    
    TH1D *ptrHist = new TH1D("hist","",bins_Volt,0,maxyhist);
    
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
                if(all_events_same_window){
                    trace_lenght = atoi(temp)*last_event_n;
                    one_window = atoi(temp);
                }
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
        if(find_peak_in_a_selected_window){
            double *peak = new double[2];  
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
            index = find_peak_fix_time(mintp, maxtp);
            peak[0] = trace_DLED[0][index];
            peak[1] = trace_DLED[1][index];
            ptrHist->Fill(peak[1]);
        }
        
//***** DCR
        if(DCR_CNT_bool || DCR_DELAYS_bool){
            DCR_cnt = DCR_cnt + find_peaks(thr_to_find_peaks,max_peak_width, min_peak_width,DCR_DELAYS_bool);
            DCR_time = DCR_time + trace_lenght*TMath::Power(10,-9);
        }
        else{
            if(drawHistAllPeaks){//all peaks but not DCR
                thr_to_find_peaks = 10; //mV
                find_peaks(thr_to_find_peaks,max_peak_width, min_peak_width,DCR_DELAYS_bool);
            }
        }
        
        
        
//***** DISPLAY        
        if(display){
            if(n_ev==0){
                c->SetGrid();
                cDLED->SetGrid();
            }
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
    
//***** HIST ALL PEAKS
    if(drawHistAllPeaks){
        if(nfile == 1){
            TCanvas *cAllPeaks1 = new TCanvas("AllPeaks1","AllPeaks1",w,h);
            cAllPeaks1-> SetGrid();
            cAllPeaks1->cd();
            ptrHistAllPeaks1->Draw("hist");
            
            if(fitHistAllPeaks){
                if(file == "20180221_HD3-2_1_DARK_34_AS_2_01.txt" or file == "20180221_HD3-2_2_DARK_34_AS_2_02.txt" or file == "20180221_HD3-2_3_DARK_34_AS_2_01.txt")
                {fit1Low = 12; fit1High = 26; fit2Low = 28; fit2High = 42;}
                if(file == "20180221_HD3-2_1_DARK_35_AS_2_01.txt" or file == "20180221_HD3-2_2_DARK_35_AS_2_02.txt" or file == "20180221_HD3-2_3_DARK_35_AS_2_01.txt")
                {fit1Low = 12; fit1High = 28; fit2Low = 32; fit2High = 46;}
                if(file == "20180221_HD3-2_1_DARK_36_AS_2_01.txt" or file == "20180221_HD3-2_2_DARK_36_AS_2_02.txt" or file == "20180221_HD3-2_3_DARK_36_AS_2_01.txt")
                {fit1Low = 12; fit1High = 28; fit2Low = 36; fit2High = 50;}
                fit_hist_all_peaks(cAllPeaks1, ptrHistAllPeaks1, fit1Low, fit1High, fit2Low, fit2High);
            }
            
            if(drawHistAllPeaksAll){
                cAllPeaks-> SetGrid();
                cAllPeaks->cd();
                ptrHistAllPeaks1->Draw("histsame");
            }
            
        }
        if(nfile == 2){
            TCanvas *cAllPeaks2 = new TCanvas("AllPeaks2","AllPeaks2",w,h);
            cAllPeaks2-> SetGrid();
            cAllPeaks2->cd();
            ptrHistAllPeaks2->Draw("hist");
            
            if(drawHistAllPeaksAll){
                cAllPeaks->cd();
                ptrHistAllPeaks2->SetLineColor(kGreen+1);
                ptrHistAllPeaks2->Draw("histsame");
            }
        }
        if(nfile == 3){
            TCanvas *cAllPeaks3 = new TCanvas("AllPeaks3","AllPeaks3",w,h);
            cAllPeaks3-> SetGrid();
            cAllPeaks3->cd();
            ptrHistAllPeaks3->Draw("hist");
            
            if(drawHistAllPeaksAll){
                cAllPeaks->cd();
                ptrHistAllPeaks3->SetLineColor(kRed+1);
                ptrHistAllPeaks3->Draw("histsame");
                
                auto legend = new TLegend(0.7,0.7,0.9,0.9);
                legend->AddEntry(ptrHistDCRthr1,"HV = 34.00 V","l");
                legend->AddEntry(ptrHistDCRthr2,"HV = 35.00 V","l");
                legend->AddEntry(ptrHistDCRthr3,"HV = 36.00 V","l");
                legend->Draw();
            }
        }
    }
    
//***** DCR from CNT 
    if(DCR_CNT_bool){
        
        DCR = DCR_cnt/DCR_time;
        cout<<"\nDCR = "<<DCR*TMath::Power(10,-6)<<" MHz"<<endl;
        cDCR->SetGrid();
        cDCR->cd();
        
        double cnt_temp, binCenter;
        for(int i=0; i<bins_DCR; i++){
            cnt_temp=0;
            for(int j=i; j<bins_DCR; j++){
                if(nfile==1) cnt_temp = cnt_temp + ptrHistAllPeaks1->GetBinContent(j);
                if(nfile==2) cnt_temp = cnt_temp + ptrHistAllPeaks2->GetBinContent(j);
                if(nfile==3) cnt_temp = cnt_temp + ptrHistAllPeaks3->GetBinContent(j);
            }
            cnt_temp = cnt_temp/DCR_time;
            
            if(nfile==1){ binCenter = ptrHistAllPeaks1->GetBinCenter(i); ptrHistDCRthr1 -> Fill(binCenter,cnt_temp);}
            if(nfile==2){ binCenter = ptrHistAllPeaks2->GetBinCenter(i); ptrHistDCRthr2 -> Fill(binCenter,cnt_temp);}
            if(nfile==3){ binCenter = ptrHistAllPeaks3->GetBinCenter(i); ptrHistDCRthr3 -> Fill(binCenter,cnt_temp);}
        }
        cDCR->cd();
        
        cDCR->SetLogy();
        if(nfile == 1){
            ptrHistDCRthr1->GetXaxis()->SetTitle("Threshold (mV)");
            ptrHistDCRthr1->GetYaxis()->SetTitle("DCR (Hz)");
            ptrHistDCRthr1->Draw("hist");
        }
        if(nfile == 2){
            ptrHistDCRthr2->SetLineColor(kGreen+1);
            ptrHistDCRthr2->Draw("histsame");
        }
        if(nfile == 3){
            ptrHistDCRthr3->SetLineColor(kRed+1);
            ptrHistDCRthr3->Draw("histsame");
            auto legend = new TLegend(0.7,0.7,0.9,0.9);
            legend->AddEntry(ptrHistDCRthr1,"HV = 34.00 V","l");
            legend->AddEntry(ptrHistDCRthr2,"HV = 35.00 V","l");
            legend->AddEntry(ptrHistDCRthr3,"HV = 36.00 V","l");
            legend->Draw();
        }

//***** CROSS-TALK
    if(CROSS_TALK_bool){
        double DCR_pe_0_5 = 0.;
        double DCR_pe_1_5 = 0.;
        double CrossTalk = 0.;
        
        if(pe_0_5==0||pe_1_5==0){
            cout<<"ERROR: CHECK THR FOR 0.5 AND 1.5 PE"<<endl;
            return 1;
        }
        
        for(int i=0; i<bins_DCR; i++){
            if(nfile==1){
                if((ptrHistAllPeaks1->GetBinCenter(i) > pe_0_5 - delta_pe) and (ptrHistAllPeaks1->GetBinCenter(i) < pe_0_5 + delta_pe)){
                    DCR_pe_0_5 = ptrHistDCRthr1->GetBinContent(i);
                    cout<<"Found DCR_pe_0_5"<<endl;
                }
                if((ptrHistAllPeaks1->GetBinCenter(i) > pe_1_5 - delta_pe) and (ptrHistAllPeaks1->GetBinCenter(i) < pe_1_5 + delta_pe)){
                    DCR_pe_1_5 = ptrHistDCRthr1->GetBinContent(i);
                    cout<<"Found DCR_pe_1_5"<<endl;
                }
            }
            if(nfile==2){
                if((ptrHistAllPeaks2->GetBinCenter(i) > pe_0_5 - delta_pe) and (ptrHistAllPeaks2->GetBinCenter(i) < pe_0_5 + delta_pe)){
                    DCR_pe_0_5 = ptrHistDCRthr2->GetBinContent(i);
                    cout<<"Found DCR_pe_0_5"<<endl;
                }
                if((ptrHistAllPeaks2->GetBinCenter(i) > pe_1_5 - delta_pe) and (ptrHistAllPeaks2->GetBinCenter(i) < pe_1_5 + delta_pe)){
                    DCR_pe_1_5 = ptrHistDCRthr2->GetBinContent(i);
                    cout<<"Found DCR_pe_1_5"<<endl;
                }
            }
            if(nfile==3){
                if((ptrHistAllPeaks3->GetBinCenter(i) > pe_0_5 - delta_pe) and (ptrHistAllPeaks3->GetBinCenter(i) < pe_0_5 + delta_pe)){
                    DCR_pe_0_5 = ptrHistDCRthr3->GetBinContent(i);
                    cout<<"Found DCR_pe_0_5"<<endl;
                }
                if((ptrHistAllPeaks3->GetBinCenter(i) > pe_1_5 - delta_pe) and (ptrHistAllPeaks3->GetBinCenter(i) < pe_1_5 + delta_pe)){
                    DCR_pe_1_5 = ptrHistDCRthr3->GetBinContent(i);
                    cout<<"Found DCR_pe_1_5"<<endl;
                }
            }
        }
        
        CrossTalk = DCR_pe_1_5/DCR_pe_0_5;
        cout<<"CrossTalk = "<<CrossTalk<<endl;
}

    }
    
//***** DCR from DELAYS
    if(DCR_DELAYS_bool){
        if(show_hists_DCR_DELAYS){
            if(nfile == 1){
                    TCanvas *cDelays1 = new TCanvas("Delays_SiPM_file1","Delays_SiPM_file1",w,h);
                    ptrHistDelays_1->GetXaxis()->SetTitle("Time (ns)");
                    ptrHistDelays_1->GetYaxis()->SetTitle("");
                    ptrHistDelays_1->Draw("hist");
            }
            if(nfile == 2){
                    TCanvas *cDelays2 = new TCanvas("Delays_SiPM_file2","Delays_SiPM_file2",w,h);
                    ptrHistDelays_2->GetXaxis()->SetTitle("Time (ns)");
                    ptrHistDelays_2->GetYaxis()->SetTitle("");
                    ptrHistDelays_2->Draw("hist");
            }
            if(nfile == 3){
                    TCanvas *cDelays3 = new TCanvas("Delays_SiPM_file3","Delays_SiPM_file3",w,h);
                    ptrHistDelays_3->GetXaxis()->SetTitle("Time (ns)");
                    ptrHistDelays_3->GetYaxis()->SetTitle("");
                    ptrHistDelays_3->Draw("hist");
            }
        }
        
        fit_hist_del(expDelLow_max, expDelHigh_max);
        
    }
    
    if(showHist_bool){
        TCanvas *cHist = new TCanvas("hist_GAIN","hist_GAIN",w,h);
        cHist->SetGrid();
        cHist->cd();
        if(SetLogyHist) cHist->SetLogy();
        ptrHist->Draw("hist");
    }
    
    delete ptrHist;
    
    first_time_main_called = false;
    
    return 0;
}



//------------------------------------------------------------------------------
//---------------------------[   OTHER FUNCTIONS   ]----------------------------
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
void help(){
    cout<<"USAGE:"<<endl;
    cout<<"Analysis(string file, int last_event_n, bool display)"<<endl;
    cout<<"Ana3(string file1, string file2, string file3, int last_event_n)"<<endl;
}

//------------------------------------------------------------------------------
int loopAnalysis(string file1, int last_event_n){
    bool display = false;
    int control = 0;
    
    pe_0_5 = 10; //mV
    pe_1_5 = 25; //mV

    
    thr_to_find_peaks = 7; //mV
    double max_thr_to_find_peaks = 40; //mV
    double gap_between_thr = 1; //mV
    
    int n_DCR;
    n_DCR = (int)((max_thr_to_find_peaks - thr_to_find_peaks)/gap_between_thr);
    
    DCR = new double*[1]; //first index: nfile => DCR[0][...] for nfile = 1 (in this case I have only one file, not as in loopAna3)
        for(int i = 0; i < 3; i++) {
            DCR[i] = new double[n_DCR];
    }
    errDCR = new double*[1]; //first index: nfile => DCR[0][...] for nfile = 1 (in this case I have only one file, not as in loopAna3)
        for(int i = 0; i < 3; i++) {
            errDCR[i] = new double[n_DCR];
    }
    
    int h = 0;
    double *DCR_thr = new double[n_DCR]; 
    
    while(thr_to_find_peaks <= max_thr_to_find_peaks){
        control = Analysis(file1,last_event_n,display);
        if(control==5)
            return 5;
        Get_DCR_temp_and_errDCR_temp(nfile);
        ResetHistsDelays();
        DCR_thr[h] = thr_to_find_peaks; //mV
        DCR[0][h] = DCR_temp[0];
        errDCR[0][h] = errDCR_temp[0];
        thr_to_find_peaks = thr_to_find_peaks + gap_between_thr; //I jump to the next thr in order to evaluate the new DRC
        h++;
    } 
    
    TGraphErrors *gDCR_1 = new TGraphErrors(n_DCR, DCR_thr, DCR[0],NULL, errDCR[0]);
    
    gDCR_1->SetLineColor(kGreen+1);
    
    gDCR_1->SetLineWidth(2);
    
    gDCR_1->SetFillColorAlpha(kGreen+1, 0.3);
    
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

    double CrossTalk = 0.;
    double errCrossTalk = 0.;
    for(int i=0; i<1; i++){
        cout<<"FILE "<<i<<endl;
        cout<<"   DCR at 0.5 pe = ("<<DCR_pe_0_5_vect[i]*TMath::Power(10,-6)<<" +- "<<errDCR_pe_0_5_vect[i]*TMath::Power(10,-6)<<") MHz"<<endl;
        cout<<"   DCR at 1.5 pe = ("<<DCR_pe_1_5_vect[i]*TMath::Power(10,-6)<<" +- "<<errDCR_pe_1_5_vect[i]*TMath::Power(10,-6)<<") MHz"<<endl;
        if(CROSS_TALK_bool){
            CrossTalk = DCR_pe_1_5_vect[i]/DCR_pe_0_5_vect[i];
            errCrossTalk = CrossTalk * TMath::Sqrt( (errDCR_pe_0_5_vect[i]/DCR_pe_0_5_vect[i])*(errDCR_pe_0_5_vect[i]/DCR_pe_0_5_vect[i]) + (errDCR_pe_1_5_vect[i]/DCR_pe_1_5_vect[i])*(errDCR_pe_1_5_vect[i]/DCR_pe_1_5_vect[i]) );
        cout<<"   Cross Talk    = ("<<CrossTalk<<" +- "<<errCrossTalk<<")"<<endl;
        }
        cout<<endl;
    }
    
    return 0;
    
}

//------------------------------------------------------------------------------
int Ana3(string file1, string file2, string file3, int last_event_n){
    
    bool display = false;
    int control = 0;
    //file1
    nfile = 1;
    pe_0_5 = 10;   pe_1_5 = 24;
    control = Analysis(file1,last_event_n,display);
    
    if(control==5)
        return 5;
    
    Get_DCR_temp_and_errDCR_temp(nfile);
    ResetHistsDelays();
    
    //file2
    nfile=2;
    pe_0_5 = 10;   pe_1_5 = 26;
    control = Analysis(file2,last_event_n,display);
    
    Get_DCR_temp_and_errDCR_temp(nfile);
    ResetHistsDelays();
    
    //file3
    nfile=3;
    pe_0_5 = 10;   pe_1_5 = 26;
    control = Analysis(file3,last_event_n,display);
    
    Get_DCR_temp_and_errDCR_temp(nfile);
    ResetHistsDelays();

    
/*          HV      0.5pe   1.5pe   FILE
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

    
    return 0;
}

//------------------------------------------------------------------------------
int loopAna3(string file1, string file2, string file3, int last_event_n){ //in order to draw the stair plot for DCR; DCR is obtained by exp Fit of peaks delay distribution
    bool display = false;
    int control = 0;
    
    thr_to_find_peaks = 7; //mV
    double max_thr_to_find_peaks = 41; //mV
    double gap_between_thr = 1; //mV
    
    int n_DCR;
    n_DCR = (int)((max_thr_to_find_peaks - thr_to_find_peaks)/gap_between_thr);
    
    DCR = new double*[3]; //first index: nfile => DCR[0][...] for nfile = 1; DCR[1][...] for nfile = 2; DCR[2][...] for nfile = 3. 
        for(int i = 0; i < 3; i++) {
            DCR[i] = new double[n_DCR];
    }
    errDCR = new double*[3]; //first index: nfile => errDCR[0][...] for nfile = 1; errDCR[1][...] for nfile = 2; errDCR[2][...] for nfile = 3. 
        for(int i = 0; i < 3; i++) {
            errDCR[i] = new double[n_DCR];
    }
    
    int h = 0;
    double *DCR_thr = new double[n_DCR]; 
    
    while(thr_to_find_peaks <= max_thr_to_find_peaks){
        control = Ana3(file1, file2, file3,last_event_n);
        if(control==5)
            return 5;
        DCR_thr[h] = thr_to_find_peaks; //mV
        DCR[0][h] = DCR_temp[0];
        DCR[1][h] = DCR_temp[1];
        DCR[2][h] = DCR_temp[2];
        errDCR[0][h] = errDCR_temp[0];
        errDCR[1][h] = errDCR_temp[1];
        errDCR[2][h] = errDCR_temp[2];
        thr_to_find_peaks = thr_to_find_peaks + gap_between_thr; //I jump to the next thr in order to evaluate the new DRC
        h++;
    } 
    
    TGraphErrors *gDCR_1 = new TGraphErrors(n_DCR, DCR_thr, DCR[0],NULL, errDCR[0]);
    TGraphErrors *gDCR_2 = new TGraphErrors(n_DCR, DCR_thr, DCR[1],NULL, errDCR[1]);
    TGraphErrors *gDCR_3 = new TGraphErrors(n_DCR, DCR_thr, DCR[2],NULL, errDCR[2]);
    
    gDCR_1->SetLineColor(kGreen+1);
    gDCR_2->SetLineColor(kBlue);
    gDCR_3->SetLineColor(kRed+1);
    
    gDCR_1->SetLineWidth(2);
    gDCR_2->SetLineWidth(2);
    gDCR_3->SetLineWidth(2);
    
    gDCR_1->SetFillColorAlpha(kGreen+1, 0.3);
    gDCR_2->SetFillColorAlpha(kBlue, 0.3);
    gDCR_3->SetFillColorAlpha(kRed+1, 0.3);
    
    
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

    double CrossTalk = 0.;
    double errCrossTalk = 0.;
    for(int i=0; i<3; i++){
        cout<<"FILE "<<i<<endl;
        cout<<"   DCR at 0.5 pe = ("<<DCR_pe_0_5_vect[i]*TMath::Power(10,-6)<<" +- "<<errDCR_pe_0_5_vect[i]*TMath::Power(10,-6)<<") MHz"<<endl;
        cout<<"   DCR at 1.5 pe = ("<<DCR_pe_1_5_vect[i]*TMath::Power(10,-6)<<" +- "<<errDCR_pe_1_5_vect[i]*TMath::Power(10,-6)<<") MHz"<<endl;
        if(CROSS_TALK_bool){
            CrossTalk = DCR_pe_1_5_vect[i]/DCR_pe_0_5_vect[i];
            errCrossTalk = CrossTalk * TMath::Sqrt( (errDCR_pe_0_5_vect[i]/DCR_pe_0_5_vect[i])*(errDCR_pe_0_5_vect[i]/DCR_pe_0_5_vect[i]) + (errDCR_pe_1_5_vect[i]/DCR_pe_1_5_vect[i])*(errDCR_pe_1_5_vect[i]/DCR_pe_1_5_vect[i]) );
        cout<<"   Cross Talk    = ("<<CrossTalk<<" +- "<<errCrossTalk<<")"<<endl;
        }
        cout<<endl;
    }
    
    
    return 0;
    
}


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
int find_peaks(double thr_to_find_peaks, int max_peak_width, int min_peak_width, bool DCR_DELAYS_bool){ //I look for every peaks in the trace, only if I'm in DARK mode
    ii=0;
    int index_peak;
    int DCR_cnt_temp = 0;
    int index_old = 0;
    int index_new = 0;
    bool found_first_peak = false;
    int peak_width = max_peak_width;
    while(ii<trace_DLED_lenght){//I find peaks after DLED
        if(trace_DLED[1][ii]>thr_to_find_peaks){//I only consider points above thr_to_find_peaks
            DCR_cnt_temp++; //I've seen a peak; if I'm in dark mode it's DCR
            
            //Now I want to see the peak amplitude.
            for(int k=ii+min_peak_width; (k<ii+max_peak_width); k++){
                if(trace_DLED[1][k]<thr_to_find_peaks){
                    peak_width = k-ii;
                    break;
                }
            }
            
            if(ii+peak_width<trace_DLED_lenght)
                index_peak = find_peak_fix_time(ii, ii+peak_width);
            else 
                index_peak = find_peak_fix_time(ii, trace_DLED_lenght);
            
            //DCR from the delay (Itzler Mark - Dark Count Rate Measure (pag 5 ss))
            if(DCR_DELAYS_bool){
                if(!found_first_peak){
                    index_old = index_peak;
                    found_first_peak = true;
                }
                else{
                    index_new=index_peak;
                    //I fill the hist with delays between 1 pe peaks (or higher):
                    if(nfile == 1) ptrHistDelays_1 -> Fill(index_new - index_old); 
                    if(nfile == 2) ptrHistDelays_2 -> Fill(index_new - index_old);
                    if(nfile == 3) ptrHistDelays_3 -> Fill(index_new - index_old);  
                    index_old = index_new;
                    
                }
            }
            
            if(nfile == 1) ptrHistAllPeaks1->Fill(trace_DLED[1][index_peak]);
            if(nfile == 2) ptrHistAllPeaks2->Fill(trace_DLED[1][index_peak]);
            if(nfile == 3) ptrHistAllPeaks3->Fill(trace_DLED[1][index_peak]);
            
            ii=ii+peak_width; //in order to jump at the following peak.
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
        y[ii] = y[ii];
        
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

//------------------------------------------------------------------------------
void ResetHistsDelays(){
    ptrHistDelays_1->Reset();
    ptrHistDelays_2->Reset();
    ptrHistDelays_3->Reset();
}

//------------------------------------------------------------------------------
void Get_DCR_temp_and_errDCR_temp(int nfile){
    DCR_temp[nfile-1] = expDel->GetParameter(0)*TMath::Power(10,9);
    errDCR_temp[nfile-1] = expDel->GetParError(0)*TMath::Power(10,9);
    
    if(CROSS_TALK_bool){
        if((int)(thr_to_find_peaks*10000)==(int)(pe_0_5*10000)){
            DCR_pe_0_5_vect[nfile-1] = DCR_temp[nfile-1];
            errDCR_pe_0_5_vect[nfile-1] = errDCR_temp[nfile-1];
        }
        else{
            if((int)(thr_to_find_peaks*10000)==(int)(pe_1_5*10000)){
                DCR_pe_1_5_vect[nfile-1] = DCR_temp[nfile-1];
                errDCR_pe_1_5_vect[nfile-1] = errDCR_temp[nfile-1];
            }
        }
    }
    
}


//------------------------------------------------------------------------------
void fit_hist_del(double expDelLow, double expDelHigh){ //fit hists filled with time delays in order to find DCR
    cout<<"Fit hist delays file "<<nfile<<endl;
    if(nfile == 1) ptrHistDelays_1 -> Fit(expDel, "q", "", expDelLow, expDelHigh);
    if(nfile == 2) ptrHistDelays_2 -> Fit(expDel, "q", "", expDelLow, expDelHigh);
    if(nfile == 3) ptrHistDelays_3 -> Fit(expDel, "q", "", expDelLow, expDelHigh);
}


//------------------------------------------------------------------------------
void fit_hist_all_peaks(TCanvas *c, TH1D *hist, double fit1Low, double fit1High, double fit2Low, double fit2High){
    
    double mean1, mean2, gain;
    double errmean1, errmean2, errgain;
    
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

