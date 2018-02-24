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
#include "TAxis.h"
#include "TPad.h"
#include "TLine.h"
#include "TSpectrum.h"
#include "TVirtualFitter.h"

//------------------------------------------------------------------------------
//-------------------------------[   FUNCTIONS   ]------------------------------
//------------------------------------------------------------------------------
void DLED(int trace_lenght, int dleddt);
int find_peak_fix_time(int mintp, int maxtp);
int find_peaks(double thr_to_find_peaks, int max_peak_width, int min_peak_width, bool delay_bool);
void average_func(int trace_lenght);
void show_trace(TCanvas* canv, double *x, double *y, int trace_lenght, double miny, double maxy, int mintp, int maxtp, bool line_bool, bool delete_bool, bool reverse);
void help();
int Ana3(string file1, string file2, string file3, int last_event_n, bool display);
void fit_hist_del();
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

double pe_0_5 = 0.010; // !!! CHECK !!! (see comments in function Ana3(...) for values)
double pe_1_5 = 0.025; // !!! CHECK !!! (see comments in function Ana3(...) for values)


double w = 1000;
double h = 800;

int bins_Volt = 204;
int bins_DCR = 206;
int bins_Delays = 200;

double maxyhistAllPeaks = .2; 
double maxyhistDCR = 200;
double maxyhistDelays = 200;

double delta_pe = 0.0005;

double expDelLow= 30.;
double expDelHigh = 160.;

int nfile = 1;

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

TH1D *ptrHistDelays_pe_0_5_1 = new TH1D("histDelays_pe_0_5_1","",bins_Delays,0,maxyhistDelays);
TH1D *ptrHistDelays_pe_0_5_2 = new TH1D("histDelays_pe_0_5_2","",bins_Delays,0,maxyhistDelays);
TH1D *ptrHistDelays_pe_0_5_3 = new TH1D("histDelays_pe_0_5_3","",bins_Delays,0,maxyhistDelays);

TF1 *expDel = new TF1("expDel","[1]*[0]*TMath::Exp(-[0]*x)",expDelLow,expDelHigh);

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
    double thr_to_find_peaks = 0.010; //thr_to_find_peaks, as seen in DLED trace (in V); it should be similar to pe_0_5
    int max_peak_width = 12; //used for find_peaks
    int min_peak_width = 5;
    double maxyhist = .2;
   
    
    bool DCR_bool = true;
    bool CROSS_TALK_bool = true; //if true, also DCR_bool must be true
    bool delay_bool = true; //if true, also DCR_bool and CROSS_TALK_bool must be true
    bool average = false;
    bool drawHistAllPeaks = true;
    bool drawHistAllPeaksAll = true;
    
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
    cDCR->SetGrid();
    
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
            DCR_cnt = DCR_cnt + find_peaks(thr_to_find_peaks,max_peak_width, min_peak_width,delay_bool);
            DCR_time = DCR_time + trace_lenght*TMath::Power(10,-9);
            //cout<<find_peaks(trace_lenght,thr_to_find_peaks,max_peak_width)<<endl;
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
        
        if(drawHistAllPeaks){
            if(nfile == 1){
                TCanvas *cAllPeaks1 = new TCanvas("AllPeaks1","AllPeaks1",w,h);
                cAllPeaks1-> SetGrid();
                cAllPeaks1->cd();
                ptrHistAllPeaks1->Draw("hist");
                
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
        
        
        DCR = DCR_cnt/DCR_time;
        cout<<"\nDCR = "<<DCR*TMath::Power(10,-6)<<" MHz"<<endl;
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
            
            if(nfile==1){ binCenter = ptrHistAllPeaks1->GetBinCenter(i); ptrHistDCRthr1 -> Fill(binCenter*1000,cnt_temp);}
            if(nfile==2){ binCenter = ptrHistAllPeaks2->GetBinCenter(i); ptrHistDCRthr2 -> Fill(binCenter*1000,cnt_temp);}
            if(nfile==3){ binCenter = ptrHistAllPeaks3->GetBinCenter(i); ptrHistDCRthr3 -> Fill(binCenter*1000,cnt_temp);}
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
    if(CROSS_TALK_bool and DCR_bool){
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
        
        //Calucalte delays if I want DCR and CrossTalk
        if(delay_bool){
            if(nfile == 1){
                TCanvas *cDelays = new TCanvas("Delays_SiPM1","Delays_SiPM1",w,h);
                ptrHistDelays_pe_0_5_1->GetXaxis()->SetTitle("Time (ns)");
                ptrHistDelays_pe_0_5_1->GetYaxis()->SetTitle("");
                //fit_hist_del();
                ptrHistDelays_pe_0_5_1->Draw("hist");
            }
            if(nfile == 2){
                TCanvas *cDelays = new TCanvas("Delays_SiPM2","Delays_SiPM2",w,h);
                ptrHistDelays_pe_0_5_2->GetXaxis()->SetTitle("Time (ns)");
                ptrHistDelays_pe_0_5_2->GetYaxis()->SetTitle("");
                ptrHistDelays_pe_0_5_2->Draw("hist");
            }
            if(nfile == 3){
                TCanvas *cDelays = new TCanvas("Delays_SiPM3","Delays_SiPM3",w,h);
                ptrHistDelays_pe_0_5_3->GetXaxis()->SetTitle("Time (ns)");
                ptrHistDelays_pe_0_5_3->GetYaxis()->SetTitle("");
                ptrHistDelays_pe_0_5_3->Draw("hist");
            }
        }
        
    }
        
        
    }
    
    cHist->cd();
    if(SetLogyHist) cHist->SetLogy();
    ptrHist->Draw();
    
        
    return 0;
}



//------------------------------------------------------------------------------
//---------------------------[   OTHER FUNCTIONS   ]----------------------------
//------------------------------------------------------------------------------

void help(){
    cout<<"USAGE:"<<endl;
    cout<<"Analysis(string file, int last_event_n, bool display)"<<endl;
    cout<<"Ana3(string file1, string file2, string file3, int last_event_n, bool display)"<<endl;
}



//------------------------------------------------------------------------------
int Ana3(string file1, string file2, string file3, int last_event_n, bool display){
    //file1
    nfile = 1;
    pe_0_5 = 0.010;   pe_1_5 = 0.025;
    Analysis(file1,last_event_n,display);
    
    //file2
    nfile=2;
    pe_0_5 = 0.010;   pe_1_5 = 0.026;
    Analysis(file2,last_event_n,display);
    
    //file3
    nfile=3;
    pe_0_5 = 0.010;   pe_1_5 = 0.028;
    Analysis(file3,last_event_n,display);
    
/*          HV      0.5pe   1.5pe   FILE
    --------------------------------------------------------------------
    SiPM1   34 V    10 mV   26 mV   20180221_HD3-2_1_DARK_34_AS_2_01.txt
            35 V    10 mV   28 mV   20180221_HD3-2_1_DARK_35_AS_2_01.txt
            36 V    10 mV   30 mV   20180221_HD3-2_1_DARK_36_AS_2_01.txt
    --------------------------------------------------------------------
    SiPM2   34 V    10 mV   25 mV   20180221_HD3-2_2_DARK_34_AS_2_02.txt
            35 V    10 mV   26 mV   20180221_HD3-2_2_DARK_34_AS_2_02.txt
            36 V    10 mV   28 mV   20180221_HD3-2_2_DARK_34_AS_2_02.txt
    --------------------------------------------------------------------
    SiPM3   34 V    10 mV   25 mV   20180221_HD3-2_3_DARK_34_AS_2_01.txt
            35 V    10 mV   26 mV   20180221_HD3-2_3_DARK_35_AS_2_01.txt
            36 V    10 mV   28 mV   20180221_HD3-2_3_DARK_36_AS_2_01.txt
*/
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
int find_peaks(double thr_to_find_peaks, int max_peak_width, int min_peak_width, bool delay_bool){ //I look for every peaks in the trace, only if I'm in DARK mode
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
            if(delay_bool){
                if(!found_first_peak){
                    index_old = index_peak;
                    found_first_peak = true;
                }
                else{
                    index_new=index_peak;
                    //I fill the hist with delays between 1 pe peaks (or higher):
                    if(nfile == 1) ptrHistDelays_pe_0_5_1 -> Fill(index_new - index_old); 
                    if(nfile == 2) ptrHistDelays_pe_0_5_2 -> Fill(index_new - index_old);
                    if(nfile == 3) ptrHistDelays_pe_0_5_3 -> Fill(index_new - index_old);  
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
void fit_hist_del(){ //fit hists filled with time delays in order to find DCR
    cout<<"Fit hist delays file "<<nfile<<endl;
    if(nfile == 1) ptrHistDelays_pe_0_5_1 -> Fit(expDel, "+", "", expDelLow, expDelHigh);
    if(nfile == 2) ptrHistDelays_pe_0_5_2 -> Fit(expDel, "+", "", expDelLow, expDelHigh);
    if(nfile == 3) ptrHistDelays_pe_0_5_3 -> Fit(expDel, "+", "", expDelLow, expDelHigh);
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

