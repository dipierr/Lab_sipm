//reads the output files on CAEN digitizer and fills a root tree with peak for the 1 adc channel.
// The tree is saved in a root file called out.root


#include <cstdio>
#include <cstdlib>
#include <string>
#include <map>
#include <TH1.h>
#include <TH1F.h>
#include <TH2.h>
#include <TF1.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TMath.h>
#include <TFile.h>
#include <TTree.h>

#include <TStyle.h>
#include <TCanvas.h>
#include <TApplication.h>
#include <TROOT.h>
#include <TMinuit.h>

#include <TLegend.h>
#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>

using namespace std;

/*
Bin interval for baseline calculation, set the CAEN digitizer horizontal offset accordingly or viceversa
*/

void ana_trace(ifstream &fileIn, double & peak_adc, int & bin_peak, double & charge, double & baseline, int baseline_start, int baseline_nbins, int search_peak_window_start, int search_peak_window_nbins, int sum_before_peak, int sum_bins){
   
	char dummy[30];
	int board_id, n_samples, channel,event_n;
  
	fileIn>>dummy>>dummy>>n_samples; 
	fileIn>>dummy>>board_id;
	fileIn>>dummy>>channel;
	fileIn>>dummy>>dummy>>event_n;
	fileIn>>dummy>>dummy;
	fileIn>>dummy>>dummy>>dummy>>dummy;
	fileIn>>dummy>>dummy>>dummy>>dummy;

	baseline = 0.;
	double n_baseline=0.;
 	
	peak_adc=0.;
	
	double data_adc;
	double trace[n_samples];
	double timebin = 1.; //ns
	double oscilloscope_res = 50.; //ohm
	
	if(!fileIn.eof()){    	
	
		for (int j=0; j<n_samples; j++){
			fileIn >> data_adc;
			trace[j] = data_adc;
			if (j >= baseline_start && j < baseline_start + baseline_nbins ){
				baseline = baseline + data_adc;
				n_baseline++;
			} 
				
			if (j >= search_peak_window_start && j < search_peak_window_start + search_peak_window_nbins ){
				if (data_adc > peak_adc){
					peak_adc = data_adc;
					bin_peak = j;
				}
			}
			
		}
		baseline = baseline/n_baseline;
	   peak_adc = peak_adc - baseline;
		charge = 0.;
		for (int j = bin_peak - sum_before_peak; j < bin_peak - sum_before_peak + sum_bins; j++){		   
			charge = charge + (trace[j]-baseline)*timebin;
		}
		charge = charge / oscilloscope_res; //1e-12 C, pC
	}
}


int main(int argc, char **argv){
  
	string configuration_file; 

	if(argc<2){
		cerr<<"ERRORE : Numero di Argomenti errato "<<endl; 
		cerr<<"Comando : ./spettri_all.exe [file_di_configurazione.cfg] "<<endl;
		cerr<<"Non si e' trovato il file di configurazione: se ne usa uno di default : spettri.cfg"<<endl;
		configuration_file="spettri.cfg";        
	}
	else{
		configuration_file=argv[1];    
	}
  
	ifstream spettri_config_file(configuration_file.c_str());
	if (!spettri_config_file.is_open()){
		cerr<<"ERRORE : Non si riesce a torvare il file di configurazione : "<<configuration_file<<endl;
		exit(EXIT_FAILURE);
	}  
  	
	char dumm[30];
	string filename, filenameout;
	int baseline_nbins,baseline_start;
	int search_peak_window_nbins,search_peak_window_start;
	int sum_before_peak,sum_bins;
	
	spettri_config_file>>dumm>>filename;
	spettri_config_file>>dumm>>filenameout;
	spettri_config_file>>dumm>>baseline_nbins;
	spettri_config_file>>dumm>>baseline_start;
	spettri_config_file>>dumm>>search_peak_window_nbins;
	spettri_config_file>>dumm>>search_peak_window_start;
	spettri_config_file>>dumm>>sum_before_peak;
	spettri_config_file>>dumm>>sum_bins;
	
	cout<<"Filename in: "<<filename<<endl;
	cout<<"Filename out: "<<filenameout<<endl;
	cout<<"Number of bins for baseline calculation: "<<baseline_nbins<<endl;
	cout<<"Starting bin for baseline calculation: "<<baseline_start<<endl;
	cout<<"Number of bins of the window for peak search: "<<search_peak_window_nbins<<endl;
	cout<<"Starting bin for peak search: "<<search_peak_window_start<<endl;
	cout<<"Number of bins before peak included in summation: "<<sum_before_peak<<endl;
	cout<<"Number of bins included in summation around peak: "<<sum_bins<<endl;
	
	double baseline;
	double ch_mv = 1000./1024.;
	double peak_adc;
	double peak_mV;
	int bin_peak;
	double charge;
	int counter = 0 ;
	
 	//open output root file//
	//TFile *rootOutFile = new TFile(filenameout.c_str(),"recreate");
	TFile *rootOutFile = new TFile("sipm_spectrum.root","recreate");
	//create the tree
 	char tname[7] = "traces";
 	TTree*  tree1 = new TTree(tname,tname);
	//make the branches//
	tree1->Branch("Event_id", &counter, "Event ID/i");	
	tree1->Branch("Baseline", &baseline, "baseline/D");
	tree1->Branch("PeakADC", &peak_adc, "peakadc/D");
	tree1->Branch("PeakmV", &peak_mV, "peakmV/D");
	tree1->Branch("Bin_Peak", &bin_peak, "Bin Peak/i");
	tree1->Branch("Charge", &charge, "Charge/D");
	
	TH1F *h1peakspectrum = new TH1F("Peak_Spectrum","Peak spectrum",300,0.,30.);
	h1peakspectrum->GetXaxis()->SetTitle("mV");
	
	TH1F *h1chargespectrum = new TH1F("Charge_Spectrum","Charge spectrum",300,0.,10.);
	h1chargespectrum->GetXaxis()->SetTitle("pC");
	
	ifstream fileIn0 (filename.c_str());
	
	while (!fileIn0.eof()){ // loop lettura file 
		
		ana_trace(fileIn0, peak_adc, bin_peak, charge, baseline, baseline_start, baseline_nbins, search_peak_window_start, search_peak_window_nbins, sum_before_peak, sum_bins );
	   
		peak_mV = peak_adc * ch_mv;
		
	   if (fileIn0.eof()) break;
	   
		tree1->Fill();
		
		h1peakspectrum->Fill(peak_mV);
		h1chargespectrum->Fill(charge);
		counter++;
		  
      cout<<"Events: "<<counter<<endl;	  
	} // end loop lettura eventi

h1peakspectrum->Write();
h1chargespectrum->Write();            
tree1->Write();
rootOutFile->Close();

}
