#include "AnaSipm.h"
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <algorithm>
#include <functional>
#include "TFile.h"
#include "TTree.h"
#include "TTask.h"
#include "TVirtualGraphPainter.h"
#include "TGraph.h"
#include "TSpectrum.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TGaxis.h"
#include "TROOT.h"
#include "TString.h"

ClassImp(Trace);
ClassImp(TraceHeader);
ClassImp(RunSipm);
ClassImp(Decode);
ClassImp(CharacterizeSipm);

//______________________________________________________________________________
Trace::Trace() : fId(0), fTrace_length(1024), fAmplitude_Array(fTrace_length, 0), fTime_Array(fTrace_length, 0)
{
   // Default constructor with only default values.
}

//______________________________________________________________________________
Trace::~Trace()
{

}

//______________________________________________________________________________
void Trace::Build(UInt_t id, Float_t *amplitude_array, Float_t *time_array)
{
   // Fill the amplitude and time vectors starting from
   // amplitude_array and time_array.

   fId = id;

   fAmplitude_Array.clear();
   fTime_Array.clear();

   // fill vectors with assign() starting from C-style array instead of 
   // using a for loop. Uses iterators (better: pointers as iterators)
   // from: https://stackoverflow.com/questions/2434196/how-to-initialize-stdvector-from-c-style-array

   fAmplitude_Array.assign(amplitude_array, amplitude_array + fTrace_length);
   fTime_Array.assign(time_array, time_array + fTrace_length);

}

//______________________________________________________________________________
void Trace::Negate()
{
   // Reverse the amplitude trace, if any
   std::transform(fAmplitude_Array.begin(), fAmplitude_Array.end(), fAmplitude_Array.begin(), std::negate<Float_t>());
}

//______________________________________________________________________________
std::vector<Float_t> Trace::GetAmplitudeArray()
{
   return fAmplitude_Array;
}

//______________________________________________________________________________
std::vector<Float_t> Trace::GetTimeArray()
{
   return fTime_Array;
}

//______________________________________________________________________________
std::vector<Float_t> Trace::AmplitudeDLED(Int_t dt)
{
   Int_t trace_dled_length = fTrace_length - dt;
   std::vector<Float_t> trace_dled;
   trace_dled.reserve(trace_dled_length);

   for (int i = 0; i < trace_dled_length; ++i)
   {
      trace_dled.push_back(fAmplitude_Array.at(i+dt) - fAmplitude_Array.at(i));
   }

   return trace_dled;
}

//______________________________________________________________________________
std::vector<Float_t> Trace::TimeDLED(Int_t dt)
{
   Int_t trace_dled_length = fTrace_length - dt;
   std::vector<Float_t> time_dled;
   time_dled.reserve(trace_dled_length);

   for (int i = 0; i < trace_dled_length; ++i)
   {
      time_dled.push_back(fTime_Array.at(i+dt));
   }

   return time_dled;
}

Float_t Trace::GetMean(std::vector<Float_t> vec)
{
   // Compute mean of a vector.
   // From: https://stackoverflow.com/questions/7616511/calculate-mean-and-standard-deviation-from-a-vector-of-samples-in-c-using-boos/12405793#12405793

   Float_t sum = std::accumulate(vec.begin(), vec.end(), 0.0);

   return (sum / vec.size());
}

Float_t Trace::GetStdDev(std::vector<Float_t> vec)
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

std::vector<Int_t> Trace::FindPeaks(std::vector<Float_t> vec, Bool_t dled_bool, Int_t dt)
{   
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
      return emptyVec;
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

   return signals;
}

//______________________________________________________________________________
//void Trace::FindPeaks()
//{
//   // Find peaks of the trace
//
//   TCanvas *c1 = new TCanvas("c1","c1",10,10,1000,900);
//
//   TH1F *traceh = new TH1F("Traceh", "Trace;Time [ns];Amplitude [mV]", fTrace_length - 1, fTime_Array.data());
//   for (int i = 0; i < fTrace_length - 1; ++i)
//   {
//      traceh->Fill(fTime_Array,fAmplitude_Array);
//   }
//
//   traceh->Draw();
//
//   TSpectrum *s = new TSpectrum(100);
//   Int_t nfound = s->Search(traceh,2);
//   TH1 *tracehb = s->Background(traceh);
//
//   if (tracehb) c1->Update();
//
//
//
//}

//______________________________________________________________________________
void Trace::Paint(Option_t *option)
{
   TVirtualGraphPainter *painter = TVirtualGraphPainter::GetPainter();
   TGraph *trace_graph = new TGraph(fTrace_length, fTime_Array.data(), fAmplitude_Array.data());
   if(painter) painter->PaintHelper(trace_graph, option);
   delete trace_graph;
}

//______________________________________________________________________________
TraceHeader::TraceHeader() :  fDate(""), fSiBrand(""), fSiType(""), fIsDark(kTRUE), fAmplifier(""), fOut(2), 
                              fHVGen(""), fHV(0.0), fVoltagePulser(0.0)
{
   // Default constructor with only default values.
}

//______________________________________________________________________________
TraceHeader::TraceHeader(std::string &date, std::string &brand, std::string &type, Bool_t isdark, std::string &amplifier,
                        UInt_t chn, std::string &hvgen, Float_t hv, Float_t vpulser) :  fDate(date), 
                        fSiBrand(brand), fSiType(type), fIsDark(isdark), fAmplifier(amplifier), fOut(chn), fHVGen(hvgen), fHV(hv), fVoltagePulser(vpulser)
{
   // Specific constructor filling each class member.
}

//______________________________________________________________________________
TraceHeader::~TraceHeader()
{

}

//______________________________________________________________________________
void TraceHeader::Build(std::string &filename)
{
   // Fill TraceHeader members from file with name filename.
   // 
   // Filename is of the form:
   // 
   // DATE_SITYPE_MODE[_VPULSER]_HVGEN_HV_AMPLIFIER[_CHN]_NUMEVENTS_SUBRUN.dat
   // 
   // where:
   // 
   // DATE: date in the format YYYYMMDD
   // SITYPE: type of silicon, from which brand can be inferred
   // MODE: it can be DARK or LED
   // VPULSER: voltage of the pulser; present only if MODE == LED, in the format FFFF (voltage * 1000)
   // HVGEN: HV generator brand (e.g. AGILENT)
   // HV: value of HV
   // AMPLIFIER: brand of amplifier (e.g. AS for AdvanSid, PD for Padua)
   // CHN: channel number; present only if AMPLIFIER == AS
   // NUMEVENTS: number of events contained in the file
   // SUBRUN: if more runs are taken with the same settings (HV, HVGEN, AMPLIFIER etc...), denotes the run number

   SetDate(filename.substr(0,8));
   
   if(filename.find("HD3") != filename.length()){
      SetBrand("FBK");
      std::size_t pos = filename.find("HD3");
      SetType(filename.substr(pos, pos+7));
   }

   std::size_t pos_is_dark = filename.length();

   if(filename.find("DARK") != filename.length()){
      SetDark(kTRUE);
      pos_is_dark = filename.find("DARK");
      SetVoltagePulser(0.0);
   }
   
   if(filename.find("LED") != filename.length()){
      SetDark(kFALSE);
      pos_is_dark = filename.find("LED");
      Float_t vpulser = std::stof(filename.substr(pos_is_dark + 4, 4))/1000.0;
      SetVoltagePulser(vpulser);
   }

   std::size_t pos_ampl_beg = filename.find("_", pos_is_dark + 4) + 1;
   std::size_t pos_ampl_end = filename.find("_", pos_ampl_beg);

   SetHVGen(filename.substr(pos_ampl_beg, pos_ampl_end - pos_ampl_beg));

   SetHV(std::stof(filename.substr(pos_ampl_end + 1, 2)));

   std::string amplifier = filename.substr(pos_ampl_end + 4, 2);
   SetAmplifier(amplifier);

   if(amplifier == "AS"){
      UInt_t chn = std::stof(filename.substr(pos_ampl_end + 7, 1));
      SetChannelOut(chn);
   }
   else
   {
      SetChannelOut(0);
   }

   return;
}

//______________________________________________________________________________
RunSipm::RunSipm(const char *name, const char *title) : TTask(name,title)
{

}

//______________________________________________________________________________
RunSipm::~RunSipm()
{
   std::cout << "RunSipm task terminating." << std::endl;
}

//______________________________________________________________________________
void RunSipm::Exec(Option_t * /*option*/)
{
   std::cout << "RunSipm task executing." << std::endl;
}

//______________________________________________________________________________
Decode::Decode(const char *name, const char *title, std::string &filename, Int_t num_events) : TTask(name,title), fOutfile(0), fTraceTree(0), fHeaderTree(0),
                                                                                               fTrace(0), fTraceHeader(0), fFilename(filename), fNumevents(num_events)
{

}

//______________________________________________________________________________
Decode::~Decode()
{
   std::cout << "Deleting objects and closing." << std::endl;
   if(fOutfile) { fOutfile->Close(); }
   if(fTrace) { delete fTrace; fTrace = 0; }
   if(fTraceHeader) { delete fTraceHeader; fTraceHeader = 0; }
}

//______________________________________________________________________________
void Decode::DecodeInit()
{

   std::cout << "Initializing Decode Task ... " << std::endl;
   std::cout << "\t * initializing output TFile, TTrees and Branches ... " << std::endl;

   fTrace         = new Trace();
   fTraceHeader   = new TraceHeader();
   std::string outfilename = fFilename;
   outfilename.replace(outfilename.end()-4, outfilename.end(), ".root", 5);
   fOutfile       = new TFile(outfilename.c_str(), "RECREATE");
   fHeaderTree     = new TTree("Trace_header", "Header containing info about the acquisition.");
   fTraceTree    = new TTree("Trace", "Tree containing traces.");

   fHeaderTree->Branch("HeaderBranch", &fTraceHeader);
   fTraceTree->Branch("TraceBranch", &fTrace);

   std::cout << "... initializing Decode Task done!" << std::endl;
}

//______________________________________________________________________________
Int_t Decode::DecodeProcess()
{

   std::cout << "Processing Decode Task ... " << std::endl;

   std::cout << "\t* filling trees ... " << std::endl;

   // fill TraceHeader object from filename and put into the tree

   fTraceHeader->Build(fFilename);

   fHeaderTree->Fill();

   // Decode a number of events equal to num_events
   // from file filename.
   // 
   // Adapted from read_binary.cpp inside DRS4 software package.

   struct FHEADER {
      char           tag[3];
      char           version;
   };   // file header

   struct THEADER {
      char           time_header[4];
   };   // time header
   
   struct BHEADER {
      char           bn[2];
      unsigned short board_serial_number;
   };   // board serial number header
   
   struct EHEADER {
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
   };   // event header
   
   struct TCHEADER {
      char           tc[2];
      unsigned short trigger_cell;
   };   // trigger cell header
   
   struct CHEADER {
      char           c[1];
      char           cn[3];
   };   // channel header


   FHEADER  file_header;
   THEADER  time_header;
   BHEADER  board_header;
   EHEADER  event_header;
   TCHEADER trigger_cell_header;
   CHEADER  channel_header;

   const Int_t m_trace_length = fTrace->Get_trace_length();
   
   unsigned int scaler;
   std::vector<unsigned short> voltage;
   voltage.resize(m_trace_length);
   float time[16][4][m_trace_length], waveform[16][4][m_trace_length];
   float bin_width[16][4][m_trace_length];
   int i, j, b, chn, n, chn_index, n_boards;
   int channel_no;
   int channel_max = 4;
   int channels_acquired = 0;
   float t1 = 0, t2 = 0, dt = 0;
   Int_t ret_value = 0;

   // binary file definition and opening
   std::ifstream binary_file;

   binary_file.open(fFilename, std::ifstream::binary);

   // file_header reading
   binary_file.read((char *)&file_header, sizeof(FHEADER));

   if(!binary_file){
      std::cerr << "\t* error in reading file header from " << fFilename << ". Exiting ... " << std::endl;
      ret_value = 0;
   }

   if (file_header.tag[0] != 'D' || file_header.tag[1] != 'R' || file_header.tag[2] != 'S') {
      std::cerr << "\t* found invalid file header in file " << fFilename << ". Exiting ... " << std::endl;
      ret_value = 0;
   }
   
   if (file_header.version != '2') {
      std::cerr << "\t* found invalid file version in file " << fFilename << ". Exiting ... " << std::endl;
      ret_value = 0;
   }

   // time_header reading 
   binary_file.read((char *)&time_header, sizeof(time_header));

   if(!binary_file){
      std::cerr << "\t* error in reading time header from " << fFilename << ". Exiting ... " << std::endl;
      ret_value = 0;
   }

   if (memcmp(time_header.time_header, "TIME", 4) != 0) {
      std::cerr << "\t* found invalid time header in file " << fFilename << ". Exiting ... " <<  std::endl;
      ret_value = 0;
   }

   // board header reading

   for (b = 0 ; ; b++) {
      // read board header
      binary_file.read((char *)&board_header, sizeof(board_header));

      if(!binary_file){
         std::cerr << "\t* error in reading board header from " << fFilename << ". Exiting ... " << std::endl;
         ret_value = 0;
      }

      if (memcmp(board_header.bn, "B#", 2) != 0) {
         // probably event header found
         binary_file.seekg(-sizeof(board_header), std::ios::cur);
         break;
      }
      
      // read time bin widths
      memset(bin_width[b], 0, sizeof(float)*16*4*m_trace_length);

      for (channel_no = 0 ; channel_no < channel_max; channel_no++){

         // read channel header
         binary_file.read((char *)&channel_header, sizeof(channel_header));

         if(!binary_file){
            std::cerr << "\t* error in reading channel header from " << fFilename << ". Exiting ... " << std::endl;
            ret_value = 0;
         }

         if (channel_header.c[0] != 'C') {
            // event header found
            binary_file.seekg(-sizeof(channel_header), std::ios::cur);
            break;
         }

         channels_acquired++;

         i = channel_header.cn[2] - '0' - 1;
         //printf("Found timing calibration for channel #%d\n", i+1);

         binary_file.read((char *)&bin_width[b][i][0], m_trace_length*sizeof(float));

         if(!binary_file){
            std::cerr << "\t* error in reading bin width from " << fFilename << ". Exiting ... " << std::endl;
            ret_value = 0;
         }

         // fix for 2048 bin mode: float channel
         if (bin_width[b][i][1023] > 10 || bin_width[b][i][1023] < 0.01) {
            for (j=0 ; j<512 ; j++)
               bin_width[b][i][j+512] = bin_width[b][i][j];
         }
      }
   }

   std::cout << "\t* number of DRS4 channels acquired: " << channels_acquired << std::endl;

   n_boards = b;
   
   // loop over all events in the data file
   for(n=0 ; n < fNumevents; n++)
   {
      // read event header
      binary_file.read((char *)&event_header, sizeof(event_header));

      if(!binary_file){
         std::cerr << "\t* error in reading event header from " << fFilename << ". Exiting ... " << std::endl;
         std::cerr << "\t* probably you are trying to read more events than the ones contained in the file ... " << std::endl;
         ret_value = 0;
         n--;
         break;
      }
      
      //std::cout << "Read ev\t" << n <<std::endl;
          
      //printf("Found event #%d %d %d\n", eh.event_serial_number, eh.second, eh.millisecond);
      
      // loop over all boards in data file
      for (b=0 ; b<n_boards ; b++) {

         binary_file.read((char *)&board_header, sizeof(board_header));

         if(!binary_file){
            std::cerr << "\t* error in reading event header from " << fFilename << ". Exiting ... " << std::endl;
            ret_value = 0;
         }
         
         if (memcmp(board_header.bn, "B#", 2) != 0) {
            std::cerr <<"Invalid board header in file"<<std::endl;
            ret_value = 0;
         }

         // read trigger cell header
         binary_file.read((char *)&trigger_cell_header, sizeof(trigger_cell_header));

         if(!binary_file){
            std::cerr << "\t* error in reading trigger cell header from " << fFilename << ". Exiting ... " << std::endl;
            ret_value = 0;
         }
         
         if (memcmp(trigger_cell_header.tc, "T#", 2) != 0) {
            std::cerr<<"\t* found invalid trigger cell header in file"<<std::endl;
            ret_value = 0;
         }

         // reach channel data
         for (chn=0 ; chn<channels_acquired ; chn++) {
            // read channel header
            binary_file.read((char *)&channel_header, sizeof(channel_header));

            if(!binary_file){
               std::cerr << "\t* error in reading channel header from " << fFilename << ". Exiting ... " << std::endl;
               ret_value = 0;
            }

            if (channel_header.c[0] != 'C') {
               // event header found
               binary_file.seekg(-sizeof(channel_header), std::ios::cur);
               break;
            }

            chn_index = channel_header.cn[2] - '0' - 1;

            binary_file.read((char *)&scaler, sizeof(int));

            if(!binary_file){
               std::cerr << "\t* error in reading scaler from " << fFilename << ". Exiting ... " << std::endl;
               ret_value = 0;
            }

            // read waveform from binary file
            binary_file.read((char *)voltage.data(), m_trace_length*sizeof(unsigned short ));

            if(!binary_file){
               std::cerr << "\t* error in reading trace from " << fFilename << ". Exiting ... " << std::endl;
               ret_value = 0;
            }
            
            for (i=0 ; i<m_trace_length ; i++) {
               // convert data to milli Volts
               waveform[b][chn_index][i] = (voltage.at(i) / 65536. + event_header.range/1000.0 - 0.5)*1000.0;
               
               // calculate time for this cell
               for (j=0,time[b][chn_index][i]=0 ; j<i ; j++)
                  time[b][chn_index][i] += bin_width[b][chn_index][(j+trigger_cell_header.trigger_cell) % m_trace_length];
            }

         } // end for loop on DRS channels
         
         // if more than 1 DRS4 channels was acquired, align cell #0 of all the other channels
         if(channels_acquired > 1)
         {
            t1 = time[b][0][(m_trace_length-trigger_cell_header.trigger_cell) % m_trace_length];
            for (chn=1 ; chn<channels_acquired ; chn++) {
               t2 = time[b][chn][(m_trace_length-trigger_cell_header.trigger_cell) % m_trace_length];
               dt = t1 - t2;
               for (i=0 ; i<m_trace_length ; i++)
                  time[b][chn][i] += dt;
            } 
         }

         // fill the tree
         for (chn = 0; chn < channels_acquired; ++chn)
         {
            fTrace->Build(n, waveform[b][chn], time[b][chn]);
            fTraceTree->Fill();
         }

         if(n % 5000 == 0){
            std::cout << "\t* read event # " << n << std::endl;
         }

         ret_value = n;

      } // end for loop on boards (if more than one)
   } // end for loop on number of events

   if(binary_file.is_open()){
      std::cout << "\t* closing binary file ... " << std::endl;
      binary_file.close();
      std::cout << "\t* binary file closed ... " << std::endl;
   }

   std::cout << "... processing Decode Task done!" << std::endl;

   return ret_value;
}

//______________________________________________________________________________
void Decode::DecodeTerminate(Int_t ret_value)
{

   std::cout << "Terminating Decode Task ... " << std::endl;

   Int_t num_events_read = fTraceTree->GetEntries();

   std::cout << "\t* " << num_events_read << " events read from " << fFilename << " ... " << std::endl;

   if (!ret_value && num_events_read == 0)
   {
      std::cerr << "\t* problem in reading binary file ... " << std::endl;
   }
   else
   {
      std::cout << "\t* writing trees to file and closing ... " << std::endl;
      fOutfile->Write(0,TObject::kOverwrite);
      fHeaderTree->Print();
      fTraceTree->Print();
   }

   fOutfile->Close();
   delete fTrace;
   fTrace = 0;
   delete fTraceHeader;
   fTraceHeader = 0;

   std::cout << "\t* file closed and objects deleted ... " << std::endl;
   std::cout << "... terminating Decode Task done!" << std::endl;

   return;
}

//______________________________________________________________________________
void Decode::Exec(Option_t * /*option*/)
{
   std::cout << "Executing Decode Task ... " << std::endl;
   DecodeInit();
   Int_t retvalue = DecodeProcess();
   DecodeTerminate(retvalue);
   std::cout << "... executing Decode Task done!" << std::endl;
}

//______________________________________________________________________________
CharacterizeSipm::CharacterizeSipm(const char *name, const char *title, std::string &filename) : TTask(name,title), fInfile(0), fTraceTree(0), fHeaderTree(0),
                                                                                                 fTrace(0), fTraceHeader(0), fFilename(filename)
{

}

//______________________________________________________________________________
CharacterizeSipm::~CharacterizeSipm()
{
   std::cout << "Deleting objects and closing." << std::endl;
   if(fInfile) { fInfile->Close(); }
   if(fTrace) { delete fTrace; fTrace = 0; }
   if(fTraceHeader) { delete fTraceHeader; fTraceHeader = 0; }
}

//______________________________________________________________________________
void CharacterizeSipm::CharacterizeSipmInit()
{

   std::cout << "Initializing CharacterizeSipm Task ... " << std::endl;
   std::cout << "\t * initializing output TFile, TTrees and Branches ... " << std::endl;

   // fTrace and fTraceHeader are not initialized with new operator since they will be filled
   // when looping through the TTree entries. So their initial value should be 0, as in the constructor

   fInfile        = TFile::Open(fFilename.c_str());
   fHeaderTree    = (TTree *)fInfile->Get("Trace_header");
   fTraceTree     = (TTree *)fInfile->Get("Trace");

   fHeaderTree->SetBranchAddress("HeaderBranch", &fTraceHeader);
   fTraceTree->SetBranchAddress("TraceBranch", &fTrace);

   std::cout << "... initializing CharacterizeSipm Task done!" << std::endl;
}

//______________________________________________________________________________
void CharacterizeSipm::CharacterizeSipmProcess(Option_t * option)
{
   TString opt = option;

   if (opt.Contains("NB")){
      gROOT->SetBatch(kFALSE);
   }
   else
   {
      gROOT->SetBatch(); 
   }

   gROOT->IsBatch();

   std::cout << "Processing CharacterizeSipm Task ... " << std::endl;

   fHeaderTree->GetEntry(); // fTraceHeader is now filled

   Bool_t dled_bool     = kFALSE;
   Bool_t reverse_bool  = kFALSE;

   std::string amplifier = fTraceHeader->GetAmplifier();

   if (amplifier.compare("AS") == 0)
   {
      dled_bool = kTRUE;
      reverse_bool = kTRUE;
   }

   if (dled_bool)
   {
      std::cout << "\t * DLED will be used on traces ... " << std::endl;
   }

   if (reverse_bool)
   {
      std::cout << "\t * traces will be reversed ... " << std::endl;
   }

   Int_t nentries = fTraceTree->GetEntries();

   std::cout << "\t * going to read " << nentries << " entries ... " << std::endl;

   std::vector<Float_t> amplitude;
   std::vector<Float_t> timet;
   std::vector<Float_t> dled_trace;
   std::vector<Float_t> dled_time;

   std::vector<Int_t> peak_pos;
   std::vector<Float_t> peak_x;
   std::vector<Float_t> peak_y;

   //Int_t trace_length = 1024;

   // find minimum and maximu of times and amplitudes

   //fTraceTree->SetEstimate(nentries*1024 + 1);
//
   //fTraceTree->Draw("fAmplitude_Array", "", "goff");
   //std::vector<Double_t> all_amplitudes;
   //all_amplitudes.assign(fTraceTree->GetV1(), fTraceTree->GetV1() + fTraceTree->GetSelectedRows());
   //Float_t min_amplitude = *std::min_element(all_amplitudes.begin(), all_amplitudes.end());
   //Float_t max_amplitude = *std::max_element(all_amplitudes.begin(), all_amplitudes.end());
//
   //fTraceTree->Draw("fTime_Array", "", "goff");
   //std::vector<Double_t> all_times;
   //all_times.assign(fTraceTree->GetV1(), fTraceTree->GetV1() + fTraceTree->GetSelectedRows());
   //Float_t min_time = *std::min_element(all_times.begin(), all_times.end());
   //Float_t max_time = *std::max_element(all_times.begin(), all_times.end());
//
   //std::cout << min_amplitude << " " << max_amplitude << std::endl;
   //std::cout << min_time << " " << max_time << std::endl;

   TCanvas *c1 = new TCanvas("c1", "c1", 5);
   c1->Divide(1,2);
   TGraph  trace_graph;
   TGraph  dled_graph;
   TGraph  peaks;

   //TGaxis xaxis_low(0.1, 0.1, 0.9, 0.1, min_time, max_time);
   //TGaxis xaxis_up(0.1, 0.9, 0.9, 0.9, min_time, max_time, 510, "U-");
   //TGaxis yaxis_left(0.1, 0.1, 0.1, 0.9, min_amplitude, max_amplitude);
   //TGaxis yaxis_right(0.9, 0.1, 0.9, 0.9, min_amplitude, max_amplitude, 510, "U+");
//
   //c1->cd(1);
   //xaxis_low.Draw();
   //xaxis_up.Draw();
   //yaxis_left.Draw();
   //yaxis_right.Draw();
   //c1->cd(2);
   //xaxis_low.Draw();
   //xaxis_up.Draw();
   //yaxis_left.Draw();
   //yaxis_right.Draw();

   // loop over Trace TTree

   for (Int_t j = 0; j < nentries; ++j)
   {
      if (j%1000 == 0)
      {
         std::cout << "\t * processing " << j << "th event ... " << std::endl;
      }
      
      fTraceTree->GetEntry(j);

      // reverse trace if needed e.g. if amplifier is AdvanSid

      if (reverse_bool)
      {
         fTrace->Negate();
      }

      amplitude = fTrace->GetAmplitudeArray();
      timet     = fTrace->GetTimeArray();

      // apply DLED to trace if needed e.g. if amplifier is AdvanSid
      // also, find peaks
      
      if (dled_bool)
      {
         dled_trace = fTrace->AmplitudeDLED(9);
         dled_time  = fTrace->TimeDLED(9);
         peak_pos   = fTrace->FindPeaks(dled_trace,kTRUE,9);
      }
      else
      {
         peak_pos = fTrace->FindPeaks(amplitude);
      }

      // get only positive peaks and create two vector containing the
      // times and amplitudes of the peaks

      for(UInt_t i = 0; i < peak_pos.size(); i++){
         if(peak_pos[i] != 1){
            continue;
         }
         else
         {         
            dled_bool ? peak_x.push_back(dled_time[i])  : peak_x.push_back(timet[i]);
            dled_bool ? peak_y.push_back(dled_trace[i]) : peak_y.push_back(amplitude[i]);
         }
      }

      trace_graph = TGraph(amplitude.size(), timet.data(), amplitude.data());
      dled_graph  = TGraph(dled_trace.size(), dled_time.data(), dled_trace.data());
      peaks       = TGraph(peak_x.size(), peak_x.data(), peak_y.data());

      peaks.SetMarkerStyle(22);
      peaks.SetMarkerColor(kRed);

      c1->Clear("D");
      c1->cd(1);
      trace_graph.Draw("ALP");
      c1->cd(2);
      dled_graph.Draw("ALP");
      peaks.Draw("Psame");
      c1->Update();

      while(!gSystem->ProcessEvents()) {
         if (c1->WaitPrimitive()==0)
         {
            break;
         }
      }

      peak_x.clear();
      peak_y.clear();
   } 

   delete c1;

   return;

}

//______________________________________________________________________________
void CharacterizeSipm::CharacterizeSipmTerminate()
{
   std::cout << "Terminating CharacterizeSipm Task ... " << std::endl;

   std::cout << "\t* closing file ... " << std::endl;

   fInfile->Close();
   delete fTrace;
   fTrace = 0;
   delete fTraceHeader;
   fTraceHeader = 0;

   std::cout << "\t* file closed and objects deleted ... " << std::endl;
   std::cout << "... terminating CharacterizeSipm Task done!" << std::endl;

   return;
}

//______________________________________________________________________________
void CharacterizeSipm::Exec(Option_t * option)
{
   std::cout << "Executing CharacterizeSipm Task ... " << std::endl;
   CharacterizeSipmInit();
   CharacterizeSipmProcess(option);
   CharacterizeSipmTerminate();
   std::cout << "... executing CharacterizeSipm Task done!" << std::endl;
}