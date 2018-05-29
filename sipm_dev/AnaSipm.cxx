#include "AnaSipm.h"
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <algorithm>
#include <functional>
#include "TFile.h"
#include "TTree.h"

ClassImp(Trace);
ClassImp(TraceHeader);
ClassImp(AnaSipm);

//______________________________________________________________________________
Trace::Trace() : fId(0), fTrace_length(1024), fAmplitude_Array(fTrace_length, 0), fTime_Array(fTrace_length, 0)
{

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

   fAmplitude_Array.assign(amplitude_array, amplitude_array + fTrace_length);
   fAmplitude_Array.assign(time_array, time_array + fTrace_length);

}

//______________________________________________________________________________
void Trace::Detail()
{
   std::cout << "Trace ID number :\t" << fId << std::endl;
   //std::cout << "Trace First voltage :\t" << fAmplitude_Array[0] << std::endl;
   //std::cout << "Trace First time :\t" << fTime_Array[0] << std::endl;
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
AnaSipm::AnaSipm()
{
   // Deafult constructor of AnaSipm class.
}

//______________________________________________________________________________
AnaSipm::~AnaSipm()
{

} 

//______________________________________________________________________________
Int_t AnaSipm::Decode(Trace *trace, TraceHeader *trace_header, std::string filename, Int_t num_events)
{

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

   const Int_t m_trace_length = trace->Get_trace_length();
   
   unsigned int scaler;
   std::vector<unsigned short> voltage;
   voltage.resize(m_trace_length);
   float time[16][4][m_trace_length], waveform[16][4][m_trace_length];
   float bin_width[16][4][m_trace_length];
   int i, j, b, chn, n, chn_index, n_boards;
   int channel_no;
   int channel_no_max = 1;
   //float t1, t2, dt;
   Int_t ret_value = 0;
   Int_t num_events_read = 0;

   std::string outfilename = filename;
   outfilename.replace(outfilename.end()-4, outfilename.end(), ".root", 5);

   TFile *trace_out         = new TFile(outfilename.c_str(), "recreate");
   TTree *trace_header_tree = new TTree("Trace_header", "Header containing info about the acquisition.");
   TTree *trace_tree        = new TTree("Trace", "Tree containing traces.");

   trace_header_tree->Branch("HeaderBranch", &trace_header);
   trace_tree->Branch("TraceBranch", "Trace", &trace);

   // fill TraceHeader object from filename and put into the tree

   trace_header->Build(filename);

   trace_header_tree->Fill();

   // binary file definition and opening
   std::ifstream binary_file;

   binary_file.open(filename, std::ifstream::binary);

   // file_header reading
   binary_file.read((char *)&file_header, sizeof(FHEADER));

   if(!binary_file){
      std::cerr << "Error in reading file header from " << filename << ". Exiting." << std::endl;
      ret_value = 0;
   }

   if (file_header.tag[0] != 'D' || file_header.tag[1] != 'R' || file_header.tag[2] != 'S') {
      std::cerr << "Found invalid file header in file " << filename << ". Exiting." << std::endl;
      ret_value = 0;
   }
   
   if (file_header.version != '2') {
      std::cerr << "Found invalid file version in file " << filename << ". Exiting." << std::endl;
      ret_value = 0;
   }

   // time_header reading 
   binary_file.read((char *)&time_header, sizeof(time_header));

   if(!binary_file){
      std::cerr << "Error in reading time header from " << filename << ". Exiting." << std::endl;
      ret_value = 0;
   }

   if (memcmp(time_header.time_header, "TIME", 4) != 0) {
      std::cerr << "Invalid time header in file " << filename << ". Exiting." <<  std::endl;
      ret_value = 0;
   }

   // board header reading

   for (b = 0 ; ; b++) {
      // read board header
      binary_file.read((char *)&board_header, sizeof(board_header));

      if(!binary_file){
         std::cerr << "Error in reading board header from " << filename << ". Exiting." << std::endl;
         ret_value = 0;
      }

      if (memcmp(board_header.bn, "B#", 2) != 0) {
         // probably event header found
         binary_file.seekg(-sizeof(board_header), std::ios::cur);
         break;
      }
      
      // read time bin widths
      memset(bin_width[b], sizeof(bin_width[0]), 0);

      for (channel_no = 0 ; channel_no < channel_no_max; channel_no++){

         // read channel header
         binary_file.read((char *)&channel_header, sizeof(channel_header));

         if(!binary_file){
            std::cerr << "Error in reading channel header from " << filename << ". Exiting." << std::endl;
            ret_value = 0;
         }

         if (channel_header.c[0] != 'C') {
            // event header found
            binary_file.seekg(-sizeof(channel_header), std::ios::cur);
            break;
         }

         i = channel_header.cn[2] - '0' - 1;
         //printf("Found timing calibration for channel #%d\n", i+1);

         binary_file.read((char *)&bin_width[b][i][0], m_trace_length*sizeof(float));

         if(!binary_file){
            std::cerr << "Error in reading bin width from " << filename << ". Exiting." << std::endl;
            ret_value = 0;
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
   for(n=0 ; n < num_events; n++)
   {
      // read event header
      binary_file.read((char *)&event_header, sizeof(event_header));

      if(!binary_file){
         std::cerr << "Error in reading event header from " << filename << ". Exiting." << std::endl;
         ret_value = 0;
      }
      
      //std::cout << "Read ev\t" << n <<std::endl;
          
      //printf("Found event #%d %d %d\n", eh.event_serial_number, eh.second, eh.millisecond);
      
      // loop over all boards in data file
      for (b=0 ; b<n_boards ; b++) {

         binary_file.read((char *)&board_header, sizeof(board_header));

         if(!binary_file){
            std::cerr << "Error in reading event header from " << filename << ". Exiting." << std::endl;
            ret_value = 0;
         }
         
         if (memcmp(board_header.bn, "B#", 2) != 0) {
            std::cerr <<"Invalid board header in file"<<std::endl;
            ret_value = 0;
         }

         // read trigger cell header
         binary_file.read((char *)&trigger_cell_header, sizeof(trigger_cell_header));

         if(!binary_file){
            std::cerr << "Error in reading trigger cell header from " << filename << ". Exiting." << std::endl;
            ret_value = 0;
         }
         
         if (memcmp(trigger_cell_header.tc, "T#", 2) != 0) {
            std::cerr<<"Invalid trigger cell header in file"<<std::endl;
            ret_value = 0;
         }

         // reach channel data
         for (chn=0 ; chn<1 ; chn++) {
            // read channel header
            binary_file.read((char *)&channel_header, sizeof(channel_header));

            if(!binary_file){
               std::cerr << "Error in reading channel header from " << filename << ". Exiting." << std::endl;
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
               std::cerr << "Error in reading scaler from " << filename << ". Exiting." << std::endl;
               ret_value = 0;
            }

            // read waveform from binary file
            binary_file.read((char *)voltage.data(), m_trace_length*sizeof(unsigned short ));

             if(!binary_file){
               std::cerr << "Error in reading trace from " << filename << ". Exiting." << std::endl;
               ret_value = 0;
            }
            
            for (i=0 ; i<m_trace_length ; i++) {
               // convert data to milli Volts
               waveform[b][chn_index][i] = (voltage.at(i) / 65536. + event_header.range/1000.0 - 0.5)*1000.0;
               
               // calculate time for this cell
               for (j=0,time[b][chn_index][i]=0 ; j<i ; j++)
                  time[b][chn_index][i] += bin_width[b][chn_index][(j+trigger_cell_header.trigger_cell) % m_trace_length];
            }

            trace->Build(n, waveform[b][chn_index], time[b][chn_index]);
            
            if(n % 5000 == 0){
               trace->Detail();
            }

            trace_tree->Fill();

            ret_value = n;
            num_events_read = n;

         }
         
         // align cell #0 of all channels
         //t1 = time[b][0][(m_trace_length-trigger_cell_header.trigger_cell) % m_trace_length];
         //for (chn=1 ; chn<4 ; chn++) {
         //   t2 = time[b][chn][(m_trace_length-trigger_cell_header.trigger_cell) % m_trace_length];
         //   dt = t1 - t2;
         //   for (i=0 ; i<m_trace_length ; i++)
         //      time[b][chn][i] += dt;
         //}
      }
   }

   if(num_events_read != 0){
      trace_out->Write(0,TObject::kOverwrite);
      trace_header_tree->Print();
      trace_tree->Print();
      trace_out->ls();
   }
      
   trace_out->Close();

   if(binary_file.is_open()){
      binary_file.close();
   }

   return ret_value;
}