#include "Sipm.h"
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include "TFile.h"
#include "TTree.h"

//using namespace std;

ClassImp(Trace);
ClassImp(TraceHeader);
ClassImp(Sipm);

//______________________________________________________________________________
Trace::Trace() : fId(0), fTrace_length(1024), fAmplitude_Array(0), fTime_Array(0)
{

}

//______________________________________________________________________________
Trace::Trace(UInt_t id): fId(id), fTrace_length(1024)
{
   fAmplitude_Array = new Float_t[fTrace_length];
   fTime_Array      = new Float_t[fTrace_length];
}

//______________________________________________________________________________
Trace::~Trace()
{

   if(fAmplitude_Array){
      delete[] fAmplitude_Array;
   }

   fAmplitude_Array = 0;

   if(fTime_Array){
      delete[] fTime_Array;
   }

   fTime_Array = 0;

}

//______________________________________________________________________________
void Trace::Build(UInt_t id, Float_t *amplitude_array, Float_t *time_array)
{

   fId              = id;

   for (Int_t i = 0; i < fTrace_length; ++i)
   {
      fAmplitude_Array[i] = amplitude_array[i];
      fTime_Array[i]      = time_array[i];
   }

}

//______________________________________________________________________________
void Trace::Detail()
{
   std::cout << "Trace ID number :\t" << fId << std::endl;
   if(fAmplitude_Array)
      std::cout << "Trace First voltage :\t" << fAmplitude_Array[0] << std::endl;
   if(fTime_Array)
      std::cout << "Trace First time :\t" << fTime_Array[0] << std::endl;
}

//______________________________________________________________________________
TraceHeader::TraceHeader() :  fDate(""), fSiBrand(""), fSiType(""), fIsDark(kTRUE), fAmplifier(""), fOut(2), 
                              fHVGen(""), fHV(0.0), fVoltagePulser(0.0)
{

}

//______________________________________________________________________________
TraceHeader::TraceHeader(std::string &date, std::string &brand, std::string &type, Bool_t isdark, std::string &amplifier,
                        UInt_t chn, std::string &hvgen, Float_t hv, Float_t vpulser) :  fDate(date), 
                        fSiBrand(brand), fSiType(type), fIsDark(isdark), fAmplifier(amplifier), fOut(chn), fHVGen(hvgen), fHV(hv), fVoltagePulser(vpulser)
{

}

//______________________________________________________________________________
TraceHeader::~TraceHeader()
{

}

//______________________________________________________________________________
void TraceHeader::Build(std::string &filename)
{
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
Sipm::Sipm()
{

}

//______________________________________________________________________________
Sipm::~Sipm()
{

} 

//______________________________________________________________________________
Int_t Sipm::Decode(Trace *trace, TraceHeader *trace_header, std::string filename, Int_t num_events)
{

   // Decode a number of events equal to num_events
   // from file filename

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
   //unsigned short voltage[m_trace_length];
   float time[16][4][m_trace_length], waveform[16][4][m_trace_length];
   float bin_width[16][4][m_trace_length];
   int i, j, b, chn, n, chn_index, n_boards;
   int channel_no;
   int channel_no_max = 4;
   float t1, t2, dt;

   std::string outfilename = filename;
   outfilename.replace(outfilename.end()-4, outfilename.end(), ".root", 5);

   TFile *trace_out         = new TFile(outfilename.c_str(), "recreate");
   TTree *trace_header_tree = new TTree("Trace_header", "Header containing info about the acquisition.");
   TTree *trace_tree        = new TTree("Trace", "Tree containing traces.");

   trace_header_tree->Branch("HeaderBranch", &trace_header);
   trace_tree->Branch("TraceBranch", &trace);

   trace_header->Build(filename);

   trace_header_tree->Fill();

   // binary file reading
   std::ifstream binary_file;

   binary_file.open(filename, std::ifstream::binary);

   // file_header reading 
   binary_file.read((char *)&file_header, sizeof(FHEADER));

   if(!binary_file){
      std::cerr << "Error in reading file header from " << filename << ". Exiting." << std::endl;
      binary_file.close();
      return 0;
   }

   if (file_header.tag[0] != 'D' || file_header.tag[1] != 'R' || file_header.tag[2] != 'S') {
      std::cerr << "Found invalid file header in file " << filename << ". Exiting." << std::endl;
      binary_file.close();
      return 0;
   }
   
   if (file_header.version != '2') {
      std::cerr << "Found invalid file version in file " << filename << ". Exiting." << std::endl;
      binary_file.close();
      return 0;
   }

   // time_header reading 
   binary_file.read((char *)&time_header, sizeof(time_header));

   if(!binary_file){
      std::cerr << "Error in reading time header from " << filename << ". Exiting." << std::endl;
      binary_file.close();
      return 0;
   }

   if (memcmp(time_header.time_header, "TIME", 4) != 0) {
      std::cerr << "Invalid time header in file " << filename << ". Exiting." <<  std::endl;
      binary_file.close();
      return 0;
   }

   // board header reading

   for (b = 0 ; ; b++) {
      // read board header

      binary_file.read((char *)&board_header, sizeof(board_header));

      if(!binary_file){
         std::cerr << "Error in reading board header from " << filename << ". Exiting." << std::endl;
         binary_file.close();
         return 0;
      }

      if (memcmp(board_header.bn, "B#", 2) != 0) {
         // probably event header found
         binary_file.seekg(-sizeof(board_header), std::ios::cur);
         break;
      }
      
      // read time bin widths
      memset(bin_width[b], sizeof(bin_width[0]), 0);

      for (channel_no = 0 ; channel_no < channel_no_max; channel_no++){

         binary_file.read((char *)&channel_header, sizeof(channel_header));

         if(!binary_file){
            std::cerr << "Error in reading channel header from " << filename << ". Exiting." << std::endl;
            binary_file.close();
            return 0;
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
            binary_file.close();
            return 0;
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
         binary_file.close();
         return 0;
      }
      
      std::cout << "Read ev\t" << n <<std::endl;
          
      //printf("Found event #%d %d %d\n", eh.event_serial_number, eh.second, eh.millisecond);
      
      // loop over all boards in data file
      for (b=0 ; b<n_boards ; b++) {

         binary_file.read((char *)&board_header, sizeof(board_header));

         if(!binary_file){
            std::cerr << "Error in reading event header from " << filename << ". Exiting." << std::endl;
            binary_file.close();
            return 0;
         }
         
         if (memcmp(board_header.bn, "B#", 2) != 0) {
            //printf("Invalid board header in file \'%s\', aborting.\n", filename);
//             return 0;
            std::cerr <<"Invalid board header in file"<<std::endl;
         }

         binary_file.read((char *)&trigger_cell_header, sizeof(trigger_cell_header));

         if(!binary_file){
            std::cerr << "Error in reading trigger cell header from " << filename << ". Exiting." << std::endl;
            binary_file.close();
            return 0;
         }
         
         if (memcmp(trigger_cell_header.tc, "T#", 2) != 0) {
            //printf("Invalid trigger cell header in file \'%s\', aborting.\n", filename);
            std::cerr<<"Invalid trigger cell header in file"<<std::endl;
//             return 0;
         }

         // reach channel data
         for (chn=0 ; chn<4 ; chn++) {
            // read channel header

            binary_file.read((char *)&channel_header, sizeof(channel_header));

            if(!binary_file){
               std::cerr << "Error in reading channel header from " << filename << ". Exiting." << std::endl;
               binary_file.close();
               return 0;
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
               binary_file.close();
               return 0;
            }

            binary_file.read((char *)voltage.data(), m_trace_length*sizeof(unsigned short ));

             if(!binary_file){
               std::cerr << "Error in reading trace from " << filename << ". Exiting." << std::endl;
               binary_file.close();
               return 0;
            }
            
            for (i=0 ; i<m_trace_length ; i++) {
               // convert data to volts
               waveform[b][chn_index][i] = (voltage.at(i) / 65536. + event_header.range/1000.0 - 0.5);
               
               // calculate time for this cell
               for (j=0,time[b][chn_index][i]=0 ; j<i ; j++)
                  time[b][chn_index][i] += bin_width[b][chn_index][(j+trigger_cell_header.trigger_cell) % m_trace_length];
            }

            trace->Build(n, waveform[b][chn_index], time[b][chn_index]);
            trace->Detail();
            trace_tree->Fill();

         }
         
         // align cell #0 of all channels
         t1 = time[b][0][(m_trace_length-trigger_cell_header.trigger_cell) % m_trace_length];
         for (chn=1 ; chn<4 ; chn++) {
            t2 = time[b][chn][(m_trace_length-trigger_cell_header.trigger_cell) % m_trace_length];
            dt = t1 - t2;
            for (i=0 ; i<m_trace_length ; i++)
               time[b][chn][i] += dt;
         }
      }
   }

   trace_out->Write();
   trace_header_tree->Print();
   trace_tree->Print();
   trace_out->Close();

   binary_file.close();
   return num_events;
}