#ifndef ANASIPM_H
#define ANASIPM_H

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// Sipm                                                                 //
//                                                                      //
// Description of the Sipm class                                        //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TObject.h"
#include "TFile.h"
#include "TTree.h"
#include "TTask.h"
#include <string>
 
class Trace : public TObject
{
   private:
      UInt_t                  fId;                // id of the event
      Int_t                   fTrace_length;      // length of the trace
      std::vector<Float_t>    fAmplitude_Array;   // vector containing amplitudes
      std::vector<Float_t>    fTime_Array;        // vector containing times

   public:
      Trace();
      ~Trace();

      void Build(UInt_t id, Float_t *amplitude_array, Float_t *time_array);
      void Negate();
      //void FindPeaks();
      void Paint(Option_t *option);

      std::vector<Float_t> GetAmplitudeArray();
      std::vector<Float_t> GetTimeArray();

      std::vector<Float_t> AmplitudeDLED(Int_t dt);
      std::vector<Float_t> TimeDLED(Int_t dt);

      Float_t GetMean(std::vector<Float_t> vec);
      Float_t GetStdDev(std::vector<Float_t> vec);

      std::vector<Int_t> FindPeaks(std::vector<Float_t> vec, Bool_t dled_bool = kFALSE, Int_t dt = 0);     

      Int_t Get_trace_length() { return fTrace_length; }

      ClassDef(Trace,1);
};

class TraceHeader : public TObject
{
   private:
      std::string fDate;            // date in the format YYYYMMDD
      std::string fSiBrand;         // brand of silicon (FBK, Hamamatsu, SensL)
      std::string fSiType;          // type of silicon (e.g. HD3-2_1)
      Bool_t      fIsDark;          // acquisition in DARK or LED
      std::string fAmplifier;       // amplifier type
      UInt_t      fOut;             // out channel (e.g. 2 for AdvanSid board)
      std::string fHVGen;           // brand of HV generator
      Float_t     fHV;              // value of HV
      Float_t     fVoltagePulser;   // voltage of pulser (only in LED mode)
      //std::string fDateFirstEvent;  // time of first event of the acquisition
      //std::string fDateLastEvent;   // time of last event of the acquisition

   public:
      TraceHeader();
      TraceHeader(std::string &date, std::string &brand, std::string &type, Bool_t isdark, std::string &amplifier, UInt_t chn,
                  std::string &hvgen, Float_t hv, Float_t vpulser);
      ~TraceHeader();

      std::string GetDate()           { return fDate; };
      std::string GetBrand()          { return fSiBrand; };
      std::string GetType()           { return fSiType; };
      Bool_t      IsDark()            { return fIsDark; };
      std::string GetAmplifier()      { return fAmplifier; };
      UInt_t      GetChannelOut()     { return fOut; };
      std::string GetHVGen()          { return fHVGen; };
      Float_t     GetHV()             { return fHV; };
      Float_t     GetVoltagePulser()  { return fVoltagePulser; };
      //std::string GetDateFirstEvent() { return fDateFirstEvent; };
      //std::string GetDateLastEvent()  { return fDateLastEvent; };

      void SetDate(std::string date)                       { fDate = date; };
      void SetBrand(std::string brand)                     { fSiBrand = brand; };
      void SetType(std::string type)                       { fSiType = type; };
      void SetDark(Bool_t dark)                            { fIsDark = dark; };
      void SetAmplifier(std::string amplifier)             { fAmplifier = amplifier; };
      void SetChannelOut(UInt_t chn)                       { fOut = chn; };
      void SetHVGen(std::string hvgen)                     { fHVGen = hvgen; };
      void SetHV(Float_t hv)                               { fHV = hv; };
      void SetVoltagePulser(Float_t vpulser)               { fVoltagePulser = vpulser; };
      //void SetDateFirstEvent(std::string date_first_event) { fDateFirstEvent = date_first_event; };
      //void SetDateLastEvent(std::string date_last_event)   { fDateLastEvent = date_last_event; };
      
      void Build(std::string &filename);

      ClassDef(TraceHeader,1);
};

class RunSipm : public TTask
{
   public:
      RunSipm() {;}
      RunSipm(const char *name, const char *title);
      virtual ~RunSipm();
      void Exec(Option_t *option="");
      ClassDef(RunSipm,1);   // Run Reconstruction task
};

class Decode : public TTask {

   private:            
      TFile                *fOutfile;        // output TFile        
      TTree                *fTraceTree;      // tree for traces
      TTree                *fHeaderTree;     // tree for trace header
      Trace                *fTrace;          // Trace object
      TraceHeader          *fTraceHeader;    // TraceHeader object
      std::string          fFilename;        // name of input file
      Int_t                fNumevents;       // number of events to read

   public:
      Decode() {;}
      Decode(const char *name, const char *title, std::string &filename, Int_t num_events);
      virtual ~Decode();
      void Exec(Option_t *option="");
      void DecodeInit();
      Int_t DecodeProcess();
      void DecodeTerminate(Int_t ret_value);
      ClassDef(Decode,1);   // Geometry initialisation task
};

class CharacterizeSipm : public TTask {

   private:
      TFile                *fInfile;         // TFile to open
      TTree                *fTraceTree;      // tree for traces
      TTree                *fHeaderTree;     // tree for trace header
      Trace                *fTrace;          // Trace object
      TraceHeader          *fTraceHeader;    // TraceHeader object
      std::string          fFilename;        // name of input file

   public:
      CharacterizeSipm() {;}
      CharacterizeSipm(const char *name, const char *title, std::string &filename);
      virtual ~CharacterizeSipm();
      void Exec(Option_t *option="");
      void CharacterizeSipmInit();
      void CharacterizeSipmProcess(Option_t *option="");
      void CharacterizeSipmTerminate();
      ClassDef(CharacterizeSipm,1);
};

//class HistoManager {
//
//private:
//   TH1F  *fNtrack;
//   TH1F  *fNseg;
//   TH1F  *fTemperature;
//   TH1F  *fPx;
//   TH1F  *fPy;
//   TH1F  *fPz;
//   TH1F  *fRandom;
//   TH1F  *fMass2;
//   TH1F  *fBx;
//   TH1F  *fBy;
//   TH1F  *fMeanCharge;
//   TH1F  *fXfirst;
//   TH1F  *fXlast;
//   TH1F  *fYfirst;
//   TH1F  *fYlast;
//   TH1F  *fZfirst;
//   TH1F  *fZlast;
//   TH1F  *fCharge;
//   TH1F  *fNpoint;
//   TH1F  *fValid;
//
//public:
//   HistogramManager(TDirectory *dir);
//   virtual ~HistogramManager();
//
//   void Hfill(Event *event);
//
//   ClassDef(HistogramManager,1)  //Manages all histograms
//};
 
#endif