#ifndef SIPM_H
#define SIPM_H

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// Sipm                                                                 //
//                                                                      //
// Description of the Sipm class                                        //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TObject.h"
#include <string>
 
class Trace : public TObject
{
   private:
      UInt_t       fId;                // id of the event
      Int_t        fTrace_length;      // length of the trace
      Float_t     *fAmplitude_Array;   //[fTrace_length]
      Float_t     *fTime_Array;        //[fTrace_length]

   public:
      Trace();
      Trace(UInt_t id);
      ~Trace();

      void Build(UInt_t id, Float_t *amplitude_array, Float_t *time_array);
      void Detail();

      int Get_trace_length() { return fTrace_length; }

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

class Sipm : public TObject
{
   public:
      Sipm();
      ~Sipm(); 

      Int_t Decode(Trace *trace, TraceHeader *trace_header, std::string filename, Int_t num_events);

      ClassDef(Sipm,1);

};
 
#endif