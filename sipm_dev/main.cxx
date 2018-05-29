#include <iostream>
#include <cstdlib>

#include <AnaSipm.h>

int main(int argc, char const *argv[])
{
   if (argc > 3)
   {
      std::cerr << "Too many arguments. The program requires only two command line argument." << std::endl;
      return 1;
   }

   std::string filename(argv[1]);
   Int_t num_events = atoi(argv[2]);

   AnaSipm ana_sipm;
   Trace *trace = new Trace();
   TraceHeader *trace_header = new TraceHeader();

   Int_t ret = ana_sipm.Decode(trace, trace_header, filename, num_events);

   if (!ret)
   {
      std::cerr << "Problem in reading binary file." << std::endl;
   }

   delete trace;
   delete trace_header;
   return 0;
}