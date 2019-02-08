#include <iostream>
#include <cstdlib>

#include <AnaSipm.h>
#include "TTask.h"

int main(int argc, char const *argv[])
{
   if (argc > 3)
   {
      std::cerr << "Too many arguments. The program requires only two command line argument." << std::endl;
      return 1;
   }

   std::string filename(argv[1]);
   Int_t num_events = atoi(argv[2]);

   TTask *anasipm = new RunSipm("anasipm", "Process data from SiPM.");
   TTask *decode  = new Decode("decode", "Decode data from input file.", filename, num_events);

   anasipm->Add(decode);
   anasipm->ExecuteTask();

   delete anasipm;
   return 0;
}