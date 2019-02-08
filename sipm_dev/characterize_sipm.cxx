#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstring>

#include <AnaSipm.h>
#include "TTask.h"
#include "TApplication.h"

void usage(){

   std::cout << "###############################################################################################" << std::endl;
   std::cout << "#                                                                                             #" << std::endl;
   std::cout << "#  This is characterize_sipm executable. It is used to charcaterize SiPMs                     #" << std::endl;
   std::cout << "#                                                                                             #" << std::endl;
   std::cout << "###############################################################################################" << std::endl;
   std::cout << "                                                                                               " << std::endl;
   std::cout << "                                                                                               " << std::endl;
   std::cout << "  Usage:                                                                                       " << std::endl;
   std::cout << "                                                                                               " << std::endl;
   std::cout << "    characterize_sipm --help [-h] : shows this help                                            " << std::endl;
   std::cout << "                                                                                               " << std::endl;
   std::cout << "    characterize_sipm filename : read events from filename.                                    " << std::endl;
   std::cout << "                                                                                               " << std::endl;

}

int main(int argc, char* argv[])
{

   std::string filename;

   // command line arguments check

   std::string batch("NB");
   Bool_t kBatch = kFALSE;

   if (argc == 2){
      if ((strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0))
      {
         usage();
         return 0;
      }
      else
      {
         filename = std::string(argv[1]);
      }
   }
   else if (argc == 3)
   {
      if (strcmp(argv[1], "-b") == 0)
      {
         batch = "B";
         kBatch = kTRUE;
      }

      filename = std::string(argv[2]);
   }
   else
   {
      std::cerr << "\n";
      std::cerr << "Too many/too few arguments! The program requires two command line arguments.\n" << std::endl;
      usage();
      return 1; 
   }

   TApplication *app = new TApplication("app", &argc, argv);

   // initialize main task and decode subtask, then run them

   TTask *anasipm = new RunSipm("anasipm", "Process data from SiPM.");
   TTask *characterize = new CharacterizeSipm("characterize", "characterize sipm", filename);

   anasipm->Add(characterize);
   anasipm->ExecuteTask(batch.c_str());

   delete anasipm;

   if (!kBatch)
   {
      app->Run();  
   }

   return 0;
}