#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstring>

#include <AnaSipm.h>
#include "TTask.h"

void usage(){

   std::cout << "###############################################################################################" << std::endl;
   std::cout << "#                                                                                             #" << std::endl;
   std::cout << "#  This is decode_file executable. It decodes a file (if needed, e.g. for DRS4 binary files)  #" << std::endl;
   std::cout << "#  and saves the traces in a ROOT file.                                                       #" << std::endl;
   std::cout << "#                                                                                             #" << std::endl;
   std::cout << "###############################################################################################" << std::endl;
   std::cout << "                                                                                               " << std::endl;
   std::cout << "                                                                                               " << std::endl;
   std::cout << "  Usage:                                                                                       " << std::endl;
   std::cout << "                                                                                               " << std::endl;
   std::cout << "    decode_file --help [-h] : shows this help                                                  " << std::endl;
   std::cout << "                                                                                               " << std::endl;
   std::cout << "    decode_file filename num_traces : read num_traces events from filename. Output ROOT file   " << std::endl;
   std::cout << "                                      will have the same name of filename with .root extension " << std::endl;
   std::cout << "                                                                                               " << std::endl;

}

int main(int argc, char const *argv[])
{

   std::string filename;
   Int_t num_events;

   // command line arguments check

   if ( argc == 2 && (strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0) )
   {
      usage();
      return 0;
   }
   else if ( argc == 3 )
   {
      filename = std::string(argv[1]);

      std::ifstream infile(filename.c_str());
      
      // check if file exists and is not corrupted
      if (!infile.good())
      {
         std::cerr << "The filename you specified does not exist or is corrupted. Exiting." << std::endl;
         return 1; 
      } 

      // check that 2nd argument is a positive integer
      try{
         num_events = std::stoi(std::string(argv[2]));        
      }
      catch(std::exception const & e)
      {
         std::cerr << "You did not specify an integer, but something else. Exiting." << std::endl;
         return 1;
      }

      if (num_events < 0)
      {
         std::cerr << "A positive integer number is expected. Exiting." << std::endl;
         return 1;
      }
          
   }
   else
   {
      std::cerr << "\n";
      std::cerr << "Too many/too few arguments! The program requires two command line arguments.\n" << std::endl;
      usage();
      return 1; 
   }

   // initialize main task and decode subtask, then run them

   TTask *anasipm = new RunSipm("anasipm", "Process data from SiPM.");
   TTask *decode  = new Decode("decode", "Decode data from input file.", filename, num_events);

   anasipm->Add(decode);
   anasipm->ExecuteTask();

   delete anasipm;
   return 0;
}