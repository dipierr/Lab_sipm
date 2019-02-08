/********************************************************************\

  Name:         drs_exam.cpp
  Created by:   Stefan Ritt

  Contents:     Simple example application to read out a DRS4
                evaluation board

  $Id: drs_exam.cpp 21308 2014-04-11 14:50:16Z ritt $

\********************************************************************/

#include <math.h>

#ifdef _MSC_VER

#include <windows.h>

#elif defined(OS_LINUX)

#define O_BINARY 0

#include <unistd.h>
#include <ctype.h>
#include <sys/ioctl.h>
#include <errno.h>

#define DIR_SEPARATOR '/'

#endif

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <vector>

#include "strlcpy.h"
#include "DRS.h"

/* ROOT includes */

#include "TFile.h"
#include "TTree.h"

/*------------------------------------------------------------------*/
/*
int find_remove_spikes(float *waveform_1, float *waveform_2){

   std::vector<int> spike_pos;

   for (int j = 5; j < 1024 - 5; ++j)
   {
      if(-waveform_2[j] + waveform_2[j+1] + waveform_2[j+2] - waveform_2[j+3] > 20.){
         spike_pos.push_back(j+1);
         if(waveform_2[j+2]>10.){
            spike_pos.push_back(j+2);
         }
      }
   }

   for (std::vector<int>::iterator it = spike_pos.begin(); it != spike_pos.end(); ++it)
   {
      if (*(it+1) - *(it) == 1)
      {
         waveform[*it] = (waveform[*it+2] - waveform[*it-1])*(*it - 1);
      }
   }


}

*/


int RemoveSpikes(bool cascading, int trigger_cell, float waveform[8][1024])
{
   int i, j, k, l, c, ofs, nChn;
   double x, y;
   int sp[8][10];
   int rsp[10], rot_sp[10];
   int n_sp[8], n_rsp;
   float wf[8][2048];
   int  nNeighbor, nSymmetric;

   bool rotated = true;
   int spos[1024];

   nChn = cascading ? 8 : 4;

   /* rotate waveform back relative to cell #0 */
   if (cascading) {
      if (rotated) {
         for (i=0 ; i<nChn ; i++)
            for (j=0 ; j<1024 ; j++)
               wf[i][(j+trigger_cell) % 1024] = waveform[i/2][(i%2)*1024+j];
      } else {
         for (i=0 ; i<nChn ; i++)
            for (j=0 ; j<1024 ; j++)
               wf[i][j] = waveform[i/2][(i%2)*1024+j];
      }
   } else {
      if (rotated) {
         for (i=0 ; i<nChn ; i++)
            for (j=0 ; j<1024 ; j++)
               wf[i][(j+trigger_cell) % 1024] = waveform[i][j];
      } else {
         for (i=0 ; i<nChn ; i++)
            for (j=0 ; j<1024 ; j++)
               wf[i][j] = waveform[i][j];
      }
   }


   memset(sp, 0, sizeof(sp));
   memset(n_sp, 0, sizeof(n_sp));
   memset(rsp, 0, sizeof(rsp));
   n_rsp = 0;

   /* find spikes with special high-pass filter */
   for (j=0 ; j<1024 ; j++) {
      for (i=0 ; i<nChn ; i++) {
         if (-wf[i][j]+wf[i][(j+1) % 1024]+wf[i][(j+2) % 1024]-wf[i][(j+3) % 1024] > 20) {
            if (n_sp[i] < 10) // record maximum of 10 spikes
               sp[i][n_sp[i]++] = j;
            else
               return 1;        // too many spikes -> something wrong
            spos[j]++;
         }
         if (-wf[i][j]+wf[i][(j+1) % 1024]+wf[i][(j+2) % 1024]-wf[i][(j+3) % 1024] < -20) {
            if (n_sp[i] < 10) // record maximum of 10 spikes
               sp[i][n_sp[i]++] = j;
            else
               return 1;        // too many spikes -> something wrong
            spos[j]++;
         }
      }
   }

   /* find spikes at cell #0 and #1023 */
   for (i=0 ; i<nChn ; i++) {
      if (wf[i][0]+wf[i][1]-2*wf[i][2] > 20) {
         if (n_sp[i] < 10)
            sp[i][n_sp[i]++] = 0;
      }
      if (-2*wf[i][1021]+wf[i][1022]+wf[i][1023] > 20) {
         if (n_sp[i] < 10)
            sp[i][n_sp[i]++] = 1020;
      }
   }

   /* go through all spikes and look for symmetric spikes and neighbors */
   for (i=0 ; i<nChn ; i++) {
      for (j=0 ; j<n_sp[i] ; j++) {
         /* check if this spike has a symmetric partner in any channel */
         for (k=nSymmetric=0 ; k<nChn ; k++) {
               for (l=0 ; l<n_sp[k] ; l++)
                  if (sp[i][j] == (1020-sp[k][l]+1024) % 1024) {
                     nSymmetric++;
                     break;
                  }
            }

         /* check if this spike has same spike in any other channels */
         for (k=nNeighbor=0 ; k<nChn ; k++)
            if (i != k) {
               for (l=0 ; l<n_sp[k] ; l++)
                  if (sp[i][j] == sp[k][l]) {
                     nNeighbor++;
                     break;
                  }
            }

         if (nSymmetric + nNeighbor >= 2) {
            /* if at least two matching spikes, treat this as a real spike */
            for (k=0 ; k<n_rsp ; k++)
               if (rsp[k] == sp[i][j])
                  break;
            if (n_rsp < 10 && k == n_rsp)
               rsp[n_rsp++] = sp[i][j];
         }
      }
   }

   /* rotate spikes according to trigger cell */
   if (rotated) {
      for (i=0 ; i<n_rsp ; i++)
         rot_sp[i] = (rsp[i] - trigger_cell + 1024) % 1024;
   } else {
      for (i=0 ; i<n_rsp ; i++)
         rot_sp[i] = rsp[i];
   }

   if(n_sp[1] != 0 ||  n_sp[2] != 0 || n_sp[3] != 0){
      return n_sp[1];
   }

   /* recognize spikes if at least one channel has it */
   for (k=0 ; k<n_rsp ; k++) {
      for (i=0 ; i<nChn ; i++) {

         if (cascading) {
            c = i/2;
            ofs = (i%2)*1024;
         } else {
            c = i;
            ofs = 0;
         }
         if (k < n_rsp-1 && rsp[k] == 0 && rsp[k+1] == 1020) {
            /* remove double spike */
            j = rot_sp[k] > rot_sp[k+1] ? rot_sp[k+1] : rot_sp[k];
            x = waveform[c][ofs+(j+1) % 1024];
            y = waveform[c][ofs+(j+6) % 1024];
            if (fabs(x-y) < 15) {
               waveform[c][ofs+(j+2) % 1024] = x + 1*(y-x)/5;
               waveform[c][ofs+(j+3) % 1024] = x + 2*(y-x)/5;
               waveform[c][ofs+(j+4) % 1024] = x + 3*(y-x)/5;
               waveform[c][ofs+(j+5) % 1024] = x + 4*(y-x)/5;
            } else {
               waveform[c][ofs+(j+2) % 1024] -= 50.0f;
               waveform[c][ofs+(j+3) % 1024] -= 50.0f;
               waveform[c][ofs+(j+4) % 1024] -= 50.0f;
               waveform[c][ofs+(j+5) % 1024] -= 50.0f;
            }
         } else {
            /* remove single spike */
            x = waveform[c][ofs+rot_sp[k]];
            y = waveform[c][ofs+(rot_sp[k]+3) % 1024];

            if (fabs(x-y) < 15) {
               waveform[c][ofs+(rot_sp[k]+1) % 1024] = x + 1*(y-x)/3;
               waveform[c][ofs+(rot_sp[k]+2) % 1024] = x + 2*(y-x)/3;
            } else {
               waveform[c][ofs+(rot_sp[k]+1) % 1024] -= 50.0f;
               waveform[c][ofs+(rot_sp[k]+2) % 1024] -= 50.0f;
            }
         }
      }
      if (k < n_rsp-1 && rsp[k] == 0 && rsp[k+1] == 1020)
         k++; // skip second half of double spike
   }

   /* uncomment to show unfixed spikes
   m_skipDisplay = true;
   for (i=0 ; i<1024 ; i++)
      for (j=0 ; j<nChn ; j++) {
         if (waveform[j][i] > 10 || waveform[j][i] < -10) {
            m_skipDisplay = false;
            break;
         }
   }
   */

   return 0;
}


int main()
{
   int i, j, nBoards;
   DRS *drs;
   DRSBoard *b;
   float time_array[8][1024];
   float wave_array[8][1024];
   //FILE  *f;

   // Initialize TFile to save ROOT file and TTree and branches

   TFile *file_drs = new TFile("20180420_drs4_all_channels_spike_removal_test_3.root", "recreate");
   TTree *tree_drs = new TTree("Waveforms", "Tree containing waveforms acquired from the DRS4");

   // When filling, each branch will take 1024 values for each waveform

   
   tree_drs->Branch("time_ch_1", time_array[0], "time_ch_1[1024]/F");
   tree_drs->Branch("volt_ch_1", wave_array[0], "volt_ch_1[1024]/F");
   tree_drs->Branch("time_ch_2", time_array[1], "time_ch_2[1024]/F");
   tree_drs->Branch("volt_ch_2", wave_array[1], "volt_ch_2[1024]/F");
   tree_drs->Branch("time_ch_3", time_array[2], "time_ch_3[1024]/F");
   tree_drs->Branch("volt_ch_3", wave_array[2], "volt_ch_3[1024]/F");
   tree_drs->Branch("time_ch_4", time_array[3], "time_ch_4[1024]/F");
   tree_drs->Branch("volt_ch_4", wave_array[3], "volt_ch_4[1024]/F");

   /* do initial scan */
   drs = new DRS();

   /* show any found board(s) */
   for (i=0 ; i<drs->GetNumberOfBoards() ; i++) {
      b = drs->GetBoard(i);
      printf("Found DRS4 evaluation board, serial #%d, firmware revision %d\n", 
         b->GetBoardSerialNumber(), b->GetFirmwareVersion());
   }

   /* exit if no board found */
   nBoards = drs->GetNumberOfBoards();
   if (nBoards == 0) {
      printf("No DRS4 evaluation board found\n");
      return 0;
   }

   /* continue working with first board only */
   b = drs->GetBoard(0);

   /* initialize board */
   b->Init();

   /* set sampling frequency */
   b->SetFrequency(5, true);

   /* enable transparent mode needed for analog trigger */
   b->SetTranspMode(1);

   /* set input range to -0.5V ... +0.5V */
   b->SetInputRange(0);

   /* use following line to set range to 0..1V */
   //b->SetInputRange(0.5);
   
   /* use following line to turn on the internal 100 MHz clock connected to all channels  */
   /* if set to 1, it takes in the internal clock waveforms */
   b->EnableTcal(0);

   /* use following lines to enable hardware trigger on CH1 at 50 mV positive edge */
   if (b->GetBoardType() >= 8) {        // Evaluaiton Board V4&5
      b->EnableTrigger(1, 0);           // enable hardware trigger
      b->SetTriggerSource(1<<0);        // set CH1 as source
   } else if (b->GetBoardType() == 7) { // Evaluation Board V3
      b->EnableTrigger(0, 1);           // lemo off, analog trigger on
      b->SetTriggerSource(0);           // use CH1 as source
   }
   b->SetTriggerLevel(0.01);            // 0.05 V
   b->SetTriggerPolarity(false);        // positive edge
   
   /* use following lines to set individual trigger elvels */
   //b->SetIndividualTriggerLevel(1, 0.1);
   //b->SetIndividualTriggerLevel(2, 0.2);
   //b->SetIndividualTriggerLevel(3, 0.3);
   //b->SetIndividualTriggerLevel(4, 0.4);
   //b->SetTriggerSource(15);
   
   b->SetTriggerDelayNs(0);             // zero ns trigger delay
   
   /* use following lines to enable the external trigger */
   //if (b->GetBoardType() == 8) {     // Evaluaiton Board V4
   //   b->EnableTrigger(1, 0);           // enable hardware trigger
   //   b->SetTriggerSource(1<<4);        // set external trigger as source
   //} else {                          // Evaluation Board V3
   //   b->EnableTrigger(1, 0);           // lemo on, analog trigger off
   // }

   /* open file to save waveforms */
   //f = fopen("data.txt", "w");
   //if (f == NULL) {
   //   perror("ERROR: Cannot open file \"data.txt\"");
   //   return 1;
   //}
   
   /* repeat ten times */
   for (j=0 ; j<30000 ; j++) {

      if (j%1000 == 0)
      {
         printf("Event %d read.\n", j);
      }

      /* start board (activate domino wave) */
      b->StartDomino();

      /* wait for trigger */
      //printf("Waiting for trigger...");
      
      fflush(stdout);
      while (b->IsBusy());

      /* read all waveforms */
      b->TransferWaves(0, 8);

      /* read time (X) array of first channel in ns */
      b->GetTime(0, 0, b->GetTriggerCell(0), time_array[0]);

      /* decode waveform (Y) array of first channel in mV */
      b->GetWave(0, 0, wave_array[0]);

      /* read time (X) array of second channel in ns
       Note: On the evaluation board input #1 is connected to channel 0 and 1 of
       the DRS chip, input #2 is connected to channel 2 and 3 and so on. So to
       get the input #2 we have to read DRS channel #2, not #1. */
      b->GetTime(0, 2, b->GetTriggerCell(0), time_array[1]);

      /* decode waveform (Y) array of second channel in mV */
      b->GetWave(0, 2, wave_array[1]);

      // third channel
      b->GetTime(0, 4, b->GetTriggerCell(0), time_array[2]);
      b->GetWave(0, 4, wave_array[2]);

      // fourth channel
      b->GetTime(0, 6, b->GetTriggerCell(0), time_array[3]);
      b->GetWave(0, 6, wave_array[3]);

      int spike;
      spike = RemoveSpikes(false, b->GetTriggerCell(0), wave_array);

      if (spike)
      {
         --j;
         continue;
      }

      /* Save waveform: X=time_array[i], Yn=wave_array[n][i] */
      //printf("Event #%d ----------------------\n  t1[ns]  u1[mV]  t2[ns] u2[mV]\n", j);
      //for (i=0 ; i<1024 ; i++)
      //   fprintf(f, "%7.3f %7.1f %7.3f %7.1f\n", time_array[0][i], wave_array[0][i], time_array[1][i], wave_array[1][i]);

      tree_drs->Fill();

      /* print some progress indication */
      //printf("\rEvent #%d read successfully\n", j);
   }

   tree_drs->Write();
   tree_drs->Print();
   file_drs->Close();
   delete file_drs;
   //fclose(f);
   
   /* delete DRS object -> close USB connection */
   delete drs;
}
