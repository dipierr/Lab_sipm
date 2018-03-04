# Ana_Traces_SiPM_ReadMe.md


## Devices


* Traces acquired by AGILENT MS06054A or by Digitizer CAEN DT 5751 (set the variables in OPTIONS section)
  * Files from Digitizer CAEN DT 5751 are like:
  >Record Length: 1029  
  >BoardID: 31  
  >Channel: 0  
  >Event Number: 0  
  >Pattern: 0x0000  
  >Trigger Time Stamp: 149244  
  >DC offset (DAC): 0x2E14  
  >800  
  >[...]  
  >Record Length: 1029  
  >[...]  

* change the code (READ section) to open other files

---

## Calculations


* DCR            (Hamamatsu - MPPC Characterization)
* Cross Talk     (Hamamatsu - MPPC Characterization pag 44)
* GAIN           (Hamamatsu - MPPC Characterization pag 36)

---

## Tecniques

* DLED for peak detections

---
Set the correct variables (OPTIONS section) to do what you need

---

### HOW TO COMPILE:
>$ root -l  
>root[0] .L Ana_Traces_SiPM.C++

### HOW TO USE
  * Analyze a SiPM at a certain HV and view some useful plots:  
  > root[1] Ana1("file1", last_event_n, display_one_ev_param, LED_bool);  

  * Dark Count Rate and Cross Talk for a SiPM at a certain HV:
  > root[1] DCR_CT_1SiPM_1HVs("file1", last_event_n)
    * file1 (string), see above 
    * last_event_n (int), see above

  * Dark Count Rate and Cross Talk for a SiPM at HV = 34 V, 35 V and 36 V:
  > root[1] DCR_CT_1SiPM_3HVs("file1", "file2", "file3", last_event_n)
    * file1 (string), see above
    * file2 (string), see above
    * file3 (string), see above
    * last_event_n (int), see above
    
  * WHERE:
   *  file1 (string) is the name of the file you want to analyze (e.g.: if you want to analyze a file named wave0.txt, just type "wave0.txt")
   * file2 (string), as file1
   * file3 (string), as file1
   * last_event_n (int) is the last event you want to analyze (e.g.: if you want to analyze the first 1000 events stored in the file "file1", just type 1000)
   * display_one_ev_param (bool)
      *  true if you want to display a single event (trace from file and trace after the DLED procedure)
      *  false if you want to see all the traces (oscilloscope mode)
   *   LED_bool (bool)
       *  true in case of a LED measure (remember to set the window in which data will be analyzed)
       *  false in Dark Mode
---

See the code for more information.
