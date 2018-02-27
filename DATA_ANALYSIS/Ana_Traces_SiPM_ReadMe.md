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
1. PREDEFINED FUNCTIONS
  * Dark Count Rate and Cross Talk for a SiPM at a certain HV:
  > root[1] DCR_CT_1SiPM_1HVs("file1", last_event_n)

  * Dark Count Rate and Cross Talk for a SiPM at HV = 34 V, 35 V and 36 V:
  > root[1] DCR_CT_1SiPM_3HVs("file1", "file2", "file3", last_event_n)

  * Gain for a SiPM at a certain HV:
  > root[1] GAIN_1SiPM_1HV("file1", last_event_n)

  * Analyze a SiPM at a certain HV and view some useful plots:
  > root[1] Ana1("file1", last_event_n)

2. CUSTOM FUNCTIONS
  (remember to set in the correct way the variables in OPTIONS section):

  * To analyze a single file
  > root[1] Analysis("file1", last_event_n, display)  

  * To loop on a single file (useful for plots, e.g. DCR)
  > root[1] loopAnalysis("file1", last_event_n)

  * To analyze 3 files  
  > root[1] Ana3("file1", "file2", "file3", last_event_n)

  * To loop on 3 files: (useful for DCR plots)
  >root[1] loopAna3("file1", "file2", "file3", last_event_n)

NOTE:
* "file1" = name of the file you want to analyze
* last_event_n = last event to analyze (integer)
* display = bool (true or false)

See the code for more information.
