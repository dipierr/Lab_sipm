# Ana_Traces_SiPM_ReadMe.txt


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
(remember to set in the correct way the variables in OPTIONS section):

* To analyze a single file
> root[0] Analysis(string file1,int last_event_n,bool display)  

* To loop on a single file (useful for plots, e.g. DCR)
> root[0] loopAnalysis(string file1, int last_event_n)

* To analyze 3 files  
> root[0] Ana3(string file1, string file2, string file3, int last_event_n)

* To loop on 3 files: (useful for DCR plots)
>root[0] loopAna3(string file1, string file2, string file3, int last_event_n)
