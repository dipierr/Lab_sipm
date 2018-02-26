Ana_Traces_SiPM_ReadMe.txt
==========================

Devices
-------

*Traces acquired by AGILENT MS06054A or by Digitizer CAEN DT 5751 (set the variables in OPTIONS section)
*Files from Digitizer CAEN DT 5751 are like:
    * Record Length: 1029
    * BoardID: 31
    * Channel: 0
    * Event Number: 0
    * Pattern: 0x0000
    * Trigger Time Stamp: 149244
    * DC offset (DAC): 0x2E14
    * 800
    * [...]
    * Record Length: 1029
    * [...]
    
*change the code (READ section) to open other files

---

Calculations
------------

*DCR            (Hamamatsu - MPPC Characterization)
*Cross Talk     (Hamamatsu - MPPC Characterization pag 44)

---

##Tecniques
-----------
*DLED for peak detections

---
Set the correct variables (OPTIONS section) to do what you need
---

##HOW TO COMPILE:
>$ root -l
>root[0] .L Ana_Traces_SiPM.C++

###EXAMPLE (remember to set in the correct way the variables in OPTIONS section):

*To analyze a single file
>root[0] Analysis("20180221_HD3-2_3_DARK_36_AS_2_01.txt",15000,false)
