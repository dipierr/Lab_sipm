# Lab_sipm
Tools for analysis of data taken in SiPM-4-Astroparticle Lab.

In this repository we want to put the codes for the analysis of the signals acquired at the Lab.
We are currently using the CAEN DT 5751 Digitizer (with its acquisition sw: "wavedump").
The signals come from different SiPMs and amplifiers under test.
These sw tools are based on ROOT (cern), which then is a prerequisite.

The tools are:
1. a ROOT macro for displaying the traces (display_traces_1ch.C)
2. a skeleton cpp code for reading the traces and extracting information (sipm_spectrum.cxx). 

Usage of display_traces_1ch.C
root$ .L display_traces_1ch.C
root$ display_traces_1ch("wave0_ExtTrg-Pd1.txt")
where "wave0_ExtTrg-Pd1.txt" is the file with input data

Usage of sipm_spectrum.cxx
compile:$ g++ sipm_spectrum.cxx -o sipm_spectrum.exe `root-config --cflags --libs` 
$./sipm_spectrum.exe sipm_spectrum.cfg
